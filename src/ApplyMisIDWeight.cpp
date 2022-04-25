// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Mon Apr 25, 2022 at 05:23 PM -0400
//
// Description: unfolding weights applyer (A)

#include <iostream>
#include <regex>
#include <string>
#include <tuple>
#include <vector>

#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <TFile.h>
#include <TH3D.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandomGen.h>
#include <TString.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>

#include <yaml-cpp/yaml.h>
#include <boost/range/join.hpp>
#include <cxxopts.hpp>

#include "kinematic.h"
#include "utils.h"

#define RAND_SEED 42

using namespace std;
using ROOT::RDataFrame;
using ROOT::Math::PxPyPzEVector;
using ROOT::Math::XYZVector;
using ROOT::RDF::RNode;

///////////////////
// Configuration //
///////////////////

typedef vector<pair<string, string>> vPStrStr;
typedef vector<pair<regex, string>>  vPRegStr;

static const vPStrStr MU_BRANCH_DEFS{
    // simple name, complex name
    {"MC15TuneV1_ProbNNpi", "MC15TuneV1_ProbNNpi"},
    {"MC15TuneV1_ProbNNk", "MC15TuneV1_ProbNNk"},
    {"MC15TuneV1_ProbNNp", "MC15TuneV1_ProbNNp"},
    {"MC15TuneV1_ProbNNe", "MC15TuneV1_ProbNNe"},
    {"MC15TuneV1_ProbNNmu", "MC15TuneV1_ProbNNmu"},
    {"MC15TuneV1_ProbNNghost", "MC15TuneV1_ProbNNghost"},
    {"DLLK", "PIDK"},
    {"DLLp", "PIDp"},
    {"DLLe", "PIDe"},
    {"DLLmu", "PIDmu"},
    {"DLLd", "PIDd"},
    {"IsMuon", "isMuon"},
    {"P", "P"},
    {"PZ", "PZ"},
    {"InMuonAcc", "InMuonAcc"}};

static const vPRegStr CUT_REPLACE_RULES{{regex("&"), "&&"},
                                        {regex("\\|"), "||"}};

static const string DST_BR_PREFIX = "dst";
static const string D0_BR_PREFIX  = "d0";
static const string B0_BR_PREFIX  = "b0";
static const string B_BR_PREFIX   = "b";
static const string MU_BR_PREFIX  = "mu";

static const string DST_TEST_BR = "dst_PX";
static const string D0_TEST_BR  = "d0_PX";

////////////////////////////////////
// Helpers for weight computation //
////////////////////////////////////

vector<TString> buildHistoWtNames(string targetParticle, YAML::Node node) {
  vector<TString> result{};

  for (auto it = node.begin(); it != node.end(); it++) {
    auto srcPtcl = it->first.as<string>();
    auto name    = srcPtcl + "TagTo" + capitalize(targetParticle) + "Tag";
    result.emplace_back(name);
  }

  return result;
}

vector<TString> buildHistoSmrWtnames(YAML::Node node) {
  vector<TString> result{};
  vector<TString> targetParticles = {"K", "Pi"};

  for (auto it = node.begin(); it != node.end(); it++) {
    for (auto tgt : targetParticles) {
      auto name = it->first.as<string>() + "TagTo" + tgt + "True";
      result.emplace_back(name);
    }
  }

  return result;
}

// Idea stolen from:
//   https://root-forum.cern.ch/t/running-rdataframes-define-in-for-loop/32484/2
RNode defineBranch(RNode df, string particle = "mu",
                   const vPStrStr& rules = MU_BRANCH_DEFS, int idx = 0) {
  // auto df = init_df.Alias(alias, particle+"_"+raw);
  if (rules.size() == idx) return df;

  auto inputBrName = rules[idx].second;
  if (particle != ""s) inputBrName = particle + "_" + inputBrName;

  return defineBranch(df.Define(rules[idx].first, inputBrName), particle, rules,
                      idx + 1);
}

pair<vPStrStr, vector<string>> genCutDirective(YAML::Node    node,
                                               const string& wtPrefix,
                                               string brPrefix = "is_misid_") {
  vPStrStr       directives{};
  vector<string> outputBrs{};
  const auto     wtTargetParticle = "MuTag";
  const auto     wtSmrParticles   = {"k", "pi"};

  vector<string> particles{};
  // first find particles and cuts
  for (auto it = node.begin(); it != node.end(); it++) {
    particles.emplace_back(it->first.as<string>());
    auto cut = it->second.as<string>();
    for (auto reg = CUT_REPLACE_RULES.begin(); reg != CUT_REPLACE_RULES.end();
         reg++) {
      cut = regex_replace(cut, reg->first, reg->second);
    }
    directives.emplace_back(pair{it->first.as<string>(), cut});
  }

  // generate the alias for output cut branches
  for (auto idx = 0; idx != particles.size(); idx++) {
    outputBrs.emplace_back(brPrefix + particles[idx]);
    directives.emplace_back(pair{brPrefix + particles[idx], particles[idx]});
  }

  // generate the automatic weight for each event based on the species of the
  // event
  auto expr  = ""s;
  auto first = true;

  for (const auto& p : particles) {
    auto wtBrName = wtPrefix + "_" + p + "TagTo" + wtTargetParticle;
    if (!first) expr += " + ";
    first = false;
    expr += brPrefix + p + "*" + wtBrName;
  }
  outputBrs.emplace_back(wtPrefix);
  directives.emplace_back(pair{wtPrefix, expr});
  cout << "  " << wtPrefix << " = " << expr << endl;

  // generate the DiF smearing weight for each event
  vector<string> brSmrNames{};
  for (const auto& tgt : wtSmrParticles) {
    expr  = ""s;
    first = true;
    for (const auto& p : particles) {
      auto wtBrName = wtPrefix + "_" + p + "TagTo" + capitalize(tgt) + "True";
      if (!first) expr += " + ";
      first = false;
      expr += brPrefix + p + "*" + wtBrName;
    }

    auto outputBr = wtPrefix + "_" + tgt + "_smr";
    brSmrNames.push_back(outputBr);
    outputBrs.emplace_back(outputBr);
    directives.emplace_back(pair{outputBr, expr});
    cout << "  " << outputBr << " = " << expr << endl;
  }

  // generate the DiF no smearing weight
  auto brNoSmr = wtPrefix + "_no_smr";
  outputBrs.emplace_back(wtPrefix + "_no_smr");

  expr = "1.0"s;
  for (const auto& smr : brSmrNames) {
    expr += " - " + smr;
  }
  cout << "  " << brNoSmr << " = " << expr << endl;
  directives.emplace_back(pair{brNoSmr, expr});

  return {directives, outputBrs};
}

/////////////////////////////////
// Rest frame variable helpers //
/////////////////////////////////

auto getRandSmrHelper(vector<vector<double>>& smr) {
  auto size = make_shared<unsigned long>(smr.size());
  auto rng  = make_shared<TRandomMixMax256>(42);

  return [&smr, size, rng] {
    unsigned long rand = rng->Uniform(0, *(size.get()));
    return smr[rand];
  };
}

vector<string> setBrPrefix(const string prefix, const vector<string> vars) {
  vector<string> result{};
  for (const auto& v : vars) result.emplace_back(prefix + "_" + v);
  return result;
}

void getSmrFac(vector<vector<double>>& result, string auxFile,
               string prefix = "k") {
  auto df = RDataFrame(prefix + "_smr", auxFile);
  df.Foreach(
      [&](double x, double y, double z) {
        result.emplace_back(vector<double>{x, y, z});
      },
      setBrPrefix(prefix + "_smr", {"x", "y", "z"}));
}

template <typename F>
RNode computeFitVars(RNode df, F& randGetter, double mMuHypo, double mB,
                     string suffix, vector<string>& outputBrs) {
  // we probably did some unnecessary copies here, but deducing those nested
  // lambdas can be quite hard so I'm just being lazy here.
  auto rebuildMu4MomHypo = [=, &randGetter](PxPyPzEVector v4Mu) {
    vector<double> smr = randGetter();
    return rebuildMu4Mom(v4Mu, smr, mMuHypo);
  };
  auto estB4MomHypo = [=](PxPyPzEVector v4BReco, XYZVector v3BFlight) {
    return estB4Mom(v4BReco, v3BFlight, mB);
  };

  outputBrs.emplace_back("mm2" + suffix);
  outputBrs.emplace_back("q2" + suffix);
  outputBrs.emplace_back("el" + suffix);

  // FIXME: Not sure if 'el' is defined this way, as we probably don't have a
  // lepton at all.
  return df.Define("v4_mu" + suffix, rebuildMu4MomHypo, {"v4_mu"})
      .Define("v4_b_reco" + suffix, "v4_mu" + suffix + " + v4_d")
      .Define("v4_b_est" + suffix, estB4MomHypo,
              {"v4_b_reco" + suffix, "v3_b_dir"})
      .Define("mm2" + suffix, m2Miss,
              {"v4_b_est" + suffix, "v4_b_reco" + suffix})
      .Define("q2" + suffix, q2, {"v4_b_est" + suffix, "v4_d"})
      .Define("el" + suffix, el, {"v4_b_est" + suffix, "v4_mu" + suffix});
}

template <typename F1, typename F2>
pair<RNode, vector<string>> defRestFrameVars(RNode df, TTree* tree,
                                             F1& randKGetter,
                                             F2& randPiGetter) {
  vector<string> outputBrs{};
  string         dMeson = ""s;
  string         bMeson = ""s;
  double         mB;

  if (branchExists(tree, DST_TEST_BR)) {
    dMeson = DST_BR_PREFIX;
    bMeson = B0_BR_PREFIX;
    mB     = B0_M;
  } else if (branchExists(tree, D0_TEST_BR)) {
    dMeson = D0_BR_PREFIX;
    bMeson = B_BR_PREFIX;
    mB     = B_M;
  } else {
    cout << "No known branch found for D0 nor D*. Exit now..." << endl;
    exit(1);
  }

  // Basic vectors that we use and not going to change
  df = df.Define(
             "v4_d",
             [](double px, double py, double pz, double e) {
               return PxPyPzEVector(px, py, pz, e);
             },
             setBrPrefix(dMeson, {"PX", "PY", "PZ", "PE"}))
           .Define(
               "v4_mu",
               [](double px, double py, double pz, double e) {
                 return PxPyPzEVector(px, py, pz, e);
               },
               setBrPrefix(MU_BR_PREFIX, {"PX", "PY", "PZ", "PE"}))
           .Define("v3_b_dir", buildBFlightDir,
                   setBrPrefix(bMeson, {"ENDVERTEX_X", "OWNPV_X", "ENDVERTEX_Y",
                                        "OWNPV_Y", "ENDVERTEX_Z", "OWNPV_Z"}));

  // Replace mass hypo and compute fit vars
  df = computeFitVars(df, randPiGetter, PI_M, mB, "_smr_pi", outputBrs);
  df = computeFitVars(df, randKGetter, K_M, mB, "_smr_k", outputBrs);

  return {df, outputBrs};
}

//////////
// Main //
//////////

int main(int argc, char** argv) {
  cxxopts::Options argOpts("ApplyMisIDWeight",
                           "unfolding weihgts applyer (A).");

  // clang-format off
  argOpts.add_options()
    // general
    ("h,help", "print help")
    ("Y,year", "sample year",
     cxxopts::value<string>()->default_value("2016"))
    // I/O
    ("i,input", "specify input ntuple", cxxopts::value<string>())
    ("x,aux", "specify auxiliary ntuple", cxxopts::value<string>())
    ("o,output", "specify output ntuple", cxxopts::value<string>())
    ("c,config", "specify input YAML config file",
     cxxopts::value<string>())
    // flags (typically don't change these)
    ("a,alias", "apply aliases.",
     cxxopts::value<bool>()->default_value("false"))
    ("p,particle", "specify alias particle",
     cxxopts::value<string>()->default_value("mu"))
  ;
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  // get options
  auto ntpIn      = parsedArgs["input"].as<string>();
  auto ntpOut     = parsedArgs["output"].as<string>();
  auto ntpAux     = parsedArgs["aux"].as<string>();
  auto particle   = parsedArgs["particle"].as<string>();
  auto applyAlias = parsedArgs["alias"].as<bool>();

  // parse YAML config
  auto ymlFile    = parsedArgs["config"].as<string>();
  auto ymlConfig  = YAML::LoadFile(ymlFile);
  auto year       = parsedArgs["year"].as<string>();
  auto ymlDirPath = absDirPath(parsedArgs["config"].as<string>());
  auto weightBrs  = ymlConfig["weight_brs"][year];
  auto filePrefix = absDirPath(ymlFile);

  TFile*        ntpHisto;
  vector<TH3D*> histos;
  // snapshot option
  auto writeOpts  = ROOT::RDF::RSnapshotOptions{};
  writeOpts.fMode = "UPDATE";
  bool firstTree  = true;
  for (auto entry : weightBrs) {
    auto histoPrefix    = TString(entry["prefix"].as<string>());
    auto histoFile      = entry["file"].as<string>();
    auto treeName       = entry["tree"].as<string>();
    auto weightBrPrefix = entry["name"].as<string>();
    histoFile           = filePrefix + "/" + histoFile;

    auto ntpInTest = new TFile(TString(ntpIn));
    auto treeTest  = dynamic_cast<TTree*>(ntpInTest->Get(TString(treeName)));
    if (treeTest == nullptr) {
      cout << treeName << " doesn't exist in " << ntpIn << " skipping..."
           << endl;
      continue;
    }

    cout << "Handling tree " << treeName << " from histos of prefix "
         << histoPrefix << " from file " << histoFile << endl;

    ntpHisto              = new TFile(TString(histoFile), "READ");
    auto           dfInit = RDataFrame(treeName, ntpIn);
    RNode          df     = static_cast<RNode>(dfInit);
    vector<string> outputBrNames{"runNumber", "eventNumber"};
    vector<string> weightTagBrNames{};

    if (applyAlias) {
      df = defineBranch(dfInit, particle);
      // compute ETA
      df = df.Define("ETA",
                     [](double& p, double& pz) {
                       return 0.5 * TMath::Log((p + pz) / (p - pz));
                     },
                     {"P", "PZ"});
    }

    auto histoWtNames    = buildHistoWtNames(particle, ymlConfig["tags"]);
    auto histoSmrWtNames = buildHistoSmrWtnames(ymlConfig["tags"]);
    cout << "Generate transfer factors/DiF smearing wieghts for all species"
         << endl;
    for (const auto h : boost::join(histoWtNames, histoSmrWtNames)) {
      auto histoName = histoPrefix + "__" + h;
      auto histoWt   = static_cast<TH3D*>(ntpHisto->Get(histoName));
      histos.emplace_back(histoWt);
      cout << "  Loading histo " << histoName << endl;

      auto brName = weightBrPrefix + "_" + string(h);
      weightTagBrNames.emplace_back(brName);
      cout << "  Generating " << brName << "..." << endl;
      df = df.Define(brName,
                     [histoWt](double& x, double& y, double& z) {
                       auto binIdx = histoWt->FindFixBin(x, y, z);
                       return histoWt->GetBinContent(binIdx);
                     },
                     {"P", "ETA", "nTracks"});
      outputBrNames.emplace_back(brName);
    }

    // apply tagged species cuts
    auto [directives, addOutputBrs] =
        genCutDirective(ymlConfig["tags"], weightBrPrefix);
    df = defineBranch(df, ""s, directives);
    for (auto br : addOutputBrs) outputBrNames.emplace_back(br);

    // read smearing factors from aux ntuple
    auto smrFacK  = vector<vector<double>>{};
    auto smrFacPi = vector<vector<double>>{};
    getSmrFac(smrFacK, ntpAux);
    getSmrFac(smrFacPi, ntpAux, "pi");

    // we can call these functions directly to get random smearing factors
    auto randSmrFacK  = getRandSmrHelper(smrFacK);
    auto randSmrFacPi = getRandSmrHelper(smrFacPi);

    // Recompute fit vars
    auto [dfOut, moreOutputBrs] =
        defRestFrameVars(df, treeTest, randSmrFacK, randSmrFacPi);
    for (auto br : moreOutputBrs) outputBrNames.emplace_back(br);

    cout << "Writing to " << ntpOut << endl;
    if (firstTree) {
      dfOut.Snapshot(treeName, ntpOut, outputBrNames);
      firstTree = false;
    } else
      dfOut.Snapshot(treeName, ntpOut, outputBrNames, writeOpts);

    // cleanups
    for (auto& h : histos) delete h;
    delete ntpHisto;

    delete treeTest;
    delete ntpInTest;
  }
}
