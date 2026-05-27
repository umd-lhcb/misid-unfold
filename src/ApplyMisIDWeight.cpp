// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Sun Sep 25, 2022 at 09:43 PM -0400
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
#include <cxxopts.hpp>

#include "kinematic.h"
#include "utils.h"

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

constexpr int    RAND_SEED            = 42;
constexpr double PRE_SCALE_CORRECTION = 10.0;

static const vPStrStr MU_BRANCH_DEFS{
    // PIDCalib name, DaVinci name w/o particle name
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
    {"InMuonAcc", "InMuonAcc"},
    ////
    {"Brunel_MC15TuneV1_ProbNNpi", "MC15TuneV1_ProbNNpi"},
    {"Brunel_MC15TuneV1_ProbNNk", "MC15TuneV1_ProbNNk"},
    {"Brunel_MC15TuneV1_ProbNNp", "MC15TuneV1_ProbNNp"},
    {"Brunel_MC15TuneV1_ProbNNe", "MC15TuneV1_ProbNNe"},
    {"Brunel_MC15TuneV1_ProbNNmu", "MC15TuneV1_ProbNNmu"},
    {"Brunel_MC15TuneV1_ProbNNghost", "MC15TuneV1_ProbNNghost"},
    {"Brunel_DLLK", "PIDK"},
    {"Brunel_DLLp", "PIDp"},
    {"Brunel_DLLe", "PIDe"},
    {"Brunel_DLLmu", "PIDmu"},
    {"Brunel_DLLd", "PIDd"},
    {"Brunel_IsMuon", "isMuon"},
    {"Brunel_InMuonAcc", "InMuonAcc"},
    {"Brunel_P", "P"},
    // special aliases for compute eta
    {"P", "P"},
    {"PZ", "PZ"},
};

static const vPRegStr CUT_REPLACE_RULES{{regex("&"), "&&"},
                                        {regex("\\|"), "||"}};

static const string DST_BR_PREFIX = "dst";
static const string D0_BR_PREFIX  = "d0";
static const string B0_BR_PREFIX  = "b0";
static const string B_BR_PREFIX   = "b";
static const string MU_BR_PREFIX  = "mu";

static const string DST_TEST_BR = "dst_PX";
static const string D0_TEST_BR  = "d0_PX";

////////////////////////
// RDataFrame helpers //
////////////////////////

// Idea stolen from:
//   https://root-forum.cern.ch/t/running-rdataframes-define-in-for-loop/32484/2
RNode defineBranch(RNode df, const string& particle = "mu",
                   const vPStrStr& rules = MU_BRANCH_DEFS,
                   const bool& debug = false, unsigned idx = 0) {
  if (rules.size() == idx) return df;

  auto inputBrName = rules[idx].second;
  if (particle != ""s) inputBrName = particle + "_" + inputBrName;
  if (debug)
    cout << "Define " << rules[idx].first << " as " << inputBrName << endl;

  return defineBranch(df.Define(rules[idx].first, inputBrName), particle, rules,
                      debug, idx + 1);
}

///////////////////////////////////////////
// Helpers for apply tagged species cuts //
///////////////////////////////////////////

pair<vPStrStr, vector<string>> genTaggedCutDirective(
    YAML::Node node, const string& brPrefix = "is_misid_") {
  vPStrStr       directives{};
  vector<string> outputBrs{};

  vector<string> particles{};
  // find particles and cuts
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
  for (unsigned idx = 0; idx != particles.size(); idx++) {
    outputBrs.emplace_back(brPrefix + particles[idx]);
    directives.emplace_back(pair{brPrefix + particles[idx], particles[idx]});
  }

  return {directives, outputBrs};
}

////////////////////////////////////
// Helpers for weight computation //
////////////////////////////////////

vector<TString> buildHistoWtNames(const string&          targetParticle,
                                  const vector<TString>& skims,
                                  const string& year, YAML::Node node) {
  vector<TString> result{};
  for (const auto& skim : skims) {
    for (auto it = node.begin(); it != node.end(); it++) {
      auto srcPtcl = it->first.as<string>();
      auto name    = srcPtcl + "TagTo" + capitalize(targetParticle) + "Tag_" +
                     year.substr(2, 2) + "_" + skim;
      result.emplace_back(name);
      for (auto itTrue = node.begin(); itTrue != node.end(); itTrue++) {
        auto srcPtclSingle = itTrue->first.as<string>();
        auto nameSingle    = name + "_" + srcPtclSingle + "TrueOnly";
        result.emplace_back(nameSingle);
      }
    }
  }
  return result;
}

tuple<RNode, vector<string>, vector<TH3D*>> applyWtFromHistos(
    RNode df, TFile* ntpHisto, const string& histoPrefix,
    const string& weightBrPrefix, const vector<TString>& histoWtNames,
    const bool& debug = false) {
  auto outputBrs = vector<string>{};
  auto histos    = vector<TH3D*>{};

  for (const auto& h : histoWtNames) {
    const TString histoName = histoPrefix + "__" + h;
    if (debug) cout << "  Loading histo " << histoName << endl;
    auto histoWt = static_cast<TH3D*>(ntpHisto->Get(histoName));
    if (!histoWt)
      throw runtime_error("Could not get histogram with name " + histoName);
    histos.emplace_back(histoWt);

    double prescale = 1.0;
    if (TString(h).Contains("MuTag")) {
      if (debug) {
        cout << "  " << h
             << " is a transfer-factor weight, apply 10x enhance due to "
                "prescale."
             << endl;
      }
      prescale = PRE_SCALE_CORRECTION;
    }

    // Relies on the fact that year comes before skim
    const auto    year_idx = h.First("1");
    const TString h_noyear = TString(h).Replace(year_idx - 1, 3, 0);
    auto          brName   = weightBrPrefix + "_" + h_noyear;
    if (debug) cout << "  Generating " << brName << "..." << endl;
    df = df.Define(brName,
                   [histoWt, prescale](double& x, double& y, double& z) {
                     auto binIdx = histoWt->FindFixBin(x, y, z);
                     return histoWt->GetBinContent(binIdx) * prescale;
                   },
                   {"P", "ETA", "nTracks"});
    outputBrs.emplace_back(brName);

    // Also produce weights with smeared muon kinematics
    if (histoName.EndsWith("_kTrueOnly")) {
      const auto brNameSmr = brName + "_smr_k";
      if (histoName.Contains("_vmu_")) {
        if (debug)
          cout << "  Generating " << brNameSmr << " with vmu smearing..."
               << endl;
        df = df.Define(brNameSmr,
                       [histoWt, prescale](double& x, double& y, double& z) {
                         auto binIdx = histoWt->FindFixBin(x, y, z);
                         return histoWt->GetBinContent(binIdx) * prescale;
                       },
                       {"mu_p_smr_k_vmu", "mu_eta_smr_k_vmu", "nTracks"});
      } else {
        if (debug) cout << "  Generating " << brNameSmr << "..." << endl;
        df = df.Define(brNameSmr,
                       [histoWt, prescale](double& x, double& y, double& z) {
                         auto binIdx = histoWt->FindFixBin(x, y, z);
                         return histoWt->GetBinContent(binIdx) * prescale;
                       },
                       {"mu_p_smr_k", "mu_eta_smr_k", "nTracks"});
      }
      outputBrs.emplace_back(brNameSmr);
    } else if (histoName.EndsWith("_piTrueOnly")) {
      const auto brNameSmr = brName + "_smr_pi";
      if (histoName.Contains("_vmu_")) {
        if (debug)
          cout << "  Generating " << brNameSmr << " with vmu smearing..."
               << endl;
        df = df.Define(brNameSmr,
                       [histoWt, prescale](double& x, double& y, double& z) {
                         auto binIdx = histoWt->FindFixBin(x, y, z);
                         return histoWt->GetBinContent(binIdx) * prescale;
                       },
                       {"mu_p_smr_pi_vmu", "mu_eta_smr_pi_vmu", "nTracks"});
      } else {
        if (debug) cout << "  Generating " << brNameSmr << "..." << endl;
        df = df.Define(brNameSmr,
                       [histoWt, prescale](double& x, double& y, double& z) {
                         auto binIdx = histoWt->FindFixBin(x, y, z);
                         return histoWt->GetBinContent(binIdx) * prescale;
                       },
                       {"mu_p_smr_pi", "mu_eta_smr_pi", "nTracks"});
      }
      outputBrs.emplace_back(brNameSmr);
    }
  }

  return {df, outputBrs, histos};
}

pair<vPStrStr, vector<string>> genWtDirective(
    YAML::Node node, const string& wtPrefix, const vector<TString>& skims,
    const bool& debug = false, const string& brPrefix = "is_misid_") {
  vPStrStr       directives{};
  vector<string> outputBrs{};
  const auto     wtTargetParticle = "MuTag";

  vector<string> particles{};
  // first find particles
  for (auto it = node.begin(); it != node.end(); it++)
    particles.emplace_back(it->first.as<string>());

  for (const auto& skim : skims) {
    // generate the automatic weight for each event based on the species of the
    // event
    auto expr  = ""s;
    auto first = true;

    for (const auto& p : particles) {
      auto wtBrName =
          wtPrefix + "_" + p + "TagTo" + wtTargetParticle + "_" + skim;
      if (!first) expr += " + ";
      first = false;
      expr += brPrefix + p + "*" + wtBrName;
    }
    auto wtName = wtPrefix + "_" + skim;
    outputBrs.emplace_back(wtName);
    directives.emplace_back(pair{wtName, expr});
    if (debug) cout << "  " << wtName << " = " << expr << endl;

    // Also generate weights assuming single true type
    for (const auto& pTrue : particles) {
      expr  = ""s;
      first = true;
      for (const auto& p : particles) {
        auto wtBrName = wtPrefix + "_" + p + "TagTo" + wtTargetParticle + "_" +
                        skim + "_" + pTrue + "TrueOnly";
        if (!first) expr += " + ";
        first = false;
        expr += brPrefix + p + "*" + wtBrName;
      }
      auto wtNameSingleTrue = wtName + "_" + pTrue;
      outputBrs.emplace_back(wtNameSingleTrue);
      directives.emplace_back(pair{wtNameSingleTrue, expr});
      if (debug) cout << "  " << wtNameSingleTrue << " = " << expr << endl;

      // Also produce weights with smeared muon kinematics
      if (wtNameSingleTrue.EndsWith("_k") || wtNameSingleTrue.EndsWith("_pi")) {
        const TString smr_suffix =
            wtNameSingleTrue.EndsWith("_k") ? "_smr_k" : "_smr_pi";
        expr  = ""s;
        first = true;
        for (const auto& p : particles) {
          const auto wtBrName = wtPrefix + "_" + p + "TagTo" +
                                wtTargetParticle + "_" + skim + "_" + pTrue +
                                "TrueOnly" + smr_suffix;
          if (!first) expr += " + ";
          first = false;
          expr += brPrefix + p + "*" + wtBrName;
        }
        auto wtNameSingleTrueSmr = wtNameSingleTrue + smr_suffix;
        outputBrs.emplace_back(wtNameSingleTrueSmr);
        directives.emplace_back(pair{wtNameSingleTrueSmr, expr});
        if (debug) cout << "  " << wtNameSingleTrueSmr << " = " << expr << endl;
      }
    }
  }

  return {directives, outputBrs};
}

/////////////////////////////////
// Rest frame variable helpers //
/////////////////////////////////

auto getRandSmrHelper(const vector<vector<double>>& smr) {
  auto size = make_shared<unsigned long>(smr.size());
  auto rng  = make_shared<TRandomMixMax256>(RAND_SEED);

  return [&smr, size, rng] {
    unsigned long rand = rng->Uniform(0, *(size.get()));
    return smr[rand];
  };
}

vector<string> setBrPrefix(const string& prefix, const vector<string>& vars) {
  vector<string> result{};
  result.reserve(vars.size());
  for (const auto& v : vars) result.emplace_back(prefix + "_" + v);
  return result;
}

void getSmrFac(vector<vector<double>>& result, const string& auxFile,
               const string& prefix = "k_smr") {
  auto df = RDataFrame(prefix, auxFile);
  df.Foreach(
      [&](double x, double y, double z, double rP, double dTheta, double dPhi) {
        result.emplace_back(vector<double>{x, y, z, rP, dTheta, dPhi});
      },
      setBrPrefix(prefix, {"x", "y", "z", "rP", "dTheta", "dPhi"}));
}

template <typename F>
RNode computeDiFVars(RNode df, F& randGetter, const string& suffix,
                     vector<string>& outputBrs, const string& smr_mode) {
  // we probably did some unnecessary copies here, but deducing those nested
  // lambdas can be quite hard so I'm just being lazy here.
  auto rebuildMu4MomPartial = [=, &randGetter](PxPyPzEVector v4Mu) {
    vector<double> smr = randGetter();
    return rebuildMu4Mom(v4Mu, smr, smr_mode);
  };
  auto estB4MomPartial = [=](PxPyPzEVector v4BReco, XYZVector v3BFlight) {
    return estB4Mom(v4BReco, v3BFlight);
  };

  vector<string> brNames = {"mm2",   "q2",     "el",  "b_m",
                            "v4_mu", "mu_eta", "mu_p"};
  for (auto& n : brNames) outputBrs.emplace_back(n + suffix);

  return df.Define("v4_mu" + suffix, rebuildMu4MomPartial, {"v4_mu"})
      .Define("v4_b_reco" + suffix, "v4_mu" + suffix + " + v4_d")
      .Define("v4_b_est" + suffix, estB4MomPartial,
              {"v4_b_reco" + suffix, "v3_b_dir"})
      .Define("mm2" + suffix, m2Miss,
              {"v4_b_est" + suffix, "v4_b_reco" + suffix})
      .Define("q2" + suffix, q2, {"v4_b_est" + suffix, "v4_d"})
      .Define("el" + suffix, el, {"v4_b_est" + suffix, "v4_mu" + suffix})
      .Define("b_m" + suffix, calcBM, {"v4_b_reco" + suffix})
      .Define("mu_eta" + suffix, eta, {"v4_mu" + suffix})
      .Define("mu_p" + suffix, P, {"v4_mu" + suffix});
}

template <typename F1, typename F2>
pair<RNode, vector<string>> defRestFrameVars(RNode df, TTree* tree,
                                             F1& randKGetter, F2& randPiGetter,
                                             const string& smr_mode,
                                             const string& brSuffix = "") {
  vector<string> outputBrs{};
  string         dMeson = ""s;
  string         bMeson = ""s;

  if (branchExists(tree, DST_TEST_BR)) {
    dMeson = DST_BR_PREFIX;
    bMeson = B0_BR_PREFIX;
  } else if (branchExists(tree, D0_TEST_BR)) {
    dMeson = D0_BR_PREFIX;
    bMeson = B_BR_PREFIX;
  } else {
    cout << "No known branch found for D0 nor D*. Exit now..." << endl;
    exit(1);
  }

  // Basic vectors that we use and not going to change
  df = df.Define(
             "v4_d" + brSuffix,
             [](double px, double py, double pz, double e) {
               return PxPyPzEVector(px, py, pz, e);
             },
             setBrPrefix(dMeson, {"PX", "PY", "PZ", "PE"}))
           .Define(
               "v4_mu" + brSuffix,
               [](double px, double py, double pz, double e) {
                 return PxPyPzEVector(px, py, pz, e);
               },
               setBrPrefix(MU_BR_PREFIX, {"PX", "PY", "PZ", "PE"}))
           .Define("v3_b_dir" + brSuffix, buildBFlightDir,
                   setBrPrefix(bMeson, {"ENDVERTEX_X", "OWNPV_X", "ENDVERTEX_Y",
                                        "OWNPV_Y", "ENDVERTEX_Z", "OWNPV_Z"}));

  // Replace mass hypo and compute fit vars
  df = computeDiFVars(df, randPiGetter, "_smr_pi" + brSuffix, outputBrs,
                      smr_mode);
  df =
      computeDiFVars(df, randKGetter, "_smr_k" + brSuffix, outputBrs, smr_mode);

  return {df, outputBrs};
}

//////////
// Main //
//////////

int main(int argc, char** argv) {
  cxxopts::Options argOpts("ApplyMisIDWeight",
                           "Unfolding weights applier (A).");

  // clang-format off
  argOpts.add_options()
    // general
    ("h,help", "print help")
    ("d,debug", "enable debug mode",
     cxxopts::value<bool>()->default_value("false"))
    ("Y,year", "sample year",
     cxxopts::value<string>()->default_value("2016"))
    // I/O
    ("i,input", "specify input ntuple", cxxopts::value<string>())
    ("x,aux", "specify auxiliary ntuple", cxxopts::value<string>())
    ("o,output", "specify output ntuple", cxxopts::value<string>())
    ("c,config", "specify input YAML config file",
     cxxopts::value<string>())
    // flags (typically don't change these)
    ("a,alias", "apply aliases.")
    ("p,particle", "specify alias particle",
     cxxopts::value<string>()->default_value("mu"))
    ("kSmrBrName", "specify K-smear branch name",
     cxxopts::value<string>()->default_value("k_smr"))
    ("piSmrBrName", "specify pi-smear branch name",
     cxxopts::value<string>()->default_value("pi_smr"))
    ("smrMode", "specify strategy for misid smearing (options are PThetaPhi or PxPyPz)",
     cxxopts::value<string>()->default_value("PThetaPhi"))
  ;
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  // get options
  const auto ntpIn       = parsedArgs["input"].as<string>();
  const auto ntpOut      = parsedArgs["output"].as<string>();
  const auto ntpAux      = parsedArgs["aux"].as<string>();
  const auto particle    = parsedArgs["particle"].as<string>();
  const auto applyAlias  = parsedArgs["alias"].as<bool>();
  const auto kSmrBrName  = parsedArgs["kSmrBrName"].as<string>();
  const auto piSmrBrName = parsedArgs["piSmrBrName"].as<string>();
  const auto smr_mode    = parsedArgs["smrMode"].as<string>();
  const auto ymlFile     = parsedArgs["config"].as<string>();
  const auto debug       = parsedArgs["debug"].as<bool>();

  // Get correct smearing for misid validation sample
  const string kSmrBrNameVmu  = kSmrBrName + "_ubdt_veto";
  const string piSmrBrNameVmu = piSmrBrName + "_ubdt_veto";

  // Get yml file
  const string ymlName = fileNameFromPath(ymlFile);

  // Dump some information
  cout << "\nApplyMisIDWeight configuration:\n";
  cout << boolalpha;
  cout << "- Input ntuple: " << ntpIn << "\n";
  cout << "- Output ntuple: " << ntpOut << "\n";
  cout << "- Auxiliary ntuple: " << ntpAux << "\n";
  cout << "- YAML file: " << ymlFile << "\n";
  cout << "- Particle: " << particle << "\n";
  cout << "- applyAlias: " << applyAlias << "\n";
  cout << "- K smear branch name: " << kSmrBrName << "\n";
  cout << "- Pi smear branch name: " << piSmrBrName << "\n";
  cout << "- K vmu smear branch name: " << kSmrBrNameVmu << "\n";
  cout << "- Pi vmu smear branch name: " << piSmrBrNameVmu << "\n";
  cout << "- Smearing strategy: " << smr_mode << endl;

  // parse YAML config
  const auto ymlConfig       = YAML::LoadFile(ymlFile);
  const auto year            = parsedArgs["year"].as<string>();
  const auto ymlDirPath      = absDirPath(parsedArgs["config"].as<string>());
  const auto outputDirective = ymlConfig["weight_brs"];
  const auto filePrefix      = absDirPath(ymlFile);

  // snapshot option
  auto writeOpts  = ROOT::RDF::RSnapshotOptions{};
  writeOpts.fMode = "UPDATE";

  // Produce list of skims. Currently, only tagged yields are different.
  // For vmu, no skim cuts are applied.
  // TODO recover this from ymlConfig["skims"]
  const vector<TString> skims = {"iso", "1os", "2os", "dd", "vmu", "prot"};

  // Generate names of histograms to be imported from unfolded.root
  const vector<TString> histoWtNames =
      buildHistoWtNames(particle, skims, year, ymlConfig["tags"]);
  if (debug) {
    cout << "\nEfficiency histograms:\n";
    for (const auto& h : histoWtNames) {
      cout << "\t" << h << "\n";
    }
    cout << endl;
  }

  for (auto it = outputDirective.begin(); it != outputDirective.end(); it++) {
    if (debug) cout << "--------" << endl;
    const auto treeName    = it->first.as<string>();
    const auto ntpInTest   = new TFile(ntpIn.data());
    auto       ntpsToClean = vector<TFile*>{ntpInTest};

    const auto treeTest = dynamic_cast<TTree*>(ntpInTest->Get(treeName.data()));
    if (treeTest == nullptr) {
      cout << treeName << " doesn't exist in " << ntpIn << ". skipping..."
           << endl;
      continue;
    }

    // build a dataframe from input ntuple
    auto           df = static_cast<RNode>(RDataFrame(treeName, ntpIn));
    vector<string> outputBrNames{"runNumber", "eventNumber"};
    if (applyAlias) {
      df = defineBranch(df, particle, MU_BRANCH_DEFS, debug);
      // compute ETA
      df = df.Define("ETA",
                     [](double& p, double& pz) {
                       return 0.5 * TMath::Log((p + pz) / (p - pz));
                     },
                     {"P", "PZ"});
    }

    // read smearing factors from aux ntuple
    auto smrFacK     = vector<vector<double>>{};
    auto smrFacPi    = vector<vector<double>>{};
    auto smrFacKVmu  = vector<vector<double>>{};
    auto smrFacPiVmu = vector<vector<double>>{};
    getSmrFac(smrFacK, ntpAux, kSmrBrName);
    getSmrFac(smrFacPi, ntpAux, piSmrBrName);
    getSmrFac(smrFacKVmu, ntpAux, kSmrBrNameVmu);
    getSmrFac(smrFacPiVmu, ntpAux, piSmrBrNameVmu);

    // we can call these functions directly to get random smearing factors
    auto randSmrFacK     = getRandSmrHelper(smrFacK);
    auto randSmrFacPi    = getRandSmrHelper(smrFacPi);
    auto randSmrFacKVmu  = getRandSmrHelper(smrFacKVmu);
    auto randSmrFacPiVmu = getRandSmrHelper(smrFacPiVmu);

    // Recompute fit vars
    auto [dfOut, outputBrsFitVars] =
        defRestFrameVars(df, treeTest, randSmrFacK, randSmrFacPi, smr_mode);
    for (auto br : outputBrsFitVars) outputBrNames.emplace_back(br);
    df = dfOut;

    // Recompute fit vars for vmu
    auto [dfOutVmu, outputBrsFitVarsVmu] = defRestFrameVars(
        df, treeTest, randSmrFacKVmu, randSmrFacPiVmu, smr_mode, "_vmu");
    for (auto br : outputBrsFitVarsVmu) outputBrNames.emplace_back(br);
    df = dfOutVmu;

    // add species tags
    auto [directivesTags, outputBrsTags] =
        genTaggedCutDirective(ymlConfig["tags"]);
    df = defineBranch(df, ""s, directivesTags, debug);
    for (const auto& br : outputBrsTags) outputBrNames.emplace_back(br);

    // add all kinds of weights
    for (const auto& entry : it->second) {
      const auto histoPrefix    = entry["prefix"].as<string>();
      const auto histoFile      = filePrefix + "/" + entry["file"].as<string>();
      const auto weightBrPrefix = entry["name"].as<string>();

      if (debug) {
        cout << "\n\nHandling tree " << treeName << " from histos of prefix "
             << histoPrefix << " from file " << histoFile << endl;
      }

      // add weights required by misID weights
      auto ntpHisto = new TFile(histoFile.data(), "READ");
      ntpsToClean.emplace_back(ntpHisto);

      if (debug) {
        cout << "Generate transfer factors/DiF smearing weights for all species"
             << endl;
      }
      auto [dfHistos, outputBrsHistos, histos] = applyWtFromHistos(
          df, ntpHisto, histoPrefix, weightBrPrefix, histoWtNames, debug);
      for (auto br : outputBrsHistos) outputBrNames.emplace_back(br);

      // add the actual misID weights
      auto [directives, outputBrsWts] =
          genWtDirective(ymlConfig["tags"], weightBrPrefix, skims, debug);
      df = defineBranch(dfHistos, ""s, directives, debug);
      for (auto br : outputBrsWts) outputBrNames.emplace_back(br);
    }
    if (debug) {
      cout << "\n\nWriting branches to " << ntpOut << ":" << endl;
      for (auto br : outputBrNames) cout << "\t-" << br << endl;
    }
    df.Snapshot(treeName, ntpOut, outputBrNames, writeOpts);

    for (auto& n : ntpsToClean) delete n;
  }
}
