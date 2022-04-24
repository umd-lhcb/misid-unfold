// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Sun Apr 24, 2022 at 05:30 PM -0400
//
// Description: unfolding weights applyer (A)

#include <filesystem>
#include <iostream>
#include <regex>
#include <string>
#include <tuple>
#include <vector>

#include <TFile.h>
#include <TH3D.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>

#include <yaml-cpp/yaml.h>
#include <boost/range/join.hpp>
#include <cxxopts.hpp>

#include "utils.h"

using namespace std;
using ROOT::RDataFrame;
using ROOT::RDF::RNode;

///////////////////
// Configuration //
///////////////////

typedef vector<pair<string, string>> vPStrStr;
typedef vector<pair<regex, string>>  vPRegStr;

static vPStrStr MU_BRANCH_DEFS{
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

static vPRegStr CUT_REPLACE_RULES{{regex("&"), "&&"}, {regex("\\|"), "||"}};

/////////////////////
// General helpers //
/////////////////////

string absDirPath(string pathRaw) {
  auto path    = filesystem::path(pathRaw);
  auto dirPath = path.parent_path();
  return filesystem::absolute(dirPath).string();
}

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

////////////////////////////
// Helpers for event loop //
////////////////////////////

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
  auto           wtTargetParticle = "MuTag";
  auto           wtSmrParticles   = {"k", "pi"};

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
    delete treeTest;
    delete ntpInTest;

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

    cout << "Writing to " << ntpOut << endl;
    if (firstTree) {
      df.Snapshot(treeName, ntpOut, outputBrNames);
      firstTree = false;
    } else
      df.Snapshot(treeName, ntpOut, outputBrNames, writeOpts);

    // cleanups
    for (auto& h : histos) delete h;
    delete ntpHisto;
  }
}
