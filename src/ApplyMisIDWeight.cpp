// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Mon Apr 18, 2022 at 07:47 PM -0400
//
// Description: unfolding weights applyer (A)

#include <filesystem>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include <TMath.h>
#include <ROOT/RDataFrame.hxx>

#include <yaml-cpp/yaml.h>
#include <cxxopts.hpp>

using namespace std;
using ROOT::RDataFrame;
using ROOT::RDF::RNode;

///////////////////
// Configuration //
///////////////////

typedef vector<pair<string, string>> vPStrStr;

static vPStrStr BRANCH_ALIASES{
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

/////////////////////
// General helpers //
/////////////////////

string absDirPath(string pathRaw) {
  auto path    = filesystem::path(pathRaw);
  auto dirPath = path.parent_path();
  return filesystem::absolute(dirPath).string();
}

////////////////////////////
// Helpers for event loop //
////////////////////////////

// Idea stolen from:
// https://root-forum.cern.ch/t/running-rdataframes-define-in-for-loop/32484/2
auto setBranchAlias(RNode df, string particle = "mu",
                    const vPStrStr& rules = BRANCH_ALIASES, int idx = 0) {
  // auto df = init_df.Alias(alias, particle+"_"+raw);
  if (rules.size() == idx) return df;
  return setBranchAlias(
      df.Alias(rules[idx].first, particle + "_" + rules[idx].second), particle,
      rules, idx + 1);
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

  for (auto it = weightBrs.begin(); it != weightBrs.end(); it++) {
    auto histoPrefix  = it->first.as<string>();
    auto histoFile    = it->second["file"].as<string>();
    auto treeName     = it->second["tree"].as<string>();
    auto weightBrName = it->second["name"].as<string>();
    histoFile         = filePrefix + "/" + histoFile;

    cout << "Handling tree " << treeName << " with output branch name "
         << weightBrName << " from histos of prefix " << histoPrefix
         << " from file " << histoFile << endl;

    auto dfInit = RDataFrame(treeName, ntpIn);

    RNode df = static_cast<RNode>(dfInit);
    if (applyAlias) df = setBranchAlias(dfInit, particle);

    // compute ETA
    if (applyAlias)
      df = df.Define("ETA",
                     [](double& p, double& pz) {
                       return 0.5 * TMath::Log((p + pz) / (p - pz));
                     },
                     {"P", "PZ"});
  }
}
