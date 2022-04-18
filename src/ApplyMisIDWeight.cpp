// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Sun Apr 17, 2022 at 08:13 PM -0400
//
// Description: unfolding weights applyer (A)

#include <iostream>
#include <string>
#include <tuple>
#include <vector>

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
    {"InMuonAcc", "InMuonAcc"}};

static vector<string> DEFAULT_TREES{"tree"};

/////////////////////
// Aliases helpers //
/////////////////////

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

/////////////////////////////////
// Helpers for event selection //
/////////////////////////////////

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
    // I/O
    ("h,histo", "specify ntuple containing transfer factor histos.",
     cxxopts::value<string>())
    ("i,input", "specify input ntuple", cxxopts::value<string>())
    ("o,input", "specify output ntuple", cxxopts::value<string>())
    // flags
    ("t,trees", "specify trees to process.",
     cxxopts::value<vector<string>>())
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

  auto ntpIn  = parsedArgs["input"].as<string>();
  auto ntpOut = parsedArgs["output"].as<string>();
  auto ntpFac = parsedArgs["histo"].as<string>();

  auto trees = DEFAULT_TREES;
  if (parsedArgs.count("trees"))
    trees = parsedArgs["trees"].as<vector<string>>();

  for (const auto t : trees) {
    cout << "Handling tree " << t << endl;
    auto dfInit = RDataFrame(t, ntpIn);
  }
}
