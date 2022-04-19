// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Mon Apr 18, 2022 at 10:55 PM -0400
//
// Description: unfolding weights applyer (A)

#include <filesystem>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include <TFile.h>
#include <TH3D.h>
#include <TMath.h>
#include <TString.h>
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

/////////////////////
// General helpers //
/////////////////////

string absDirPath(string pathRaw) {
  auto path    = filesystem::path(pathRaw);
  auto dirPath = path.parent_path();
  return filesystem::absolute(dirPath).string();
}

string capitalize(string str) {
  for (auto& s : str) {
    s = toupper(s);
    break;
  }
  return str;
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

////////////////////////////
// Helpers for event loop //
////////////////////////////

// Idea stolen from:
// https://root-forum.cern.ch/t/running-rdataframes-define-in-for-loop/32484/2
auto defineBranch(RNode df, string particle = "mu",
                  const vPStrStr& rules = MU_BRANCH_DEFS, int idx = 0) {
  // auto df = init_df.Alias(alias, particle+"_"+raw);
  if (rules.size() == idx) return df;

  auto inputBrName = rules[idx].second;
  if (particle != ""s) inputBrName = particle + "_" + inputBrName;

  return defineBranch(df.Define(rules[idx].first, inputBrName), particle, rules,
                      idx + 1);
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
  for (auto it = weightBrs.begin(); it != weightBrs.end(); it++) {
    auto histoPrefix  = TString(it->first.as<string>());
    auto histoFile    = it->second["file"].as<string>();
    auto treeName     = it->second["tree"].as<string>();
    auto weightBrName = it->second["name"].as<string>();
    histoFile         = filePrefix + "/" + histoFile;

    cout << "Handling tree " << treeName << " from histos of prefix "
         << histoPrefix << " from file " << histoFile << endl;

    ntpHisto           = new TFile(TString(histoFile), "READ");
    auto  histoWtNames = buildHistoWtNames(particle, ymlConfig["tags"]);
    auto  dfInit       = RDataFrame(treeName, ntpIn);
    RNode df           = static_cast<RNode>(dfInit);
    vector<string> outputBrNames{"runNumber", "eventNumber"};

    if (applyAlias) {
      df = defineBranch(dfInit, particle);
      // compute ETA
      df = df.Define("ETA",
                     [](double& p, double& pz) {
                       return 0.5 * TMath::Log((p + pz) / (p - pz));
                     },
                     {"P", "PZ"});
    }

    for (const auto h : histoWtNames) {
      auto histoName = histoPrefix + "__" + h;
      auto histoWt   = static_cast<TH3D*>(ntpHisto->Get(histoName));
      histos.emplace_back(histoWt);
      cout << "  Loading histo " << histoName << endl;

      auto brName = weightBrName + "_" + string(h);
      outputBrNames.emplace_back(brName);
      cout << "  Generating " << brName << "..." << endl;
      df = df.Define(brName,
                     [histoWt](double x, double y, double z) {
                       auto binIdx = histoWt->FindFixBin(x, y, z);
                       return histoWt->GetBinContent(binIdx);
                     },
                     {"P", "ETA", "nTracks"});
    }

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
