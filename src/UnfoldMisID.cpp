// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Mon Mar 28, 2022 at 12:38 PM -0400
//
// Description: unfolding efficiency calculator (U)

#include <cctype>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <cxxopts.hpp>

using namespace std;

/////////////////////
// General helpers //
/////////////////////

string capitalize(string str) {
  for (auto& s : str) {
    s = toupper(s);
    break;
  }
  return str;
}

////////////////////////
// Histo name helpers //
////////////////////////

vector<string> getTagNames(YAML::Node cfgTagged) {
  vector<string> result{};

  for (auto it = cfgTagged.begin(); it != cfgTagged.end(); it++) {
    auto keyRaw = it->first.as<string>();
    result.emplace_back(regex_replace(keyRaw, regex("^misid_"),
                                      ""));  // remove the heading 'misid_'
  }

  return result;
}

vector<string> getMeaYldHistoNames(vector<string> ptcl) {
  vector<string> result{};

  for (auto pt : ptcl) {
    result.emplace_back(pt + "Tag");
  }

  return result;
}

vector<vector<string>> getEffHistoNames(vector<string> ptcl) {
  vector<vector<string>> result{};

  // Here we (indirectly) define the rows and columns of the response matrix, be
  // careful!
  for (auto ptTag : ptcl) {
    vector<string> row{};
    for (auto ptTrue : ptcl) {
      row.emplace_back(ptTrue + "TrueTo" + capitalize(ptTag) + "Tag");
    }
    result.emplace_back(row);
  }

  return result;
}

///////////////////
// Histo loaders //
///////////////////

//////////
// Main //
//////////

int main(int argc, char** argv) {
  cxxopts::Options argOpts("UnfoldMisID",
                           "unfolding efficiency calculator (U).");

  // clang-format off
  argOpts.add_options()
    // general
    ("h,help", "print help")
    ("d,dryRun", "parse config and load histos, w/o unfolding",
     cxxopts::value<bool>()->default_value("false"))
    // input/output
    ("e,effHisto", "specify input ntuple containing efficiency histos",
     cxxopts::value<string>())
    ("y,yldHisto", "specify input ntuple containing measured yield histos",
     cxxopts::value<string>())
    ("c,config", "specify input YAML config file",
     cxxopts::value<string>())
    ("o,output", "specify output folder")
    // flags (typically don't configure these)
    ("targetParticle", "specify target particle for unfolding",
     cxxopts::value<string>()->default_value("mu"))
    ("outputHisto", "specify output histo name",
     cxxopts::value<string>()->default_value("unfolded.root"))
    ;
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  // parse YAML config
  auto ymlConfig    = YAML::LoadFile(parsedArgs["config"].as<string>());
  auto ptclTagged   = getTagNames(ymlConfig["tags"]);
  auto histoNameYld = getMeaYldHistoNames(ptclTagged);
  auto histoNameEff = getEffHistoNames(ptclTagged);

  // dry run
  if (parsedArgs["dryRun"].as<bool>()) {
    cout << "The tagged species are:" << endl;
    for (const auto p : ptclTagged) cout << "  " << p << endl;

    cout << "The measured yields are stored in these histos:" << endl;
    for (const auto h : histoNameYld) cout << "  " << h << endl;

    cout << "The response matrix will be built from these histos:" << endl;
    for (const auto row : histoNameEff) {
      cout << "  ";
      for (const auto elem : row) cout << elem << "\t";
      cout << endl;
    }

    return 0;
  }

  return 0;
}
