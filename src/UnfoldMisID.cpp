// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Mon Mar 28, 2022 at 12:24 AM -0400
//
// Description: unfolding efficiency calculator (U)

#include <iostream>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>
#include <cxxopts.hpp>

using namespace std;

////////////////////////
// Histo name helpers //
////////////////////////

vector<string> getTagNames(YAML::Node cfgTagged) {
  vector<string> result{};

  for (const auto& p : cfgTagged) {
    auto node   = p.first;
    auto keyRaw = node.as<string>();
    result.push_back(keyRaw);
  }

  return result;
}

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
    // input/output
    ("e,effHisto", "specify input ntuple containing efficiency histos")
    ("y,yldHisto", "specify input ntuple containing raw yield histos")
    ("c,config", "specify input YAML config file")
    ("o,output", "specify output folder")
    // flags
    ("targetParticle", "specify target particle for unfolding",
     cxxopts::value<string>()->default_value("mu"))
    ("outputHisto", "specify output histo name",
     cxxopts::value<string>()->default_value("unfolded.root"))
    // misc
    ("d,debug", "enable debug mode",
     cxxopts::value<bool>()->default_value("false"))
    ;
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  auto ymlConfig  = YAML::LoadFile(parsedArgs["config"].as<string>());
  auto ptclTagged = getTagNames(ymlConfig["tags"]);

  // debug output
  if (parsedArgs["debug"].as<bool>()) {
    cout << "The tagged species are:" << endl;
    for (const auto p : ptclTagged) cout << "  " << p << endl;
  }

  return 0;
}
