// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Sun Mar 27, 2022 at 11:29 PM -0400
//
// Description: unfolding efficiency calculator (U)

#include <cxxopts.hpp>
#include <string>

using namespace std;

////////////////////////
// Histo name helpers //
////////////////////////

//////////
// Main //
//////////

int main(int argc, char **argv) {
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
    ("targetParticle", "specify target partcle for unfolding",
     cxxopts::value<string>()->default_value("mu"))
    ("outputHisto", "specify output histo name",
     cxxopts::value<string>()->default_value("unfolded.root"))
    // misc
    ("d,debug", "enable debug mode",
     cxxopts::value<bool>()->default_value("false"))
    ;
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);

  return 0;
}
