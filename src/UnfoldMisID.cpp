// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Mon Mar 28, 2022 at 09:08 PM -0400
//
// Description: unfolding efficiency calculator (U)

#include <cctype>
#include <iostream>
#include <map>
#include <regex>
#include <string>
#include <tuple>
#include <vector>

#include <TFile.h>
#include <TH3D.h>

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

////////////////////
// Config helpers //
////////////////////

vector<string> getKeyNames(YAML::Node node, string prefix = "",
                           string suffix = "", bool repl = true,
                           string replRegex = "^misid_") {
  vector<string> result{};

  for (auto it = node.begin(); it != node.end(); it++) {
    auto   keyRaw = it->first.as<string>();
    string key;
    if (repl)
      key = regex_replace(keyRaw, regex(replRegex), "");
    else
      key = keyRaw;
    result.emplace_back(prefix + key + suffix);
  }

  return result;
}

vector<string> getYldHistoNames(vector<string> ptcl, string suffix = "Tag") {
  vector<string> result{};
  for (auto pt : ptcl) result.emplace_back(pt + suffix);
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

tuple<vector<vector<float>>, vector<int>> getBins(YAML::Node cfgBinning) {
  vector<vector<float>> binnings{};
  vector<int>           nbins{};

  for (auto it = cfgBinning.begin(); it != cfgBinning.end(); it++) {
    vector<float> binEdges{};
    int           counter = -1;
    for (auto elem : it->second) {
      binEdges.emplace_back(elem.as<float>());
      counter += 1;
    }
    binnings.emplace_back(binEdges);
    nbins.emplace_back(counter);
  }

  return tuple<vector<vector<float>>, vector<int>>{binnings, nbins};
}

///////////////////
// Histo helpers //
///////////////////

map<string, TH3D*> prepOutHisto(vector<string>&        names,
                                vector<vector<float>>& bins) {
  map<string, TH3D*> result{};

  for (size_t idx = 0; idx != names.size(); idx++) {
    auto histoName = names[idx].data();

    auto nbinsX = bins[0].size() - 1;
    auto nbinsY = bins[1].size() - 1;
    auto nbinsZ = bins[2].size() - 1;

    auto xBins = bins[0].data();
    auto yBins = bins[1].data();
    auto zBins = bins[2].data();

    auto histo = new TH3D(histoName, histoName, nbinsX, xBins, nbinsY, yBins,
                          nbinsZ, zBins);
    result[histoName] = histo;
  }

  return result;  // in principle these pointers need deleting
}

map<string, TH3D*> loadHisto(TFile* ntpYld, TFile* ntpEff,
                             vector<string>         nameMeaYld,
                             vector<vector<string>> nameEff,
                             vector<string>         prefix) {
  map<string, TH3D*> result{};

  for (const auto& n : nameMeaYld) {
    for (auto pref : prefix) {
      auto name = (pref + "__" + n);
      cout << "Loading " << name << endl;

      auto histo = static_cast<TH3D*>(ntpYld->Get(name.data()));
      if (!histo) {
        cout << "Histogram " << name << " doesn't exist! terminate now..."
             << endl;
        exit(1);
      }

      result[name] = histo;
    }
  }

  for (const auto& row : nameEff) {
    for (const auto& n : row) {
      cout << "Loading " << n << endl;

      auto histo = static_cast<TH3D*>(ntpEff->Get(n.data()));
      if (!histo) {
        cout << "Histogram " << n << " doesn't exist! terminate now..." << endl;
        exit(1);
      }

      result[n] = histo;
    }
  }

  cout << "All histograms loaded." << endl;
  return result;
}

////////////
// Unfold //
////////////

void unfold(map<string, TH3D*> histoIn, map<string, TH3D*> histoOut,
            vector<int> nbins, vector<string> nameMeaYld,
            vector<vector<string>> nameEff, vector<string> prefix,
            bool debug = false) {
  int totSize = nameMeaYld.size();

  for (int x = 1; x <= nbins[0]; x++) {
    for (int y = 1; y <= nbins[1]; y++) {
      for (int z = 1; z <= nbins[2]; z++) {
        for (auto pref : prefix) {
          cout << "Working on " << pref << endl;

          // build yield array
          double arrYld[totSize];
          for (int idx = 0; idx != totSize; idx++) {
            auto name  = pref + "__" + nameMeaYld[idx];
            auto histo = histoIn[name];
            if (histo)
              arrYld[idx] = histoIn[name]->GetBinContent(x, y, z);
            else
              cout << name << " is a nullptr!" << endl;
          }

          // build response array (2D matrix)
          double arrRes[totSize][totSize];
          for (int idxRow = 0; idxRow != totSize; idxRow++) {
            for (int idxCol = 0; idxCol != totSize; idxCol++)
              arrRes[idxRow][idxCol] =
                  histoIn[nameEff[idxRow][idxCol]]->GetBinContent(x, y, z);
          }

          if (debug) {
            cout << "Bin index: x=" << x << " y=" << y << " z=" << z << endl;

            cout << "The measured yields are:" << endl << "  ";
            for (const auto& n : arrYld) cout << n << "\t";
            cout << endl;

            cout << "The response matrix is:" << endl;
            for (const auto& row : arrRes) {
              cout << "  ";
              for (const auto& n : row) cout << n << "\t";
              cout << endl;
            }
          }
        }
      }
    }
  }
}

void unfoldDryRun(vector<string> ptcl, vector<string> nameMeaYld,
                  vector<string> nameUnfYld, vector<vector<string>> nameEff,
                  vector<vector<float>> binnings, vector<int> nbins) {
  cout << "The tagged species are:" << endl;
  for (const auto p : ptcl) cout << "  " << p << endl;

  cout << "The measured yields are stored in these histos:" << endl;
  for (const auto h : nameMeaYld) cout << "  " << h << endl;

  cout << "The unfolded yields will be stored in these histos:" << endl;
  for (const auto h : nameUnfYld) cout << "  " << h << endl;

  cout << "The response matrix will be built from these histos:" << endl;
  for (const auto row : nameEff) {
    cout << "  ";
    for (const auto elem : row) cout << elem << "\t";
    cout << endl;
  }

  cout << "The binning is defined as:" << endl;
  for (const auto row : binnings) {
    cout << "  ";
    for (const auto elem : row) cout << elem << "\t";
    cout << endl;
  }

  cout << "The bin sizes are:" << endl;
  cout << "  ";
  for (const auto n : nbins) cout << n << "\t";
  cout << endl;
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
    ("d,debug", "enable debug mode",
     cxxopts::value<bool>()->default_value("false"))
    ("D,dryRun", "parse config and load histos, w/o unfolding",
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
  auto ymlConfig       = YAML::LoadFile(parsedArgs["config"].as<string>());
  auto ptclTagged      = getKeyNames(ymlConfig["tags"]);
  auto prefix          = getKeyNames(ymlConfig["input_ntps"]);
  auto histoNameMeaYld = getYldHistoNames(ptclTagged);
  auto histoNameUnfYld = getYldHistoNames(ptclTagged, "True");
  auto histoNameEff    = getEffHistoNames(ptclTagged);
  auto [histoBinSpec, histoBinSize] = getBins(ymlConfig["binning"]);

  // dry run
  if (parsedArgs["dryRun"].as<bool>()) {
    unfoldDryRun(ptclTagged, histoNameMeaYld, histoNameUnfYld, histoNameEff,
                 histoBinSpec, histoBinSize);
    return 0;
  }

  // open ntuples
  auto ntpYld = new TFile(parsedArgs["yldHisto"].as<string>().data());
  auto ntpEff = new TFile(parsedArgs["effHisto"].as<string>().data());

  // prepare histograms
  auto histoOut = prepOutHisto(histoNameUnfYld, histoBinSpec);
  auto histoIn =
      loadHisto(ntpYld, ntpEff, histoNameMeaYld, histoNameEff, prefix);

  cout << histoIn["D0__piTag"]->GetBinContent(1, 1, 1);

  // unfold
  auto debug = parsedArgs["debug"].as<bool>();
  unfold(histoIn, histoOut, histoBinSize, histoNameMeaYld, histoNameEff, prefix,
         debug);

  // cleanup
  for (auto h : histoOut) delete h.second;
  for (auto h : histoIn) delete h.second;
  delete ntpYld;
  delete ntpEff;

  return 0;
}
