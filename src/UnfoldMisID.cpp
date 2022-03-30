// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Wed Mar 30, 2022 at 10:55 AM -0400
//
// Description: unfolding efficiency calculator (U)

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>

#include <RooUnfoldBayes.h>
#include <RooUnfoldResponse.h>

#include <yaml-cpp/yaml.h>
#include <cxxopts.hpp>

#define DEBUG_OUT cout << __LINE__ << endl;

using namespace std;

//////////////
// Typedefs //
//////////////

typedef vector<string>         vStr;
typedef vector<vector<string>> vStrStr;

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

vector<double> histoToProb(const TH1D* histo) {
  vector<double> result{};
  auto           normFac = histo->GetEntries();
  for (int idx = 1; idx <= histo->GetNbinsX(); idx++)
    result.emplace_back(histo->GetBinContent(idx) / normFac);
  return result;
}

////////////////////
// Config helpers //
////////////////////

vStr getKeyNames(YAML::Node node, string prefix = "", string suffix = "") {
  vStr result{};

  for (auto it = node.begin(); it != node.end(); it++) {
    auto key = it->first.as<string>();
    result.emplace_back(prefix + key + suffix);
  }

  return result;
}

vStrStr getYldHistoNames(const vStr& ptcl, const vStr& prefix,
                         string suffix = "Tag") {
  vStrStr result{};

  for (auto pref : prefix) {
    vStr row{};
    for (const auto& pt : ptcl) row.emplace_back(pref + "__" + pt + suffix);
    result.emplace_back(row);
  }

  return result;
}

// These are a mess
vStr getEffHistoNames(const string ptcl1, const vStr& ptcl2,
                      string spc1 = "True", string spc2 = "Tag",
                      string prefix = "") {
  vStr result{};
  // e.g. muTrueToPiTag, with mu fixed
  for (auto p2 : ptcl2)
    result.emplace_back(prefix + ptcl1 + spc1 + "To" + capitalize(p2) + spc2);
  return result;
}

vStr getEffHistoNames(const vStr& ptcl1, const string ptcl2,
                      string spc1 = "True", string spc2 = "Tag",
                      string prefix = "") {
  vStr result{};
  // e.g. muTrueToPiTag, with Pi fixed
  for (auto p1 : ptcl1)
    result.emplace_back(prefix + p1 + spc1 + "To" + capitalize(ptcl2) + spc2);
  return result;
}

vStrStr getEffHistoNames(const vStr& ptcl1, const vStr& ptcl2,
                         string spc1 = "True", string spc2 = "Tag",
                         string prefix = "") {
  vStrStr result{};

  // Here we (indirectly) define the rows and columns of the response matrix, be
  // careful!
  for (auto p2 : ptcl2) {  // w/ default values, we are looping over Tag
    // for each row, we hold tag unchanged
    auto row = getEffHistoNames(ptcl1, p2, spc1, spc2, prefix);
    result.emplace_back(row);
  }

  return result;
}
// end of mess

tuple<vector<vector<float>>, vector<int>> getBins(YAML::Node cfgBinning) {
  vector<vector<float>> binnings{};
  vector<int>           nbins{};

  for (auto it = cfgBinning.begin(); it != cfgBinning.end(); it++) {
    vector<float> binEdges{};
    int           counter = -1;
    for (const auto elem : it->second) {
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

map<string, TH3D*> prepOutHisto(const vStrStr&         names,
                                vector<vector<float>>& bins) {
  map<string, TH3D*> result{};

  for (const auto vec : names) {
    for (const auto n : vec) {
      cout << "Preparing output histo: " << n << endl;
      auto histoName = n.data();

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
  }

  cout << "All output histos prepared." << endl;
  return result;  // in principle these pointers need deleting
}

map<string, TH3D*> loadHisto(const map<TFile*, vStrStr>& directive) {
  map<string, TH3D*> result{};

  for (const auto& [ntp, val] : directive) {
    for (const auto& row : val) {
      for (const auto& n : row) {
        cout << "Loading " << n << endl;

        auto histo = static_cast<TH3D*>(ntp->Get(n.data()));
        if (!histo) {
          cout << "Histogram " << n << " doesn't exist! terminate now..."
               << endl;
          exit(1);
        }

        result[n] = histo;
      }
    }
  }

  cout << "All histograms loaded." << endl;
  return result;
}

TH3D* loadSingleHisto(map<string, TH3D*>& mapHisto, string nameHisto) {
  if (mapHisto.find(nameHisto) == mapHisto.end()) {
    cout << "Histo " << nameHisto << " is not defined! Terminate now..."
         << endl;
    terminate();  // this is serious!
  }
  return mapHisto[nameHisto];
}

////////////
// Unfold //
////////////

void unfold(map<string, TH3D*> histoIn, map<string, TH3D*> histoOut,
            vStrStr nameMeaYld, vStrStr nameEff, vStr nameMuEff,
            vStrStr nameUnfYld, vStrStr nameUnfEff, bool debug = false,
            int numOfIter = 4) {
  int totSize = nameMeaYld[0].size();

  // These are used to stored measured yields (a vector) and response matrix (a
  // 2D matrix)
  auto histMea = new TH1D("histMea", "histMea", totSize, 0, totSize);
  auto histRes =
      new TH2D("histRes", "histRes", totSize, 0, totSize, totSize, 0, totSize);
  auto histInv = new TH2D("histInv", "histInv", totSize, 0, totSize, totSize, 0,
                          totSize);  // conceptually inverted matrix of histRes
  auto histProb = new TH1D("histProb", "histProb", totSize, 0, totSize);

  // This is used to provide dimension info for response matrix
  auto histTrue = new TH1D("histTrue", "histTrue", totSize, 0, totSize);
  for (int i = 1; i <= totSize; i++) histTrue->SetBinContent(i, 1);

  // Figure out binning from one of the input histograms
  vector<int> nbins{};
  for (const auto pair : histoIn) {
    auto ntp = pair.second;
    nbins.emplace_back(ntp->GetNbinsX());
    nbins.emplace_back(ntp->GetNbinsY());
    nbins.emplace_back(ntp->GetNbinsZ());
    break;
  }

  // Main unfolding procedure
  for (int x = 1; x <= nbins[0]; x++) {
    for (int y = 1; y <= nbins[1]; y++) {
      for (int z = 1; z <= nbins[2]; z++) {
        for (int idxPref = 0; idxPref != nameMeaYld.size(); idxPref++) {
          // build yield vector
          for (int idx = 0; idx != totSize; idx++) {
            auto histo = loadSingleHisto(histoIn, nameMeaYld[idxPref][idx]);
            histMea->SetBinContent(idx + 1, histo->GetBinContent(x, y, z));
          }

          // build response matrix (2D matrix)
          for (int idxTag = 0; idxTag != totSize; idxTag++) {
            for (int idxTrue = 0; idxTrue != totSize; idxTrue++) {
              auto histo = loadSingleHisto(histoIn, nameEff[idxTag][idxTrue]);
              auto eff   = histo->GetBinContent(x, y, z);
              if (isnan(eff) || isinf(eff)) eff = 0.0;
              histRes->SetBinContent(idxTag + 1, idxTrue + 1, eff);
            }
          }

          // perform unfolding to get unfolded ("true") yield
          RooUnfoldResponse resp(nullptr, histTrue, histRes);
          RooUnfoldBayes    unfoldWorker(&resp, histMea, numOfIter);
          auto              histUnf = static_cast<TH1D*>(unfoldWorker.Hreco());

          // Save unfolded yields
          for (int idx = 0; idx != totSize; idx++) {
            auto name = nameUnfYld[idxPref][idx];
            auto yld  = histUnf->GetBinContent(idx + 1);
            if (isnan(yld) || isinf(yld)) {
              cout << "Warning: naN or inf detected for " << name << endl;
              yld = 0;
            }
            auto histo = loadSingleHisto(histoOut, name);
            cout << "Writing unfolded yield to " << name << endl;
            ;
            histo->SetBinContent(x, y, z, yld);
          }

          // Compute unfolded ("true") probability
          auto probTag  = histoToProb(histMea);
          auto probTrue = histoToProb(histUnf);

          for (int idxTag = 0; idxTag != totSize; idxTag++) {
            auto wtTagToMuTag = 1.0;

            for (int idxTrue = 0; idxTrue != totSize; idxTrue++) {
              // from pidcalib we have true -> tag
              auto probTrueToTag =
                  histRes->GetBinContent(idxTag + 1, idxTrue + 1);
              if (isnan(probTrueToTag)) probTrueToTag = 0.0;

              // probability: tag -> true
              auto probTagToTrue =
                  probTrueToTag * probTrue[idxTrue] / probTag[idxTag];
              if (isnan(probTagToTrue)) probTagToTrue = 0.0;
              histInv->SetBinContent(idxTrue + 1, idxTag + 1, probTrueToTag);

              // now contract with the mu misID eff (true -> mu tag)
              auto histo = loadSingleHisto(histoIn, nameMuEff[idxTrue]);
              auto probTrueToMuTag = histo->GetBinContent(x, y, z);
              if (isnan(probTrueToMuTag)) probTrueToMuTag = 0.0;
              wtTagToMuTag = probTagToTrue * probTrueToMuTag;
            }

            // we use idxTag as the second index, which checks out
            auto name  = nameUnfEff[idxPref][idxTag];
            auto histo = loadSingleHisto(histoOut, name);
            cout << "Writing final contracted efficiency to " << name << endl;
            histo->SetBinContent(x, y, z, wtTagToMuTag);
          }

          if (debug) {
            cout << "Bin index: x=" << x << " y=" << y << " z=" << z << endl;

            cout.precision(4);
            cout << fixed;

            cout << "The yields are (top: measured; bot: unfolded):" << endl;
            for (int idx = 1; idx <= totSize; idx++)
              cout << setw(12) << histMea->GetBinContent(idx);
            cout << endl;
            for (int idx = 1; idx <= totSize; idx++)
              cout << setw(12) << histUnf->GetBinContent(idx);
            cout << endl;

            cout << "The response matrix is:" << endl;
            for (int idxRow = 1; idxRow <= totSize; idxRow++) {
              for (int idxCol = 1; idxCol <= totSize; idxCol++)
                cout << setw(8) << histRes->GetBinContent(idxRow, idxCol);
              cout << endl;
            }

            cout << "The tag -> true efficiencies are (row: fixed tag; col: "
                    "fixed true):"
                 << endl;
            for (int idxRow = 1; idxRow <= totSize; idxRow++) {
              for (int idxCol = 1; idxCol <= totSize; idxCol++)
                cout << setw(8) << histInv->GetBinContent(idxRow, idxCol);
              cout << endl;
            }
          }
        }
      }
    }
  }

  // cleanups
  delete histMea;
  delete histRes;
  delete histInv;
  delete histTrue;
}

void unfoldDryRun(vStr ptcl, vStrStr nameMeaYld, vStrStr nameUnfYld,
                  vStrStr nameEff, vStrStr nameUnfEff, vStr nameMuEff,
                  vector<vector<float>> binnings, vector<int> nbins) {
  cout << "The tagged species are:" << endl;
  for (const auto& p : ptcl) cout << "  " << p << endl;

  cout << "The measured yields are stored in these histos:" << endl;
  for (const auto& row : nameMeaYld) {
    cout << "  ";
    for (const auto& elem : row) cout << setw(12) << elem;
    cout << endl;
  }

  cout << "The unfolded yields will be stored in these histos:" << endl;
  for (const auto& row : nameUnfYld) {
    cout << "  ";
    for (const auto& elem : row) cout << setw(12) << elem;
    cout << endl;
  }

  cout << "The unfolded efficiencies will be stored in these histos:" << endl;
  for (const auto& row : nameUnfEff) {
    cout << "  ";
    for (const auto& elem : row) cout << setw(19) << elem;
    cout << endl;
  }

  cout << "The response matrix will be built from these histos:" << endl;
  for (const auto& row : nameEff) {
    cout << "  ";
    for (const auto& elem : row) cout << setw(16) << elem;
    cout << endl;
  }

  cout << "The Mu efficiencies will be loaded from these histos:" << endl;
  for (const auto& h : nameMuEff) cout << "  " << h << endl;

  cout << "The binning is defined as:" << endl;
  for (const auto& row : binnings) {
    cout << "  ";
    for (const auto& elem : row) cout << elem << "\t";
    cout << endl;
  }

  cout << "The bin sizes are:" << endl;
  cout << "  ";
  for (const auto& n : nbins) cout << n << "\t";
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
    ("o,output", "specify output folder", cxxopts::value<string>())
    // flags (typically don't configure these)
    ("i,iteration", "specify number of unfolding iterations",
     cxxopts::value<int>()->default_value("4"))
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
  auto ymlConfig  = YAML::LoadFile(parsedArgs["config"].as<string>());
  auto ptclTarget = parsedArgs["targetParticle"].as<string>();
  auto ptclList   = getKeyNames(ymlConfig["tags"]);
  auto prefix     = getKeyNames(ymlConfig["input_ntps"]);

  auto histoNameMeaYld = getYldHistoNames(ptclList, prefix);
  auto histoNameUnfYld = getYldHistoNames(ptclList, prefix, "True");
  auto histoNameEff    = getEffHistoNames(ptclList, ptclList);
  auto histoNameMuEff  = getEffHistoNames(
      ptclList, ptclTarget);  // e.g. kTrueToMuTag, where Mu is fixed
  auto [histoBinSpec, histoBinSize] = getBins(ymlConfig["binning"]);

  vStrStr histoNameUnfEff{};
  for (auto pref : prefix)
    histoNameUnfEff.emplace_back(
        getEffHistoNames(ptclList, ptclTarget, "Tag", "Tag", pref + "__"));

  // dry run
  if (parsedArgs["dryRun"].as<bool>()) {
    unfoldDryRun(ptclList, histoNameMeaYld, histoNameUnfYld, histoNameEff,
                 histoNameUnfEff, histoNameMuEff, histoBinSpec, histoBinSize);
    return 0;
  }

  // open ntuples
  auto ntpYld = new TFile(parsedArgs["yldHisto"].as<string>().data());
  auto ntpEff = new TFile(parsedArgs["effHisto"].as<string>().data());

  auto outputFilename = parsedArgs["output"].as<string>() + "/" +
                        parsedArgs["outputHisto"].as<string>();
  auto ntpOut = new TFile(outputFilename.data(), "RECREATE");

  // prepare histograms
  vStrStr histoNameOut{};
  set_union(histoNameUnfYld.begin(), histoNameUnfYld.end(),
            histoNameUnfEff.begin(), histoNameUnfEff.end(),
            back_inserter(histoNameOut));
  auto histoOut = prepOutHisto(histoNameOut, histoBinSpec);

  vStrStr histoNameEffComb(histoNameEff);
  histoNameEffComb.emplace_back(histoNameMuEff);
  auto histoIn =
      loadHisto({{ntpYld, histoNameMeaYld}, {ntpEff, histoNameEffComb}});
  // ^^the argumnets are maps of TFile* -> vStrStr, so we have to combine names
  // from the same file!

  // unfold
  auto debug     = parsedArgs["debug"].as<bool>();
  auto numOfIter = parsedArgs["iteration"].as<int>();
  unfold(histoIn, histoOut, histoNameMeaYld, histoNameEff, histoNameMuEff,
         histoNameUnfYld, histoNameUnfEff, debug, numOfIter);

  // save output
  for (const auto& pair : histoOut)
    ntpOut->WriteObject(pair.second, pair.first.data());

  // cleanup
  for (auto& h : histoOut) delete h.second;
  for (auto& h : histoIn) delete h.second;
  delete ntpOut;
  delete ntpYld;
  delete ntpEff;

  return 0;
}
