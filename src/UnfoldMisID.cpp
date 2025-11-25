// Author: Yipeng Sun
// License: BSD 2-clause
// Last Change: Sun Apr 24, 2022 at 05:30 PM -0400
//
// Description: unfolding efficiency calculator (U)

#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
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

#include "utils.h"

#define DEBUG_OUT cout << __LINE__ << endl;

using namespace std;

//////////////
// Typedefs //
//////////////

typedef vector<string>        vStr;
typedef vector<vector<float>> vFltFlt;

/////////////////////
// General helpers //
/////////////////////

vector<double> histoToProb(const TH1D* histo) {
  double normFac = 0;
  for (int idx = 1; idx <= histo->GetNbinsX(); idx++) {
    normFac += histo->GetBinContent(idx);
  }
  // we loop twice. it's stupid but it works
  vector<double> result{};
  result.reserve(histo->GetNbinsX());
  for (int idx = 1; idx <= histo->GetNbinsX(); idx++) {
    result.emplace_back(histo->GetBinContent(idx) / normFac);
  }
  return result;
}

////////////////////
// Config helpers //
////////////////////
// ptcl: particle

vStr getKeyNames(YAML::Node node, string prefix = "", string suffix = "") {
  vStr result{};

  for (auto it = node.begin(); it != node.end(); it++) {
    auto key = it->first.as<string>();
    result.emplace_back(prefix + key + suffix);
  }

  return result;
}

// Input yield histograms read from tagged.root
string getHistoName(const string& prefix, const string& ptcl,
                    const string& skim, const string& year,
                    const string descr = "Tag") {
  return ""s + prefix + "__" + ptcl + descr + "__" + skim + "__" + year;
}

// Input efficiecy histograms read from merged.root
string getHistoEffName(const string& ptcl1, const string& ptcl2,
                       const string& year, const string& descr1 = "True",
                       const string& descr2 = "Tag") {
  return ""s + ptcl1 + descr1 + "To" + capitalize(ptcl2) + descr2 + "_" + year;
}

// Output histograms to be saved in unfolded.root
string getHistoEffName(const string& prefix, const string& ptcl1,
                       const string& ptcl2, const string& descr1,
                       const string& descr2, const string& skim,
                       const string& year) {
  const string skim_suffix = skim.size() > 0 ? "_" + skim : "";
  return ""s + prefix + "__" +
         getHistoEffName(ptcl1, ptcl2, year, descr1, descr2) + skim_suffix;
}

tuple<vFltFlt, vector<int>> getBins(YAML::Node cfgBinning) {
  vFltFlt     binnings{};
  vector<int> nbins{};

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

  return tuple<vFltFlt, vector<int>>{binnings, nbins};
}

///////////////////
// Histo helpers //
///////////////////

void printResGeneric(const TH2D*                             res,
                     function<double(const TH2D*, int, int)> getter) {
  auto           nbinsX = res->GetNbinsX();
  auto           nbinsY = res->GetNbinsY();
  vector<string> badElem{};

  cout.precision(4);
  cout << fixed;

  for (int idxCol = 1; idxCol <= nbinsX; idxCol++) {
    for (int idxRow = 1; idxRow <= nbinsY; idxRow++) {
      auto elem = getter(res, idxCol, idxRow);
      cout << setw(8) << elem;
      if (elem < 0)
        badElem.emplace_back("(" + to_string(idxRow) + ", " +
                             to_string(idxCol) + ") = " + to_string(elem));
    }
    cout << endl;
  }

  if (badElem.size() > 0) {
    cout << "  Matrix contains negative element(s)" << endl;
    for (const auto& err : badElem) cout << "    " << err << endl;
  }
}

void printResVal(const TH2D* res) {
  auto getter = [](const TH2D* res, int x, int y) {
    return res->GetBinContent(x, y);
  };
  printResGeneric(res, getter);
}

void printResErr(const TH2D* res) {
  auto getter = [](const TH2D* res, int x, int y) {
    return res->GetBinError(x, y);
  };
  printResGeneric(res, getter);
}

auto getHistoInHelper(TFile* ntpYld, TFile* ntpEff) {
  auto mapHisto = make_shared<map<string, shared_ptr<TH3D>>>();

  // NOTE: We have to capture by copy, s.t. the memory is not deallocated at
  //       after the wrapper function finishes executing
  return [=](const string& key, const string& skim, bool debug = false) {
    const bool use_vmu_eff =
        (skim == "vmu") && (key.find("MuTag") != string::npos);
    const string histo_name = use_vmu_eff ? key + "_vmu" : key;
    if (!mapHisto->count(histo_name)) {
      TH3D* histo;
      if (key.find("__") != string::npos) {
        // this is a yield histo
        if (debug) cout << "DEBUG Getting yield histo " << histo_name << endl;
        histo = static_cast<TH3D*>(ntpYld->Get(histo_name.data()));
      } else {
        if (debug) {
          cout << "DEBUG Getting " << skim << " efficiency histo " << histo_name
               << endl;
        }
        histo = static_cast<TH3D*>(ntpEff->Get(histo_name.data()));
      }

      if (histo == nullptr) {
        cout << "ERROR Histogram " << histo_name
             << " doesn't exist! terminate now..." << endl;
        terminate();
      }

      auto histoPtr = shared_ptr<TH3D>(histo);
      mapHisto->emplace(histo_name, histoPtr);
      return histoPtr;
    }

    // the key already exists
    if (debug)
      cout << "DEBUG Getting cached eff histo " << histo_name
           << " (skim = " << skim << ")" << endl;
    return mapHisto->at(histo_name);
  };
}

auto getHistoOutHelper(map<string, shared_ptr<TH3D>>* mapHisto,
                       const vFltFlt&                 binnings) {
  return [=](const string& key) {
    if (!mapHisto->count(key)) {
      auto nbinsX = binnings[0].size() - 1;
      auto nbinsY = binnings[1].size() - 1;
      auto nbinsZ = binnings[2].size() - 1;

      auto xBins = binnings[0].data();
      auto yBins = binnings[1].data();
      auto zBins = binnings[2].data();

      auto histoPtr = make_shared<TH3D>(TH3D(
          key.data(), key.data(), nbinsX, xBins, nbinsY, yBins, nbinsZ, zBins));
      mapHisto->emplace(key, histoPtr);
      return histoPtr;
    }

    // the key already exists
    return mapHisto->at(key);
  };
}

////////////
// Unfold //
////////////

void ensureUnitarity(TH2D* res, vStr ptcls, bool debug = false) {
  auto nbinsX = res->GetNbinsX();
  auto nbinsY = res->GetNbinsY();

  cout.precision(4);
  cout << fixed;

  if (debug) {
    cout << "The raw true -> tag matrix is (row: fixed tag; col: fixed true):"
         << endl;
    printResVal(res);
    cout << "The raw true -> tag matrix error is (row: fixed tag; col: fixed "
            "true):"
         << endl;
    printResErr(res);
  }

  for (int y = 1; y <= nbinsY; y++) {
    double prob = 0.0;
    for (int x = 1; x < nbinsX; x++) {
      prob += res->GetBinContent(x, y);
    }
    if (abs((1 - prob) - res->GetBinContent(nbinsX, y)) > 1e-6) {
      cout << "INFO Changing " << ptcls[y - 1] << "TrueTo" << ptcls[nbinsX - 1]
           << "Tag eff from " << res->GetBinContent(nbinsX, y) << " to "
           << 1 - prob << endl;
    }

    res->SetBinContent(nbinsX, y, 1 - prob);
  }

  if (debug) {
    cout << "The fixed true -> tag matrix is (row: fixed tag; col: fixed true):"
         << endl;
    printResVal(res);
  }
}

template <typename F1, typename F2>
void unfold(const vStr& prefix, const vStr& ptcls, const vStr& skims,
            const vector<int>& nbins, F1& histoInGetter, F2& histoOutGetter,
            const string& year, map<string, int>& count_empty_tagged,
            bool debug = false, int numOfIter = 4) {
  const int totSize = ptcls.size();

  // These are used to stored measured yields (a vector) and response matrix (a
  // 2D matrix)
  TH1D histMea(("histMea_" + year).c_str(), ("histMea_" + year).c_str(),
               totSize, 0, totSize);
  TH2D histRes(("histRes_" + year).c_str(), ("histRes_" + year).c_str(),
               totSize, 0, totSize, totSize, 0, totSize);
  TH2D histInv(("histInv_" + year).c_str(), ("histInv_" + year).c_str(),
               totSize, 0, totSize, totSize, 0,
               totSize);  // conceptually inverted matrix of histRes
  TH1D histProb(("histProb_" + year).c_str(), ("histProb_" + year).c_str(),
                totSize, 0, totSize);

  // This is used to provide dimension info for response matrix only
  TH1D histDim(("histDim_" + year).c_str(), ("histDim_" + year).c_str(),
               totSize, 0, totSize);
  for (int i = 1; i <= totSize; i++) histDim.SetBinContent(i, 1);

  // Main unfolding procedure
  for (int x = 1; x <= nbins[0]; x++) {
    for (int y = 1; y <= nbins[1]; y++) {
      for (int z = 1; z <= nbins[2]; z++) {
        for (auto& pref : prefix) {
          for (const auto& skim : skims) {
            cout << "\nINFO Unfolding bin " << x << " " << y << " " << z
                 << " of " << pref << " " << skim << " " << year << endl;

            // build yield vector
            for (int idx = 0; idx != totSize; idx++) {
              auto name  = getHistoName(pref, ptcls[idx], skim, year);
              auto histo = histoInGetter(name, skim, debug);
              histMea.SetBinContent(idx + 1, histo->GetBinContent(x, y, z));
            }

            // build response matrix (2D matrix)
            // NOTE: for a 2D array, the indexing is this:
            //         array[x][y]
            //       In our case, true -> tag translates to:
            //         true -> y index
            //         tag  -> x index
            for (int idxTag = 0; idxTag != totSize; idxTag++) {
              for (int idxTrue = 0; idxTrue != totSize; idxTrue++) {
                auto name =
                    getHistoEffName(ptcls[idxTrue], ptcls[idxTag], year);
                auto histo = histoInGetter(name, skim, debug);
                auto eff   = histo->GetBinContent(x, y, z);
                auto err   = histo->GetBinError(x, y, z);

                if (debug)
                  cout << "  Loading efficiency from " << name << ", got "
                       << eff << endl;

                if (isnan(eff) || isinf(eff)) {
                  cout << "WARNING Invalid efficiency (" << eff
                       << "). Setting it to 0." << endl;
                  eff = 0.;
                }

                if (eff < 0.) {
                  cout << "WARNING Negative efficiency (" << eff
                       << "). Setting it to 0." << endl;
                  eff = 0.;
                }

                // Avoid signed zero in log file, which may be confused with a
                // small negative efficiency. Note: floating point +0, -0 and 0
                // are numerically equivalent. See
                // https://en.wikipedia.org/wiki/Signed_zero
                if (eff == 0.) eff = abs(eff);

                histRes.SetBinContent(idxTag + 1, idxTrue + 1, eff);
                histRes.SetBinError(idxTag + 1, idxTrue + 1, err);
              }
            }
            ensureUnitarity(&histRes, ptcls, debug);

            TH1D* histUnf = nullptr;
            if (histMea.Integral() > 0) {
              // Perform unfolding to get unfolded ("true") yields
              RooUnfoldResponse resp(nullptr, &histDim, &histRes);
              RooUnfoldBayes    unfoldWorker(&resp, &histMea, numOfIter);
              histUnf = static_cast<TH1D*>(unfoldWorker.Hreco());
            } else {
              // There were no tagged events, so skip the unfolding and produce
              // an "unfolded" histogram by hand with equal yields for all tags
              // to produce a uniform distribution for the prior instead
              cout << "WARNING Total tagged yield is 0. Assuming uniform "
                      "distribution."
                   << endl;
              const string key = pref + "_" + skim + "_" + year;
              count_empty_tagged[key]++;
              const string dummy_name = key + "_" + to_string(x) + "_" +
                                        to_string(y) + "_" + to_string(z);
              histUnf = new TH1D(dummy_name.c_str(), dummy_name.c_str(),
                                 totSize, 0, totSize);
              for (int idx = 0; idx < totSize; idx++) {
                histUnf->SetBinContent(idx + 1, 0.2);
              }
            }

            if (debug) {
              // Print measured and unfolded yields
              cout << "\nDEBUG Measured yields vs unfolded vs predicted ("
                   << pref << " " << x << " " << y << " " << z << ")" << endl;
              for (int idx = 0; idx != totSize; idx++) {
                double pred = 0.;
                for (int idy = 0; idy != totSize; idy++) {
                  pred += histRes.GetBinContent(idx + 1, idy + 1) *
                          histUnf->GetBinContent(idy + 1);
                }
                cout << "\t" << ptcls[idx] << ": "
                     << histMea.GetBinContent(idx + 1) << " | "
                     << histUnf->GetBinContent(idx + 1) << " | " << pred
                     << endl;
              }
              cout << endl;
            } else {
              cout << "INFO Tagged yields are ";
              for (int idx = 0; idx < totSize; idx++) {
                const int tagged = histMea.GetBinContent(idx + 1);
                cout << tagged << " (" << ptcls[idx] << ")";
                if (idx < totSize - 1) cout << ", ";
              }
              cout << endl;

              cout << "INFO Unfolded yields are " << setprecision(2);
              for (int idx = 0; idx < totSize; idx++) {
                const double unfolded = histUnf->GetBinContent(idx + 1);
                cout << unfolded << " (" << ptcls[idx] << ")";
                if (idx < totSize - 1) cout << ", ";
              }
              cout << setprecision(5) << endl;
            }

            // Save unfolded yields
            for (int idx = 0; idx != totSize; idx++) {
              auto name = getHistoName(pref, ptcls[idx], skim, year, "True");
              auto yld  = histUnf->GetBinContent(idx + 1);
              if (isnan(yld) || isinf(yld)) {
                cout << "WARNING NaN or inf detected for " << name << endl;
                yld = 0;
              }
              auto histo = histoOutGetter(name);
              if (debug) cout << "Writing unfolded yield to " << name << endl;
              histo->SetBinContent(x, y, z, yld);
            }

            // Compute unfolded ("true") probability
            auto probTrue = histoToProb(histUnf);

            for (int idxTag = 0; idxTag != totSize; idxTag++) {
              auto wtTagToMuTag    = 0.0;
              auto probTrueNormFac = 0.0;
              // Store weights for each true species
              vector<double> wtTagToMuTag_singleTrue(5, 0.);

              // Compute the shared normalization factor
              for (int idxTrue = 0; idxTrue != totSize; idxTrue++)
                probTrueNormFac +=
                    probTrue[idxTrue] *
                    histRes.GetBinContent(idxTag + 1, idxTrue + 1);

              for (int idxTrue = 0; idxTrue != totSize; idxTrue++) {
                // from pidcalib we have true -> tag
                auto probTrueToTag =
                    histRes.GetBinContent(idxTag + 1, idxTrue + 1);
                if (isnan(probTrueToTag)) probTrueToTag = 0.0;

                // probability: tag -> true
                auto probTagToTrue =
                    probTrueToTag * probTrue[idxTrue] / probTrueNormFac;
                if (isnan(probTagToTrue)) probTagToTrue = 0.0;
                auto nameTagToTrue =
                    getHistoEffName(pref, ptcls[idxTag], ptcls[idxTrue], "Tag",
                                    "True", skim, year);
                auto histoTagToTrue = histoOutGetter(nameTagToTrue);
                histoTagToTrue->SetBinContent(x, y, z, probTagToTrue);

                if (debug) {
                  cout << "  idxTrue: " << idxTrue << " idxTag: " << idxTag
                       << endl;
                  cout << "  prob = "
                       << getHistoEffName(ptcls[idxTrue], ptcls[idxTag], year)
                       << " * "
                       << getHistoName(pref, ptcls[idxTrue], skim, year)
                       << " / "
                       << "normalization" << endl;
                  cout << "       = " << probTrueToTag << " * "
                       << probTrue[idxTrue] << " / " << probTrueNormFac << " = "
                       << probTagToTrue << endl;
                }
                histInv.SetBinContent(idxTrue + 1, idxTag + 1, probTagToTrue);

                // now contract with the mu misID eff (true -> mu tag)
                auto nameTrueToMuTag =
                    getHistoEffName(ptcls[idxTrue], "Mu", year);
                auto histo = histoInGetter(nameTrueToMuTag, skim, debug);
                auto effTrueToMuTag = histo->GetBinContent(x, y, z);

                if (isnan(effTrueToMuTag)) effTrueToMuTag = 0.0;
                auto wtTagToMuTagElem = probTagToTrue * effTrueToMuTag;

                if (debug)
                  cout << "  trans. fac. = " << nameTagToTrue << " * "
                       << nameTrueToMuTag << " = " << probTagToTrue << " * "
                       << effTrueToMuTag << " = " << wtTagToMuTagElem << endl;
                wtTagToMuTag += wtTagToMuTagElem;
                wtTagToMuTag_singleTrue[idxTrue] = wtTagToMuTagElem;
              }

              // we use idxTag as the second index, which checks out
              auto name  = getHistoEffName(pref, ptcls[idxTag], "Mu", "Tag",
                                           "Tag", skim, year);
              auto histo = histoOutGetter(name);
              histo->SetBinContent(x, y, z, wtTagToMuTag);
              if (debug)
                cout << "Writing final contracted transfer factor = "
                     << wtTagToMuTag << " to " << name << endl;

              // Save transfer factors considering single true type
              for (int idxTrue = 0; idxTrue != totSize; idxTrue++) {
                auto name_single  = name + "_" + ptcls[idxTrue] + "TrueOnly";
                auto histo_single = histoOutGetter(name_single);
                histo_single->SetBinContent(x, y, z,
                                            wtTagToMuTag_singleTrue[idxTrue]);
              }
            }

            if (debug) {
              cout << "Prefix: " << pref << "; Skim: " << skim << endl;
              cout << "Bin index: x=" << x << " y=" << y << " z=" << z << endl;

              cout.precision(4);
              cout << fixed;

              cout << "The yields are (top: tag; bot: true):" << endl;
              for (int idx = 1; idx <= totSize; idx++)
                cout << setw(12) << histMea.GetBinContent(idx);
              cout << endl;
              for (int idx = 1; idx <= totSize; idx++)
                cout << setw(12) << histUnf->GetBinContent(idx);
              cout << endl;

              cout << "The prior true probabilities are:" << endl;
              for (const auto p : probTrue) cout << setw(8) << p;
              cout << endl;

              cout << "The true -> tag matrix is (row: fixed tag; "
                      "col: fixed true):"
                   << endl;
              printResVal(&histRes);

              cout << "The tag -> true efficiency matrix is (row: fixed true; "
                      "col: fixed tag):"
                   << endl;
              printResVal(&histInv);
              cout << "------" << endl;
            }
          }
        }
      }
    }
  }
}

template <typename F1, typename F2>
void unfoldDryRun(const vStr& prefix, const vStr& ptcls, const vStr& skims,
                  const vFltFlt& binnings, const vector<int>& nbins,
                  F1& histoInGetter, F2& histoOutGetter, const string& year) {
  for (const auto& pref : prefix) {
    cout << pref << ": The measured yields are stored in these histos:" << endl;
    for (const auto& pTag : ptcls) {
      for (auto skim : skims) {
        cout << skim << endl;
        auto name = getHistoName(pref, pTag, skim, year);
        histoInGetter(name, skim, true);
        cout << "  " << name << "\n" << endl;
      }
    }

    cout << pref
         << ": The unfolded yields will be stored in these histos:" << endl;
    for (const auto& pTag : ptcls) {
      for (auto skim : skims) {
        cout << skim << endl;
        auto name = getHistoName(pref, pTag, skim, year, "True");
        histoOutGetter(name);
        cout << "  " << name << "\n" << endl;
      }
    }

    cout << pref
         << ": The response matrix will be built from these histos:" << endl;
    for (const auto& pTag : ptcls) {
      for (const auto& pTrue : ptcls) {
        for (auto skim : skims) {
          cout << skim << endl;
          auto name = getHistoEffName(pTrue, pTag, year);
          histoInGetter(name, skim, true);
          cout << name << "\n";
        }
        cout << endl;
      }
      cout << endl;
    }

    cout << pref
         << ": The unfolded efficiencies will be stored in these histos:"
         << endl;
    for (const auto& pTrue : ptcls) {
      for (const auto& pTag : ptcls) {
        cout << "  ";
        for (auto skim : skims) {
          auto name =
              getHistoEffName(pref, pTag, pTrue, "Tag", "True", skim, year);
          cout << setw(38) << name;
          histoOutGetter(name);
        }
        cout << endl;
      }
      cout << endl;
    }

    cout << pref
         << ": The Mu efficiencies will be loaded from these histos:" << endl;
    for (const auto& pTrue : ptcls) {
      auto name = getHistoEffName(pTrue, "mu", year, "True", "Tag");
      cout << "  " << name << endl;
      histoOutGetter(name);
    }
  }

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
                           "Unfolding efficiency calculator (U).");

  // clang-format off
  argOpts.add_options()
    // general
    ("h,help", "print help")
    ("d,debug", "enable debug mode",
     cxxopts::value<bool>()->default_value("false"))
    ("D,dryRun", "parse config and load histos, w/o unfolding",
     cxxopts::value<bool>()->default_value("false"))
    ("Y,years", "list of sample years",
     cxxopts::value<vector<string>>()->default_value("2016,2017,2018"))
    ("s,skims", "list of skims",
     cxxopts::value<vector<string>>()->default_value("iso,1os,2os,dd,vmu,prot"))
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

  // Check info level
  const auto debug = parsedArgs["debug"].as<bool>();
  // Get YML file name
  const string ymlFile = parsedArgs["config"].as<string>();
  const string ymlName = fileNameFromPath(ymlFile);

  // parse YAML config
  const auto ymlConfig  = YAML::LoadFile(ymlFile);
  const auto ptclTarget = parsedArgs["targetParticle"].as<string>();
  const auto ptclList   = getKeyNames(ymlConfig["tags"]);
  const auto years      = parsedArgs["years"].as<vector<string>>();

  const auto [histoBinSpec, histoBinSize] = getBins(ymlConfig["binning"]);

  // Currently, only tagged yields are different between ISO+CTRL skims.
  // For VMU, no skim cuts are applied.
  const auto skims = parsedArgs["skims"].as<vector<string>>();

  // input ntuples
  auto ntpYld = make_unique<TFile>(parsedArgs["yldHisto"].as<string>().data());
  auto ntpEff = make_unique<TFile>(parsedArgs["effHisto"].as<string>().data());

  auto histoOut       = map<string, shared_ptr<TH3D>>{};
  auto histoInGetter  = getHistoInHelper(ntpYld.get(), ntpEff.get());
  auto histoOutGetter = getHistoOutHelper(&histoOut, histoBinSpec);

  // Keep stats of bins where total tagged yields are too low
  map<string, int> count_empty_tagged;

  for (const auto& year : years) {
    cout << "INFO Starting " << year << endl;
    const string year_prefix = year.substr(2, 2);
    const auto   prefix      = getKeyNames(ymlConfig["input_ntps"][year]);

    if (parsedArgs["dryRun"].as<bool>()) {
      // dry run
      unfoldDryRun(prefix, ptclList, skims, histoBinSpec, histoBinSize,
                   histoInGetter, histoOutGetter, year_prefix);
    } else {
      // unfold
      auto numOfIter = parsedArgs["iteration"].as<int>();
      unfold(prefix, ptclList, skims, histoBinSize, histoInGetter,
             histoOutGetter, year_prefix, count_empty_tagged, debug, numOfIter);
    }
  }

  if (!count_empty_tagged.empty()) {
    cout << "\nWARNING Counting bins where tagged yields are 0" << endl;
    int total = 0;
    for (const auto& [name, count] : count_empty_tagged) {
      cout << " - " << name << ": " << count << endl;
      total += count;
    }
    cout << "Total: " << total << endl;
  }

  // save output
  auto outputFilename = parsedArgs["output"].as<string>() + "/" +
                        parsedArgs["outputHisto"].as<string>();
  TFile ntpOut(outputFilename.c_str(), "RECREATE");

  for (const auto& [key, h] : histoOut) {
    if (debug) {
      cout << "DEBUG Saving objects. Key = " << key
           << " , Object name = " << h->GetName() << endl;
    }
    ntpOut.WriteObject(h.get(), key.data());
  }

  ntpOut.Close();

  return 0;
}
