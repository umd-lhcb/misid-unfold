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
  vector<double> result{};
  double         normFac = 0;
  for (int idx = 1; idx <= histo->GetNbinsX(); idx++) {
    normFac += histo->GetBinContent(idx);
  }
  // we loop twice. it's stupid but it works
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

string getHistoName(const string& prefix, const string& ptcl,
                    const string& skim, const string descr = "Tag") {
  return ""s + prefix + "__" + ptcl + descr + "__" + skim;
}

string getHistoEffName(const string& ptcl1, const string& ptcl2,
                       const string& descr1 = "True",
                       const string& descr2 = "Tag") {
  return ""s + ptcl1 + descr1 + "To" + capitalize(ptcl2) + descr2;
}

string getHistoEffName(const string& prefix, const string& ptcl1,
                       const string& ptcl2, const string& descr1,
                       const string& descr2, const string& skim) {
  const string skim_suffix = skim.size() > 0 ? "_" + skim : "";
  return ""s + prefix + "__" + getHistoEffName(ptcl1, ptcl2, descr1, descr2) +
         skim_suffix;
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

auto getHistoInHelper(TFile* ntpYld, TFile* ntpEff, TFile* ntpEffVmu) {
  auto mapHisto = make_shared<map<string, shared_ptr<TH3D>>>();

  // NOTE: We have to capture by copy, s.t. the memory is not deallocated at
  //       after the wrapper function finishes executing
  return [=](const string& histo_name, const string& skim, bool debug = false) {
    const string key = (skim == "vmu") ? histo_name + "_vmu" : histo_name;
    if (!mapHisto->count(key)) {
      TH3D* histo;
      if (key.find("__") != key.npos) {
        // this is a yld histo
        if (debug) cout << "Getting yield histo" << endl;
        histo = static_cast<TH3D*>(ntpYld->Get(histo_name.data()));
      } else {
        if (skim == "vmu") {
          if (debug)
            cout << "Getting vmu eff histo (skim = " << skim << ")" << endl;
          histo = static_cast<TH3D*>(ntpEffVmu->Get(histo_name.data()));
        } else {
          if (debug)
            cout << "Getting nominal eff histo (skim = " << skim << ")" << endl;
          histo = static_cast<TH3D*>(ntpEff->Get(histo_name.data()));
        }
      }

      if (histo == nullptr) {
        cout << "Histogram " << histo_name << " doesn't exist! terminate now..."
             << endl;
        terminate();
      }

      auto histoPtr = shared_ptr<TH3D>(histo);
      mapHisto->emplace(key, histoPtr);
      return histoPtr;
    }

    // the key already exists
    if (debug)
      cout << "Getting cached eff histo (skim = " << skim << ")" << endl;
    return mapHisto->at(key);
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

void ensureUnitarity(TH2D* res, vStr ptcls, bool debug = true) {
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
    if (debug) {
      if (abs((1 - prob) - res->GetBinContent(nbinsX, y)) > 1e-6) {
        cout << "INFO: Changing " << ptcls[y - 1] << "TrueTo"
             << ptcls[nbinsX - 1] << "Tag eff from "
             << res->GetBinContent(nbinsX, y) << " to " << 1 - prob << endl;
      }
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
            bool debug = false, int numOfIter = 4) {
  int totSize = ptcls.size();

  // These are used to stored measured yields (a vector) and response matrix (a
  // 2D matrix)
  auto histMea = new TH1D("histMea", "histMea", totSize, 0, totSize);
  auto histRes =
      new TH2D("histRes", "histRes", totSize, 0, totSize, totSize, 0, totSize);
  auto histInv = new TH2D("histInv", "histInv", totSize, 0, totSize, totSize, 0,
                          totSize);  // conceptually inverted matrix of histRes
  auto histProb = new TH1D("histProb", "histProb", totSize, 0, totSize);

  // This is used to provide dimension info for response matrix only
  auto histDim = new TH1D("histDim", "histDim", totSize, 0, totSize);
  for (int i = 1; i <= totSize; i++) histDim->SetBinContent(i, 1);

  // Main unfolding procedure
  for (int x = 1; x <= nbins[0]; x++) {
    for (int y = 1; y <= nbins[1]; y++) {
      for (int z = 1; z <= nbins[2]; z++) {
        for (auto& pref : prefix) {
          for (const auto& skim : skims) {
            // build yield vector
            for (int idx = 0; idx != totSize; idx++) {
              // For vmu case, read vmu yields but still produce (identical)
              // iso, 1os, 2os and dd weights (As in the fit, iso wil be used
              // for vmu while other skims will be ignored)
              auto name  = getHistoName(pref, ptcls[idx], skim);
              auto histo = histoInGetter(name, skim);
              histMea->SetBinContent(idx + 1, histo->GetBinContent(x, y, z));
            }

            // build response matrix (2D matrix)
            // NOTE: for a 2D array, the indexing is this:
            //         array[x][y]
            //       In our case, true -> tag translates to:
            //         true -> y index
            //         tag  -> x index
            for (int idxTag = 0; idxTag != totSize; idxTag++) {
              for (int idxTrue = 0; idxTrue != totSize; idxTrue++) {
                auto name  = getHistoEffName(ptcls[idxTrue], ptcls[idxTag]);
                auto histo = histoInGetter(name, skim);
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

                histRes->SetBinContent(idxTag + 1, idxTrue + 1, eff);
                histRes->SetBinError(idxTag + 1, idxTrue + 1, err);
              }
            }
            ensureUnitarity(histRes, ptcls);

            // perform unfolding to get unfolded ("true") yield
            RooUnfoldResponse resp(nullptr, histDim, histRes);
            RooUnfoldBayes    unfoldWorker(&resp, histMea, numOfIter);
            auto histUnf = static_cast<TH1D*>(unfoldWorker.Hreco());

            if (debug) {
              // Print measured and unfolded yields
              cout << "\nINFO: measured vs unfolded yields vs predicted ("
                   << pref << " " << x << " " << y << " " << z << ")" << endl;
              for (int idx = 0; idx != totSize; idx++) {
                double pred = 0.;
                for (int idy = 0; idy != totSize; idy++) {
                  pred += histRes->GetBinContent(idx + 1, idy + 1) *
                          histUnf->GetBinContent(idy + 1);
                }
                cout << "\t" << ptcls[idx] << ": "
                     << histMea->GetBinContent(idx + 1) << " | "
                     << histUnf->GetBinContent(idx + 1) << " | " << pred
                     << endl;
              }
              cout << endl;
            }

            // Save unfolded yields
            for (int idx = 0; idx != totSize; idx++) {
              auto name = getHistoName(pref, ptcls[idx], skim, "True");
              auto yld  = histUnf->GetBinContent(idx + 1);
              if (isnan(yld) || isinf(yld)) {
                cout << "WARNING: naN or inf detected for " << name << endl;
                yld = 0;
              }
              auto histo = histoOutGetter(name);
              cout << "Writing unfolded yield to " << name << endl;
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
                    histRes->GetBinContent(idxTag + 1, idxTrue + 1);

              for (int idxTrue = 0; idxTrue != totSize; idxTrue++) {
                // from pidcalib we have true -> tag
                auto probTrueToTag =
                    histRes->GetBinContent(idxTag + 1, idxTrue + 1);
                if (isnan(probTrueToTag)) probTrueToTag = 0.0;

                // probability: tag -> true
                auto probTagToTrue =
                    probTrueToTag * probTrue[idxTrue] / probTrueNormFac;
                if (isnan(probTagToTrue)) probTagToTrue = 0.0;
                auto nameTagToTrue = getHistoEffName(
                    pref, ptcls[idxTag], ptcls[idxTrue], "Tag", "True", skim);
                auto histoTagToTrue = histoOutGetter(nameTagToTrue);
                histoTagToTrue->SetBinContent(x, y, z, probTagToTrue);

                if (debug) {
                  cout << "  idxTrue: " << idxTrue << " idxTag: " << idxTag
                       << endl;
                  cout << "  prob = "
                       << getHistoEffName(ptcls[idxTrue], ptcls[idxTag])
                       << " * " << getHistoName(pref, ptcls[idxTrue], skim)
                       << " / "
                       << "normalization" << endl;
                  cout << "       = " << probTrueToTag << " * "
                       << probTrue[idxTrue] << " / " << probTrueNormFac << " = "
                       << probTagToTrue << endl;
                }
                histInv->SetBinContent(idxTrue + 1, idxTag + 1, probTagToTrue);

                // now contract with the mu misID eff (true -> mu tag)
                auto nameTrueToMuTag = getHistoEffName(ptcls[idxTrue], "Mu");
                auto histo           = histoInGetter(nameTrueToMuTag, skim);
                auto effTrueToMuTag  = histo->GetBinContent(x, y, z);

                if (isnan(effTrueToMuTag)) effTrueToMuTag = 0.0;
                auto wtTagToMuTagElem = probTagToTrue * effTrueToMuTag;

                if (debug)
                  cout << "  trans. fac. = " << probTagToTrue << " * "
                       << effTrueToMuTag << " = " << wtTagToMuTagElem << endl;
                wtTagToMuTag += wtTagToMuTagElem;
                wtTagToMuTag_singleTrue[idxTrue] = wtTagToMuTagElem;
              }

              // we use idxTag as the second index, which checks out
              auto name  = getHistoEffName(pref, ptcls[idxTag], "Mu", "Tag",
                                           "Tag", skim);
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
                cout << setw(12) << histMea->GetBinContent(idx);
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
              printResVal(histRes);

              cout << "The tag -> true efficiency matrix is (row: fixed true; "
                      "col: fixed tag):"
                   << endl;
              printResVal(histInv);
              cout << "------" << endl;
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
  delete histProb;
  delete histDim;
}

template <typename F1, typename F2>
void unfoldDryRun(const vStr& prefix, const vStr& ptcls, const vStr& skims,
                  const vFltFlt& binnings, const vector<int>& nbins,
                  F1& histoInGetter, F2& histoOutGetter) {
  for (const auto& pref : prefix) {
    cout << pref << ": The measured yields are stored in these histos:" << endl;
    for (const auto& pTag : ptcls) {
      for (auto skim : skims) {
        cout << skim << endl;
        auto name = getHistoName(pref, pTag, skim);
        histoInGetter(name, skim, true);
        cout << "  " << name << "\n" << endl;
      }
    }

    cout << pref
         << ": The unfolded yields will be stored in these histos:" << endl;
    for (const auto& pTag : ptcls) {
      for (auto skim : skims) {
        cout << skim << endl;
        auto name = getHistoName(pref, pTag, skim, "True");
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
          auto name = getHistoEffName(pTrue, pTag);
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
          auto name = getHistoEffName(pref, pTag, pTrue, "Tag", "True", skim);
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
      auto name = getHistoEffName(pTrue, "mu", "True", "Tag");
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
                           "unfolding efficiency calculator (U).");

  // clang-format off
  argOpts.add_options()
    // general
    ("h,help", "print help")
    ("d,debug", "enable debug mode",
     cxxopts::value<bool>()->default_value("false"))
    ("D,dryRun", "parse config and load histos, w/o unfolding",
     cxxopts::value<bool>()->default_value("false"))
    ("Y,year", "sample year",
     cxxopts::value<string>()->default_value("2016"))
    // input/output
    ("e,effHisto", "specify input ntuple containing efficiency histos",
     cxxopts::value<string>())
    ("v,effHistoVmu", "specify input ntuple containing vmu efficiency histos",
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
  auto debug = parsedArgs["debug"].as<bool>();
  // Get YML file name
  const string ymlFile = parsedArgs["config"].as<string>();
  const string ymlName = fileNameFromPath(ymlFile);

  // parse YAML config
  auto ymlConfig  = YAML::LoadFile(ymlFile);
  auto ptclTarget = parsedArgs["targetParticle"].as<string>();
  auto ptclList   = getKeyNames(ymlConfig["tags"]);
  auto year       = parsedArgs["year"].as<string>();
  auto prefix     = getKeyNames(ymlConfig["input_ntps"][year]);
  auto [histoBinSpec, histoBinSize] = getBins(ymlConfig["binning"]);

  // input ntuples
  auto ntpYld = make_unique<TFile>(parsedArgs["yldHisto"].as<string>().data());
  auto ntpEff = make_unique<TFile>(parsedArgs["effHisto"].as<string>().data());
  auto ntpEffVmu =
      make_unique<TFile>(parsedArgs["effHistoVmu"].as<string>().data());

  auto histoOut = map<string, shared_ptr<TH3D>>{};
  auto histoInGetter =
      getHistoInHelper(ntpYld.get(), ntpEff.get(), ntpEffVmu.get());
  auto histoOutGetter = getHistoOutHelper(&histoOut, histoBinSpec);

  // Produce list of skims. Currently, only tagged yields are different.
  // For vmu, no skim cuts are applied.
  const vStr skims = {"iso", "1os", "2os", "dd", "vmu"};

  // dry run
  if (parsedArgs["dryRun"].as<bool>()) {
    unfoldDryRun(prefix, ptclList, skims, histoBinSpec, histoBinSize,
                 histoInGetter, histoOutGetter);
    return 0;
  }

  // unfold
  auto numOfIter = parsedArgs["iteration"].as<int>();
  unfold(prefix, ptclList, skims, histoBinSize, histoInGetter, histoOutGetter,
         debug, numOfIter);

  // save output
  auto outputFilename = parsedArgs["output"].as<string>() + "/" +
                        parsedArgs["outputHisto"].as<string>();
  auto ntpOut = make_unique<TFile>(outputFilename.data(), "RECREATE");

  for (const auto& [key, h] : histoOut)
    ntpOut->WriteObject(h.get(), key.data());

  return 0;
}
