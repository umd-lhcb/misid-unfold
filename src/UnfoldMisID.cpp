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
                    string descr = "Tag") {
  return ""s + prefix + "__" + ptcl + descr;
}

string getHistoEffName(const string& ptcl1, const string& ptcl2,
                       const string& descr1 = "True",
                       const string& descr2 = "Tag") {
  return ""s + ptcl1 + descr1 + "To" + capitalize(ptcl2) + descr2;
}

string getHistoEffName(const string& prefix, const string& ptcl1,
                       const string& ptcl2, const string& descr1,
                       const string& descr2) {
  return ""s + prefix + "__" + getHistoEffName(ptcl1, ptcl2, descr1, descr2);
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

  for (int idxRow = 1; idxRow <= nbinsX; idxRow++) {
    for (int idxCol = 1; idxCol <= nbinsY; idxCol++) {
      auto elem = getter(res, idxRow, idxCol);
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
  return [=](string key) {
    if (!mapHisto->count(key)) {
      TH3D* histo;
      if (key.find("__") != key.npos)
        // this is a yld histo
        histo = static_cast<TH3D*>(ntpYld->Get(key.data()));
      else
        histo = static_cast<TH3D*>(ntpEff->Get(key.data()));

      if (histo == nullptr) {
        cout << "Histogram " << key << " doesn't exist! terminate now..."
             << endl;
        terminate();
      }

      auto histoPtr = shared_ptr<TH3D>(histo);
      mapHisto->emplace(key, histoPtr);
      return histoPtr;
    }

    // the key already exists
    return mapHisto->at(key);
  };
}

auto getHistoOutHelper(map<string, shared_ptr<TH3D>>* mapHisto,
                       vFltFlt&                       binnings) {
  return [=](string key) {
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

void ensureUnitarity(TH2D* res, bool debug = true) {
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
    res->SetBinContent(nbinsX, y, 1 - prob);
  }

  if (debug) {
    cout << "The fixed true -> tag matrix is (row: fixed tag; col: fixed true):"
         << endl;
    printResVal(res);
  }
}

template <typename F1, typename F2>
void unfold(vStr& prefix, vStr& ptcls, vector<int>& nbins, F1& histoInGetter,
            F2& histoOutGetter, bool debug = false, int numOfIter = 4) {
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
          // build yield vector
          for (int idx = 0; idx != totSize; idx++) {
            auto name  = getHistoName(pref, ptcls[idx]);
            auto histo = histoInGetter(name);
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
              auto histo = histoInGetter(name);
              auto eff   = histo->GetBinContent(x, y, z);
              auto err   = histo->GetBinError(x, y, z);

              if (debug)
                cout << "  Loading efficiency from " << name << ", got " << eff
                     << endl;
              if (isnan(eff) || isinf(eff)) eff = 0.0;

              histRes->SetBinContent(idxTag + 1, idxTrue + 1, abs(eff));
              histRes->SetBinError(idxTag + 1, idxTrue + 1, err);
            }
          }
          ensureUnitarity(histRes);

          // perform unfolding to get unfolded ("true") yield
          RooUnfoldResponse resp(nullptr, histDim, histRes);
          RooUnfoldBayes    unfoldWorker(&resp, histMea, numOfIter);
          auto              histUnf = static_cast<TH1D*>(unfoldWorker.Hreco());

          // Save unfolded yields
          for (int idx = 0; idx != totSize; idx++) {
            auto name = getHistoName(pref, ptcls[idx], "True");
            auto yld  = histUnf->GetBinContent(idx + 1);
            if (isnan(yld) || isinf(yld)) {
              cout << "WARNING: naN or inf detected for " << name << endl;
              yld = 0;
            }
            auto histo = histoOutGetter(name);
            cout << "Writing unfolded yield to " << name << endl;
            ;
            histo->SetBinContent(x, y, z, yld);
          }

          // Compute unfolded ("true") probability
          auto probTrue = histoToProb(histUnf);

          for (int idxTag = 0; idxTag != totSize; idxTag++) {
            auto wtTagToMuTag    = 0.0;
            auto probTrueNormFac = 0.0;

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
                  pref, ptcls[idxTag], ptcls[idxTrue], "Tag", "True");
              auto histoTagToTrue = histoOutGetter(nameTagToTrue);
              histoTagToTrue->SetBinContent(x, y, z, probTagToTrue);

              if (debug) {
                cout << "  idxTrue: " << idxTrue << " idxTag: " << idxTag
                     << endl;
                cout << "  prob = "
                     << getHistoEffName(ptcls[idxTrue], ptcls[idxTag]) << " * "
                     << getHistoName(pref, ptcls[idxTrue]) << " / "
                     << "normalization" << endl;
                cout << "       = " << probTrueToTag << " * "
                     << probTrue[idxTrue] << " / " << probTrueNormFac << " = "
                     << probTagToTrue << endl;
              }
              histInv->SetBinContent(idxTrue + 1, idxTag + 1, probTagToTrue);

              // now contract with the mu misID eff (true -> mu tag)
              auto nameTrueToMuTag = getHistoEffName(ptcls[idxTrue], "Mu");
              auto histo           = histoInGetter(nameTrueToMuTag);
              auto effTrueToMuTag  = histo->GetBinContent(x, y, z);

              if (isnan(effTrueToMuTag)) effTrueToMuTag = 0.0;
              auto wtTagToMuTagElem = probTagToTrue * effTrueToMuTag;

              if (debug)
                cout << "  trans. fac. = " << probTagToTrue << " * "
                     << effTrueToMuTag << " = " << wtTagToMuTagElem << endl;
              wtTagToMuTag += wtTagToMuTagElem;
            }

            // we use idxTag as the second index, which checks out
            auto name =
                getHistoEffName(pref, ptcls[idxTag], "Mu", "Tag", "Tag");
            auto histo = histoOutGetter(name);
            histo->SetBinContent(x, y, z, wtTagToMuTag);
            if (debug)
              cout << "Writing final contracted transfer factor = "
                   << wtTagToMuTag << " to " << name << endl;
          }

          if (debug) {
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

  // cleanups
  delete histMea;
  delete histRes;
  delete histInv;
  delete histProb;
  delete histDim;
}

template <typename F1, typename F2>
void unfoldDryRun(vStr& prefix, vStr& ptcls, vFltFlt& binnings,
                  vector<int>& nbins, F1& histoInGetter, F2& histoOutGetter) {
  for (const auto& pref : prefix) {
    cout << pref << ": The measured yields are stored in these histos:" << endl;
    for (const auto& pTag : ptcls) {
      auto name = getHistoName(pref, pTag);
      cout << "  " << name << endl;
      histoInGetter(name);
    }

    cout << pref
         << ": The unfolded yields will be stored in these histos:" << endl;
    for (const auto& pTag : ptcls) {
      auto name = getHistoName(pref, pTag, "True");
      cout << "  " << name << endl;
      histoOutGetter(name);
    }

    cout << pref
         << ": The response matrix will be built from these histos:" << endl;
    for (const auto& pTag : ptcls) {
      cout << "  ";
      for (const auto& pTrue : ptcls) {
        auto name = getHistoEffName(pTrue, pTag);
        cout << setw(16) << name;
        histoOutGetter(name);
      }
      cout << endl;
    }

    cout << pref
         << ": The unfolded efficiencies will be stored in these histos:"
         << endl;
    for (const auto& pTrue : ptcls) {
      cout << "  ";
      for (const auto& pTag : ptcls) {
        auto name = getHistoEffName(pTag, pTrue);
        cout << setw(16) << name;
        histoOutGetter(name);
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
  auto year       = parsedArgs["year"].as<string>();
  auto prefix     = getKeyNames(ymlConfig["input_ntps"][year]);
  auto [histoBinSpec, histoBinSize] = getBins(ymlConfig["binning"]);

  // input ntuples
  auto ntpYld = make_unique<TFile>(parsedArgs["yldHisto"].as<string>().data());
  auto ntpEff = make_unique<TFile>(parsedArgs["effHisto"].as<string>().data());

  auto histoOut       = map<string, shared_ptr<TH3D>>{};
  auto histoInGetter  = getHistoInHelper(ntpYld.get(), ntpEff.get());
  auto histoOutGetter = getHistoOutHelper(&histoOut, histoBinSpec);

  // dry run
  if (parsedArgs["dryRun"].as<bool>()) {
    unfoldDryRun(prefix, ptclList, histoBinSpec, histoBinSize, histoInGetter,
                 histoOutGetter);
    return 0;
  }

  // unfold
  auto debug     = parsedArgs["debug"].as<bool>();
  auto numOfIter = parsedArgs["iteration"].as<int>();
  unfold(prefix, ptclList, histoBinSize, histoInGetter, histoOutGetter, debug,
         numOfIter);

  // save output
  auto outputFilename = parsedArgs["output"].as<string>() + "/" +
                        parsedArgs["outputHisto"].as<string>();
  auto ntpOut = make_unique<TFile>(outputFilename.data(), "RECREATE");

  for (const auto& [key, h] : histoOut)
    ntpOut->WriteObject(h.get(), key.data());

  return 0;
}
