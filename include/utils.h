// Author: Yipeng Sun, Svede Braun
// License: BSD 2-clause
// Last Change: Sun Apr 24, 2022 at 08:23 PM -0400

#pragma once

#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>

#include "TCanvas.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TString.h"
#include "TTree.h"

#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"

using std::unique_ptr, std::string, std::cout, std::endl;

using namespace RooFit;

/////////////
// General //
/////////////

// We copy the string when passing it to the function, so the original one
// remains unchanged
string capitalize(string str) {
  for (auto& s : str) {
    s = toupper(s);
    break;
  }
  return str;
}

string absDirPath(string pathRaw) {
  auto path    = std::filesystem::path(pathRaw);
  auto dirPath = path.parent_path();
  return std::filesystem::absolute(dirPath).string();
}

string fileNameFromPath(const string& path) {
  const int first = path.find_last_of("/\\");
  const int last  = path.find_last_of(".");
  return path.substr(first + 1, last - first - 1);
}

// Convert an integer number of seconds to the format Xh Ym Zs
TString format_time(const int& duration) {
  const double d = duration / 3600.;

  const int h = d;
  const int m = (d - h) * 60;
  const int s = ((d - h) * 60 - m) * 60;

  return TString::Format("%dh %dm %ds", h, m, s);
}

// Calculate square root of sum of squares
double sqrt_sum_sq(const double& x1, const double& x2) {
  return std::sqrt(x1 * x1 + x2 * x2);
}

double sqrt_sum_sq(const double& x1, const double& x2, const double& x3) {
  return std::sqrt(x1 * x1 + x2 * x2 + x3 * x3);
}

//////////////////
// ROOT-related //
//////////////////

bool branchExists(TTree* tree, string brName) {
  auto brPtr = tree->GetListOfBranches()->FindObject(brName.data());

  if (brPtr == nullptr) return false;
  return true;
}

// Print list of files added to <ch>
void print_files(const TChain& ch) {
  TObjArray* fileElements = ch.GetListOfFiles();
  for (TObject* op : *fileElements) {
    auto chainElement = static_cast<TChainElement*>(op);
    cout << "  - " << chainElement->GetTitle() << endl;
  }
}

// Check if <val> is inside <var>'s default range
bool in_var_range(const RooRealVar& var, const double& val) {
  return (val <= var.getMax()) && (val >= var.getMin());
}

// Check if <val> is between <down> and <up>
bool in_range(const double& down, const double& val, const double& up) {
  return (val <= up) && (val >= down);
}

// Create 2D pulls between RooFit pdf and data TH2F
TH2F* make_2d_pulls(const RooAbsPdf& pdf, TH2F* data, RooRealVar& x,
                    RooRealVar& y) {
  TH2F* pull = new TH2F();
  data->Copy(*pull);

  pull->SetMaximum(5.5);
  pull->SetMinimum(-5.5);
  pull->SetContour(11);

  data->SetBinErrorOption(TH1::kPoisson);

  pull->SetStats(false);

  const double x_max = data->GetXaxis()->GetXmax();
  const double x_min = data->GetXaxis()->GetXmin();

  const double y_max = data->GetYaxis()->GetXmax();
  const double y_min = data->GetYaxis()->GetXmin();

  x.setRange("hist_range", x_min, x_max);
  y.setRange("hist_range", y_min, y_max);

  unique_ptr<RooAbsReal> integral_full(pdf.createIntegral(
      RooArgSet(x, y), NormSet(RooArgSet(x, y)), Range("hist_range")));

  const double int_full_val = integral_full->getVal();
  const double n_entries    = data->Integral();

  for (int x_idx = 1; x_idx <= data->GetNbinsX(); x_idx++) {
    for (int y_idx = 1; y_idx <= data->GetNbinsY(); y_idx++) {
      const double bin_x_max = data->GetXaxis()->GetBinUpEdge(x_idx);
      const double bin_x_min = data->GetXaxis()->GetBinLowEdge(x_idx);

      const double bin_y_max = data->GetYaxis()->GetBinUpEdge(y_idx);
      const double bin_y_min = data->GetYaxis()->GetBinLowEdge(y_idx);

      x.setRange("bin_range", bin_x_min, bin_x_max);
      y.setRange("bin_range", bin_y_min, bin_y_max);

      unique_ptr<RooAbsReal> integral(pdf.createIntegral(
          RooArgSet(x, y), NormSet(RooArgSet(x, y)), Range("bin_range")));

      const double int_bin_val = integral->getVal();

      const double bin_res = data->GetBinContent(x_idx, y_idx) -
                             int_bin_val / int_full_val * n_entries;
      const double bin_error = (bin_res < 0)
                                   ? data->GetBinErrorUp(x_idx, y_idx)
                                   : data->GetBinErrorLow(x_idx, y_idx);

      if (bin_res == 0. || bin_error == 0. || std::isnan(bin_res) ||
          std::isnan(bin_error) || std::isinf(bin_res) ||
          std::isinf(bin_error)) {
        cout << "WARNING Problem in make_2d_pulls with " << data->GetName()
             << endl;
        cout << "  - x_idx = " << x_idx << endl;
        cout << "  - y_idx = " << y_idx << endl;
        cout << "  - bin_x_max = " << bin_x_max << endl;
        cout << "  - bin_x_min = " << bin_x_min << endl;
        cout << "  - bin_y_max = " << bin_y_max << endl;
        cout << "  - bin_y_min = " << bin_y_min << endl;
        cout << "  - x_max = " << x_max << endl;
        cout << "  - x_min = " << x_min << endl;
        cout << "  - y_max = " << y_max << endl;
        cout << "  - y_min = " << y_min << endl;
        cout << "  - int_bin_val = " << int_bin_val << endl;
        cout << "  - bin_res = " << bin_res << endl;
        cout << "  - bin_error = " << bin_error << endl;
        cout << "  - n_entries = " << n_entries << endl;
        cout << "  - int_full_val = " << int_full_val << endl;
      }

      pull->SetBinContent(x_idx, y_idx, bin_res / bin_error);
    }
  }

  return pull;
}
