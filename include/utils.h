// Author: Yipeng Sun, Svede Braun
// License: BSD 2-clause
// Last Change: Sun Apr 24, 2022 at 08:23 PM -0400

#pragma once

#include <filesystem>
#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TChainElement.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TString.h"
#include <TTree.h>

#include "RooDataSet.h"
#include "RooRealVar.h"

using std::string, std::cout, std::endl;

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
  const int last = path.find_last_of(".");
  return path.substr(first+1, last-first-1);
}

// Convert an integer number of seconds to the format Xh Ym Zs
TString format_time(const int &duration) {
  const double d = duration / 3600.;

  const int h = d;
  const int m = (d - h) * 60;
  const int s = ((d - h) * 60 - m) * 60;

  return TString::Format("%dh %dm %ds", h, m, s);
}

// Calculate square root of sum of squares
double sqrt_sum_sq(const double& x1, const double& x2) {
  return sqrt(x1*x1 + x2*x2);
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
void print_files(const TChain &ch) {
  TObjArray *fileElements = ch.GetListOfFiles();
  for (TObject *op : *fileElements) {
    auto chainElement = static_cast<TChainElement *>(op);
    cout << "  - " << chainElement->GetTitle() << endl;
  }
}

// Check if <val> is inside <var>'s default range
bool in_var_range(const RooRealVar &var, const double &val) {
  return (val <= var.getMax()) && (val >= var.getMin());
}

// Check if <val> is between <down> and <up>
bool in_range(const double &down, const double &val, const double &up) {
  return (val <= up) && (val >= down);
}
