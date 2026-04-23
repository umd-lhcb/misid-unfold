// Author: Lucas Meyer Garcia
// License: BSD 2-clause
//
// Description: Constants and helper functions for d0BkgDecays.cpp and
// GetMisIDCorrections.cpp.

#pragma once

#include <cmath>

#include "TCanvas.h"
#include "TColor.h"
#include "TH1.h"
#include "TLine.h"

#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"

using std::cout, std::endl;
using std::string, std::unique_ptr, std::unordered_map, std::map, std::vector,
    std::array, std::pair;

using namespace RooFit;

///////////////
// Constants //
///////////////

constexpr int Dst_ID   = 413;
constexpr int D0_ID    = 421;
constexpr int PI_ID    = 211;
constexpr int K_ID     = 321;
constexpr int MU_ID    = 13;
constexpr int NU_MU_ID = 14;
constexpr int KS_ID    = 310;

// PDG 2025
constexpr double DM   = 0.1454258;   // GeV
constexpr double D0_M = 1.86484;     // GeV
constexpr double PI_M = 0.13957039;  // GeV
constexpr double K_M  = 0.493677;    // GeV

// PIDCalib limis
constexpr double D0_M_min = 1.825;  // GeV
constexpr double D0_M_max = 1.910;  // GeV
constexpr double DM_min   = 0.141;  // GeV
constexpr double DM_max   = 0.153;  // GeV

constexpr double ONE_SIGMA = 0.682689492137;

const unordered_map<string, int> year_idx{
    {"2016", 0}, {"2017", 1}, {"2018", 2}};

// Fit variables labels
const TString d0_m_label = "#font[12]{m}(#font[12]{D}^{0})";
const TString dm_label   = "#font[12]{m}(#font[12]{D*}^{+}) - " + d0_m_label;
const TString d0_m_x_dm_label = d0_m_label + ";" + dm_label;

/////////////
// Binning //
/////////////

constexpr double BINS_P[] = {3e+3, 6e+3, 10e+3, 15.6e+3, 27e+3, 60e+3, 100e+3};
constexpr double BINS_ETA[]     = {1.7, 3.6, 5.0};
constexpr double BINS_NTRACKS[] = {0., 200., 600.};

constexpr int N_BINS_P       = sizeof(BINS_P) / sizeof(double) - 1;
constexpr int N_BINS_ETA     = sizeof(BINS_ETA) / sizeof(double) - 1;
constexpr int N_BINS_NTRACKS = sizeof(BINS_NTRACKS) / sizeof(double) - 1;

//////////////////////
// Helper functions //
//////////////////////

// Comparison function for PDG IDs.
// Sort by absolute value, and if it's the same sort by sign.
bool comp_pdg_ids(const int &a, const int &b) {
  if (std::abs(a) != std::abs(b)) {
    return std::abs(a) > std::abs(b);
  } else {
    return a > b;
  }
}

// Calculate weighted average of elements {x} with weights 1/unc^2.
template <size_t N>
constexpr double weighted_average(const array<double, N> &x,
                                  const array<double, N> &uncs) {
  double sum   = 0.;
  double sum_w = 0.;
  for (unsigned i = 0; i < x.size(); i++) {
    const double w = 1. / (uncs[i] * uncs[i]);
    sum += x[i] * w;
    sum_w += w;
  }
  return sum / sum_w;
}

// Check if fit has converged and covariant matrix is ok.
bool fit_ok(const RooFitResult *fit) {
  // 4000 probably means IMPROVE failed to find a new minimum, which is ok.
  // See section "Access to the fit status" in
  // https://root.cern.ch/doc/v624/classTH1.html#a7e7d34c91d5ebab4fc9bba3ca47dabdd
  const bool fit_converged = (fit->status() == 4000) || (fit->status() == 0);

  // Covariance matrix status (https://root.cern.ch/download/minuit.pdf)
  // - 0 Not calculated at all
  // - 1 Diagonal approximation only, not accurate
  // - 2 Full matrix, but forced positive-definite
  // - 3 Full accurate covariance matrix
  const bool cov_matrix_ok = (fit->covQual() == 3);

  return fit_converged && cov_matrix_ok;
}

// Add contents of <ds_source> to <ds_target>.
void append_to_dataset(const RooDataSet &ds_source, RooDataSet &ds_target,
                       const RooRealVar &w_var) {
  for (int i = 0; i < ds_source.numEntries(); ++i) {
    const RooArgSet  &argset   = *ds_source.get(i);
    const RooRealVar &d0_m_var = static_cast<RooRealVar &>(argset["d0_m_var"]);
    const RooRealVar &dm_var   = static_cast<RooRealVar &>(argset["dm_var"]);
    ds_target.addFast(RooArgSet(d0_m_var, dm_var, w_var), w_var.getVal());
  }
}

// Plots histograms for all observables is <ds> and saves them at <path>
void plot_dataset(const RooDataSet &ds, const TString path) {
  if (ds.numEntries() == 0) {
    cout << "\nERROR Empty dataset. Cannot print distributions for "
         << ds.GetName() << endl;
    return;
  }

  TCanvas c("c", "c", 1280, 960);
  TString hist_prefix = "h_";

  const RooArgSet *vars = ds.get(0);
  for (auto arg : *vars) {
    RooRealVar *var = static_cast<RooRealVar *>(arg);

    // Set histogram range a little beyond the roorealvar limits to avoid
    // hiding parameters at limit in overflow bin
    const double var_max = var->getMax(), var_min = var->getMin();
    const double range     = var_max - var_min;
    const double range_max = var_max + 0.001 * range;
    const double range_min = var_min - 0.001 * range;

    TH1D *hist = static_cast<TH1D *>(
        ds.createHistogram(hist_prefix + var->GetName(), *var,
                           Binning(100, range_min, range_max)));
    hist->Draw();

    c.SaveAs(path + var->GetName() + ".pdf");
    delete hist;
  }
}

///////////////////
// Color palette //
///////////////////

// See https://jfly.uni-koeln.de/color/.
const TColor col1(0.00f, 0.00f, 0.00f);  // Black
const TColor col2(0.90f, 0.60f, 0.00f);  // Orange
const TColor col3(0.35f, 0.70f, 0.90f);  // Sky blue
const TColor col4(0.00f, 0.60f, 0.50f);  // Bluish green
const TColor col5(0.95f, 0.90f, 0.25f);  // Yellow
const TColor col6(0.00f, 0.45f, 0.70f);  // Blue
const TColor col7(0.80f, 0.40f, 0.00f);  // Vermillion
const TColor col8(0.80f, 0.60f, 0.70f);  // Reddish purple
