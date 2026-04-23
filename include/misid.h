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

constexpr double RICH1_Z = 2500.;  // mm
constexpr double T3_Z    = 9000.;  // mm

const unordered_map<string, int> year_idx{
    {"2016", 0}, {"2017", 1}, {"2018", 2}};

enum Sample { ISO_CTRL, VMU, FAKE_MU };

const map<Sample, string> SAMPLES = {
    {ISO_CTRL, "iso_ctrl"}, {VMU, "vmu"}, {FAKE_MU, "fake_mu"}};

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

// Get upper d0_M fit limit according to kinematical bin.
// The PIDCalib range extends to 1.910 GeV, but we may have to exclude some of
// it due to distortions introduced by the WM cut, which are easily accounted
// for in the MC-derived PDFs but not in the combinatorial shapes.
double get_d0m_upper_limit(const TString &probe, const int &p_idx,
                           const int &eta_idx) {
  if (probe == "k") {
    if (p_idx == 5) {
      return 1.900;
    }
    if (p_idx == 4) {
      return 1.902;
    }
    if (p_idx == 3 && eta_idx == 0) {
      return 1.906;
    }
  } else if (probe == "pi") {
    if (p_idx == 2) {
      return 1.906;
    }
    if (p_idx == 1) {
      return 1.904;
    }
    if (p_idx == 0) {
      return 1.902;
    }
  }
  return 1.910;
}

//////////////////////
// Helper functions //
//////////////////////

// Check if candidate satisfies muon PID requirements.
bool check_mu_pid_run1rdx(const bool &probe_is_muon, const double &probe_DLLmu,
                          const double &probe_DLLe, const float &probe_ubdt,
                          const Sample &sample) {
  switch (sample) {
    case ISO_CTRL:
      return probe_is_muon && (probe_DLLe < 1.0) && (probe_DLLmu > 2.0) &&
             (probe_ubdt > 0.25);
    case VMU:
      return probe_is_muon && (probe_DLLe < 1.0) && (probe_DLLmu > 2.0) &&
             (probe_ubdt < 0.25);
    case FAKE_MU:
      // IMPORTANT: The FAKE_MU pid requirement is !isMuon, but here we flip it
      // and calculate the complementary efficiency so that the K/pi misid case
      // (with less stats) falls in the "passed" sub-sample. The proper
      // efficiency is calculated afterwards as 1.0 - complementary_eff.
      return probe_is_muon;
    default:
      cout << "ERROR Unexpected <sample> argument in check_mu_pid()" << endl;
      return false;
  }
}

bool check_mu_pid_run2ang(const bool &probe_is_muon, const double &probe_DLLe,
                          const float &probe_ubdt, const Sample &sample) {
  switch (sample) {
    case ISO_CTRL:
      return probe_is_muon && (probe_DLLe < 1.0) && (probe_ubdt > 0.25);
    case VMU:
      return probe_is_muon && (probe_DLLe < 1.0) && (probe_ubdt < 0.25);
    case FAKE_MU:
      // IMPORTANT: The FAKE_MU pid requirement is !isMuon, but here we flip it
      // and calculate the complementary efficiency so that the K/pi misid case
      // (with less stats) falls in the "passed" sub-sample. The proper
      // efficiency is calculated afterwards as 1.0 - complementary_eff.
      return probe_is_muon;
    default:
      cout << "ERROR Unexpected <sample> argument in check_mu_pid()" << endl;
      return false;
  }
}

bool check_mu_pid(const bool &probe_is_muon, const double &probe_DLLmu,
                  const double &probe_DLLe, const float &probe_ubdt,
                  const Sample &sample, const bool &run2ang) {
  if (run2ang) {
    return check_mu_pid_run2ang(probe_is_muon, probe_DLLe, probe_ubdt, sample);
  } else {
    return check_mu_pid_run1rdx(probe_is_muon, probe_DLLmu, probe_DLLe,
                                probe_ubdt, sample);
  }
}

// Check if track is matched to a h -> mu nu decay-in-fight chain.
bool check_dif(const int &probe_trueid, const int &probe_daughter0_trueid,
               const int &probe_daughter1_trueid, const int &probe_mc_mom_nd,
               const string &probe) {
  // Deduce hadron ID
  int h_id = 0;
  if (probe == "k") {
    h_id = K_ID;
  } else if (probe == "pi") {
    h_id = PI_ID;
  } else {
    return false;
  }

  // First case: track matched to hadron
  // DiF particles will have mu, nu_mu as daughters
  if (std::abs(probe_trueid) == h_id) {
    return (std::abs(probe_daughter0_trueid) == MU_ID) &&
           (std::abs(probe_daughter1_trueid) == NU_MU_ID);
  }

  // Second case: track matched to muon
  // Can't check if sister particle is a muon neutrino, so just check if DiF
  // vertex produced two particles
  if (std::abs(probe_trueid) == MU_ID) {
    return probe_mc_mom_nd == 2;
  }

  return false;
}

// Comparison function for PDG IDs.
// Sort by absolute value, and if it's the same sort by sign.
bool comp_pdg_ids(const int &a, const int &b) {
  if (std::abs(a) != std::abs(b)) {
    return std::abs(a) > std::abs(b);
  } else {
    return a > b;
  }
}

// Split DiF tracks according to z coordinate of DiF vertex.
// This determines whether the hadron or the muon was responsible for the RICH1
// response.
bool after_rich1(const double &dif_vtx_z) { return dif_vtx_z > RICH1_Z; }
// This determines if DiF happened before T3. If not, there is no
// resolution impact on the momentum measurement.
bool before_t3(const double &dif_vtx_z) { return dif_vtx_z < T3_Z; }

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

TLine *draw_line(const RooPlot *rp, const double &x, const EColor &color = kRed,
                 const ELineStyle &style = kDashed) {
  TLine *line = new TLine(x, rp->GetMinimum(), x, rp->GetMaximum());
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->Draw();
  return line;
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
