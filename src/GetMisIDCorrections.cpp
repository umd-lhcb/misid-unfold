// Author: Lucas Meyer Garcia
// License: BSD 2-clause
//
// Description: Calculate misid efficiencies for kaons and pions

#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TString.h"

#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooBinning.h"
#include "RooCategory.h"
#include "RooChi2Var.h"
#include "RooConstVar.h"
#include "RooCrystalBall.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooLinkedList.h"
#include "RooPlot.h"
#include "RooPowerLaw.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"

#include <yaml-cpp/yaml.h>
#include <cxxopts.hpp>

using namespace RooFit;
using namespace std::chrono;

using std::cout, std::endl;
using std::string, std::unique_ptr, std::unordered_map, std::vector;

constexpr int    Dst_ID    = 413;
constexpr double DM_min    = 0.141;      // GeV
constexpr double DM_max    = 0.153;      // GeV
constexpr double DM        = 0.1454258;  // GeV
constexpr int    D0_ID     = 421;
constexpr double D0_M_min  = 1.825;    // GeV
constexpr double D0_M_max  = 1.910;    // GeV
constexpr double D0_M      = 1.86484;  // GeV
constexpr int    PI_ID     = 211;
constexpr double PI_M      = 0.13957039;  // GeV
constexpr int    K_ID      = 321;
constexpr int    MU_ID     = 13;
constexpr double ONE_SIGMA = 0.682689492137;

constexpr double BINS_P[] = {3e+3, 6e+3, 10e+3, 15.6e+3, 27e+3, 60e+3, 100e+3};
constexpr double BINS_ETA[]     = {1.7, 3.6, 5.0};
constexpr double BINS_NTRACKS[] = {0., 200., 600.};

constexpr int N_BINS_P       = sizeof(BINS_P) / sizeof(double) - 1;
constexpr int N_BINS_ETA     = sizeof(BINS_ETA) / sizeof(double) - 1;
constexpr int N_BINS_NTRACKS = sizeof(BINS_NTRACKS) / sizeof(double) - 1;

const unordered_map<string, int> year_idx{
    {"2016", 0}, {"2017", 1}, {"2018", 2}};

void print_files(const TChain &ch) {
  TObjArray *fileElements = ch.GetListOfFiles();
  for (TObject *op : *fileElements) {
    auto chainElement = static_cast<TChainElement *>(op);
    cout << "  - " << chainElement->GetTitle() << endl;
  }
}

TString format_time(const int &duration) {
  const double d = duration / 3600.;

  const int h = d;
  const int m = (d - h) * 60;
  const int s = ((d - h) * 60 - m) * 60;

  return TString::Format("%dh %dm %ds", h, m, s);
}

bool in_range(const RooRealVar &var, const double &val) {
  return (val <= var.getMax()) && (val >= var.getMin());
}

bool fit_ok(const RooFitResult *fit) {
  // 4000 Probably means IMPROVE failed to find a new minimum, which is ok.
  // See section "Access to the fit status" in
  // https://root.cern.ch/doc/v624/classTH1.html#a7e7d34c91d5ebab4fc9bba3ca47dabdd
  const bool fit_converged = (fit->status() == 4000) || (fit->status() == 0);

  // Covariance matrix status (https://root.cern.ch/download/minuit.pdf)
  // - 0 Not calculated at all
  // - 1 Diagonal approximation only, not accurate
  // - 2 Full matrix, but forced positive-definite
  // - 3 Full accurate covariance matrix
  const bool cov_matrix_ok = fit->covQual() == 3;

  return fit_converged && cov_matrix_ok;
}

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

int main(int argc, char **argv) {
  auto start = high_resolution_clock::now();

  cxxopts::Options argOpts(
      "GetMisIDCorrections",
      "Calculate misid corrections (CrystalBall core -> core + tail).");

  // clang-format off
  argOpts.add_options()
    ("h,help", "Print help")
    ("d,debug", "Enable debug mode",
     cxxopts::value<bool>()->default_value("false"))
    ("p,particles", "Specify probed particle",
     cxxopts::value<vector<string>>()->default_value("k,pi"))
    ("c,config", "Specify input YAML config file",
     cxxopts::value<string>())
    ("years", "Specify input YAML config file",
     cxxopts::value<vector<string>>()->default_value("2016,2017,2018"))
    ("o,output", "Specify output folder",
     cxxopts::value<string>()->default_value("gen/"))
    ("s,simultaneous", "Perform full simultaneous fits of data+MC",
     cxxopts::value<bool>()->default_value("false"))
    ("n,dry-run", "Do not perform fits",
     cxxopts::value<bool>()->default_value("false"))
    ("m,mc-only", "Only run MC fits",
     cxxopts::value<bool>()->default_value("false"))
    ("vmu", "Flag misid validation PID cuts",
     cxxopts::value<bool>()->default_value("false"))
    ("fake_mu", "Flag fake muon sample PID cuts",
     cxxopts::value<bool>()->default_value("false"))
    ("minos", "Use MINOS in fits to calibration samples",
     cxxopts::value<bool>()->default_value("false"))
    ;
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  // Read arguments

  const TString opath = parsedArgs["output"].as<string>();

  const auto run_simultaneous_mc_calib = parsedArgs["simultaneous"].as<bool>();

  const auto ymlFile   = parsedArgs["config"].as<string>();
  const auto years     = parsedArgs["years"].as<vector<string>>();
  const auto particles = parsedArgs["particles"].as<vector<string>>();
  const auto dry_run   = parsedArgs["dry-run"].as<bool>();
  const auto mc_only   = parsedArgs["mc-only"].as<bool>();
  const auto vmu       = parsedArgs["vmu"].as<bool>();
  const auto fake_mu   = parsedArgs["fake_mu"].as<bool>();
  const auto use_minos = parsedArgs["minos"].as<bool>();

  if (vmu && fake_mu) {
    cout << "WARNING Both VMU and FAKE_MU flags set. The VMU flag is redundant "
            "as the FAKE MU flag implies no muBDT cut."
         << endl;
  } else if (vmu) {
    cout << "INFO Using VMU pid cuts" << endl;
  } else if (fake_mu) {
    cout << "INFO Using FAKE_MU pid cuts" << endl;
  }

  if (dry_run) {
    cout << "WARNING Dry run. Fits will not be performed." << endl;
  } else if (mc_only) {
    cout << "WARNING Only MC fits will be performed." << endl;
  } else if (run_simultaneous_mc_calib) {
    cout << "WARNING Running full simultaneous fits, which can take days!"
         << endl;
  }

  if (use_minos) cout << "INFO Using MINOS in calib fits" << endl;

  constexpr int max_fix_reattempts = 10;

  // Parse YAML config
  const auto ymlConfig = YAML::LoadFile(ymlFile)["misid_corrections"];

  // Define histogram to easily determine kinematical bins
  TH3D histo_binnning("histo_binnning", ";p;#eta;nTracks", N_BINS_P, BINS_P,
                      N_BINS_ETA, BINS_ETA, N_BINS_NTRACKS, BINS_NTRACKS);

  // Define counters for mass window efficiencies and initialize all to 0
  int count_in_mw_passed_dif[3][N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] = {
      {{{0}}}};
  int count_in_mw_passed_nondif[3][N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] = {
      {{{0}}}};
  int count_total_passed_dif[3][N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] = {
      {{{0}}}};
  int count_total_passed_nondif[3][N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] = {
      {{{0}}}};
  int count_in_mw_failed_dif[3][N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] = {
      {{{0}}}};
  int count_in_mw_failed_nondif[3][N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] = {
      {{{0}}}};
  int count_total_failed_dif[3][N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] = {
      {{{0}}}};
  int count_total_failed_nondif[3][N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] = {
      {{{0}}}};

  // Define variable to access input ntuples
  int    dst_id, d0_id, probe_trueid, probe_daughter0_trueid, ntracks;
  double dst_m, d0_m, d0_pt, probe_p, probe_pz, probe_pt, probe_eta,
      probe_dllmu, probe_dlle, probe_ghostprob, probe_mu_unbiased,
      ntracks_calib, spi_p, spi_pt, tag_p, tag_pt, k_track_chi2ndof,
      pi_track_chi2ndof, spi_track_chi2ndof, tag_probnn;
  float probe_mu_ubdt;
  bool  probe_ismuon, probe_hasmuon;

  // Other variables
  TCanvas c_single("c_single", "c_single", 1280, 960);
  TCanvas c_double("c_double", "c_double", 640, 960);
  TCanvas c_four("c_four", "c_four", 1280, 960);
  TCanvas c_mult("c_mult", "c_mult", 2560, 960);
  c_mult.Divide(4, 2);
  c_four.Divide(2, 2);
  c_double.Divide(1, 2);

  // Define output file
  const TString opath_full = opath + "/mc_misid_histos.root";
  cout << "INFO Creating output file: " << opath_full << endl;
  TFile ofile(opath_full, "RECREATE");

  // Observables
  RooRealVar d0_m_var("d0_m_var", "d0_m_var", D0_M_min, D0_M_max, "GeV");
  RooRealVar dm_var("dm_var", "dm_var", DM_min, DM_max, "GeV");
  RooBinning bins_d0_m(85, D0_M_min, D0_M_max, "bins_d0_m");
  RooBinning bins_dm(80, DM_min, DM_max, "bins_dm");
  RooBinning bins_d0_m_passed(50, D0_M_min, D0_M_max, "bins_mc_d0_passed");
  RooBinning bins_dm_passed(50, DM_min, DM_max, "bins_dm_passed");

  // Define categories
  RooCategory sample("sample", "sample");
  sample.defineType("mc_passed_dif");
  sample.defineType("mc_passed_nondif");
  sample.defineType("mc_failed_dif");
  sample.defineType("mc_failed_nondif");
  sample.defineType("calib_passed");
  sample.defineType("calib_failed");

  ///////////////////////////////
  // Individual fit components //
  ///////////////////////////////

  // D0 signal: CrystalBall + Gaussian

  // Core
  RooRealVar mean_d0_failed_dif("mean_d0_failed_dif", "mean_d0_failed_dif",
                                1.855, 1.870, "GeV");
  RooRealVar mean_d0_passed_dif("mean_d0_passed_dif", "mean_d0_passed_dif",
                                1.855, 1.870, "GeV");

  RooRealVar mean_d0_failed_nondif(
      "mean_d0_failed_nondif", "mean_d0_failed_nondif", 1.855, 1.870, "GeV");
  RooRealVar mean_d0_passed_nondif(
      "mean_d0_passed_nondif", "mean_d0_passed_nondif", 1.855, 1.870, "GeV");

  RooRealVar width_d0_dif("width_d0_dif", "width_d0_dif", 0.004, 0.020, "GeV");
  RooRealVar width_d0_nondif("width_d0_nondif", "width_d0_nondif", 0.004, 0.015,
                             "GeV");

  RooRealVar width_cb_d0_dif("width_cb_d0_dif", "width_cb_d0_dif", 0.005, 0.015,
                             "GeV");
  RooRealVar width_cb_d0_nondif("width_cb_d0_nondif", "width_cb_d0_nondif",
                                0.005, 0.015, "GeV");

  RooRealVar m_shift_d0_failed("m_shift_d0_failed", "m_shift_d0_failed", -0.004,
                               0.004, "GeV");
  RooRealVar m_shift_d0_passed("m_shift_d0_passed", "m_shift_d0_passed", -0.004,
                               0.004, "GeV");
  RooRealVar w_scale_d0("w_scale_d0", "w_scale_d0", 0.8, 1.2, "");

  RooFormulaVar mean_shifted_d0_failed_dif(
      "mean_shifted_d0_failed_dif", "mean_shifted_d0_failed_dif", "x[0] + x[1]",
      RooArgList(mean_d0_failed_dif, m_shift_d0_failed));
  RooFormulaVar mean_shifted_d0_failed_nondif(
      "mean_shifted_d0_failed_nondif", "mean_shifted_d0_failed_nondif",
      "x[0] + x[1]", RooArgList(mean_d0_failed_nondif, m_shift_d0_failed));

  RooFormulaVar mean_shifted_d0_passed_dif(
      "mean_shifted_d0_passed_dif", "mean_shifted_d0_passed_dif", "x[0] + x[1]",
      RooArgList(mean_d0_passed_dif, m_shift_d0_passed));
  RooFormulaVar mean_shifted_d0_passed_nondif(
      "mean_shifted_d0_passed_nondif", "mean_shifted_d0_passed_nondif",
      "x[0] + x[1]", RooArgList(mean_d0_passed_nondif, m_shift_d0_passed));

  RooFormulaVar width_d0_scaled_failed_dif(
      "width_d0_scaled_failed_dif", "width_d0_scaled_failed_dif", "x[0] * x[1]",
      RooArgList(width_d0_dif, w_scale_d0));
  RooFormulaVar width_d0_scaled_failed_nondif(
      "width_d0_scaled_failed_nondif", "width_d0_scaled_failed_nondif",
      "x[0] * x[1]", RooArgList(width_d0_nondif, w_scale_d0));
  RooFormulaVar width_d0_scaled_passed_dif(
      "width_d0_scaled_passed_dif", "width_d0_scaled_passed_dif", "x[0] * x[1]",
      RooArgList(width_d0_dif, w_scale_d0));
  RooFormulaVar width_d0_scaled_passed_nondif(
      "width_d0_scaled_passed_nondif", "width_d0_scaled_passed_nondif",
      "x[0] * x[1]", RooArgList(width_d0_nondif, w_scale_d0));

  RooFormulaVar width_cb_d0_scaled_failed_dif(
      "width_cb_d0_scaled_failed_dif", "width_cb_d0_scaled_failed_dif",
      "x[0] * x[1]", RooArgList(width_cb_d0_dif, w_scale_d0));
  RooFormulaVar width_cb_d0_scaled_failed_nondif(
      "width_cb_d0_scaled_failed_nondif", "width_cb_d0_scaled_failed_nondif",
      "x[0] * x[1]", RooArgList(width_cb_d0_nondif, w_scale_d0));
  RooFormulaVar width_cb_d0_scaled_passed_dif(
      "width_cb_d0_scaled_passed_dif", "width_cb_d0_scaled_passed_dif",
      "x[0] * x[1]", RooArgList(width_cb_d0_dif, w_scale_d0));
  RooFormulaVar width_cb_d0_scaled_passed_nondif(
      "width_cb_d0_scaled_passed_nondif", "width_cb_d0_scaled_passed_nondif",
      "x[0] * x[1]", RooArgList(width_cb_d0_nondif, w_scale_d0));

  // Left tail
  RooRealVar alpha_L_d0_failed_dif("alpha_L_d0_failed_dif",
                                   "alpha_L_d0_failed_dif", 1e-1, 3., "");
  RooRealVar alpha_L_d0_passed_dif("alpha_L_d0_passed_dif",
                                   "alpha_L_d0_passed_dif", 1e-3, 3., "");

  RooRealVar alpha_L_d0_failed_nondif("alpha_L_d0_failed_nondif",
                                      "alpha_L_d0_failed_nondif", 1e-1, 3., "");
  RooRealVar alpha_L_d0_passed_nondif("alpha_L_d0_passed_nondif",
                                      "alpha_L_d0_passed_nondif", 1e-1, 3., "");

  RooRealVar n_L_d0_dif("n_L_d0_dif", "n_L_d0_dif", 0.1, 20., "");
  RooRealVar n_L_d0_nondif("n_L_d0_nondif", "n_L_d0_nondif", 0.1, 20., "");

  // Right tail
  RooRealVar alpha_R_failed_d0_dif("alpha_R_failed_d0_dif",
                                   "alpha_R_failed_d0_dif", 2e-1, 3., "");
  RooRealVar alpha_R_passed_d0_dif("alpha_R_passed_d0_dif",
                                   "alpha_R_passed_d0_dif", 2e-1, 3., "");

  RooRealVar alpha_R_failed_d0_nondif("alpha_R_failed_d0_nondif",
                                      "alpha_R_failed_d0_nondif", 2e-1, 3., "");
  RooRealVar alpha_R_passed_d0_nondif("alpha_R_passed_d0_nondif",
                                      "alpha_R_passed_d0_nondif", 2e-1, 3., "");

  RooRealVar n_R_d0_dif("n_R_d0_dif", "n_R_d0_dif", 0.1, 20., "");
  RooRealVar n_R_d0_nondif("n_R_d0_nondif", "n_R_d0_nondif", 0.1, 20., "");

  // Gaussian fractions
  RooRealVar f_d0_passed_dif("f_d0_passed_dif", "f_d0_passed_dif", 0., 1.);
  RooRealVar f_d0_failed_dif("f_d0_failed_dif", "f_d0_failed_dif", 0., 1.);

  RooRealVar f_d0_passed_nondif("f_d0_passed_nondif", "f_d0_passed_nondif", 0.,
                                1.);
  RooRealVar f_d0_failed_nondif("f_d0_failed_nondif", "f_d0_failed_nondif", 0.,
                                1.);

  // MC PDFs
  RooGaussian d0_gauss_failed_mc_dif("d0_gauss_failed_mc_dif",
                                     "d0_gauss_failed_mc_dif", d0_m_var,
                                     mean_d0_failed_dif, width_d0_dif);
  RooGaussian d0_gauss_passed_mc_dif("d0_gauss_passed_mc_dif",
                                     "d0_gauss_passed_mc_dif", d0_m_var,
                                     mean_d0_passed_dif, width_d0_dif);

  RooGaussian d0_gauss_failed_mc_nondif("d0_gauss_failed_mc_nondif",
                                        "d0_gauss_failed_mc_nondif", d0_m_var,
                                        mean_d0_failed_nondif, width_d0_nondif);
  RooGaussian d0_gauss_passed_mc_nondif("d0_gauss_passed_mc_nondif",
                                        "d0_gauss_passed_mc_nondif", d0_m_var,
                                        mean_d0_passed_nondif, width_d0_nondif);

  RooCrystalBall d0_cb_failed_mc_dif(
      "d0_cb_failed_mc_dif", "d0_cb_failed_mc_dif", d0_m_var,
      mean_d0_failed_dif, width_cb_d0_dif, alpha_L_d0_failed_dif, n_L_d0_dif,
      alpha_R_failed_d0_dif, n_R_d0_dif);
  RooCrystalBall d0_cb_passed_mc_dif(
      "d0_cb_passed_mc_dif", "d0_cb_passed_mc_dif", d0_m_var,
      mean_d0_passed_dif, width_cb_d0_dif, alpha_L_d0_passed_dif, n_L_d0_dif,
      alpha_R_passed_d0_dif, n_R_d0_dif);

  RooCrystalBall d0_cb_failed_mc_nondif(
      "d0_cb_failed_mc_nondif", "d0_cb_failed_mc_nondif", d0_m_var,
      mean_d0_failed_nondif, width_cb_d0_nondif, alpha_L_d0_failed_nondif,
      n_L_d0_nondif, alpha_R_failed_d0_nondif, n_R_d0_nondif);
  RooCrystalBall d0_cb_passed_mc_nondif(
      "d0_cb_passed_mc_nondif", "d0_cb_passed_mc_nondif", d0_m_var,
      mean_d0_passed_nondif, width_cb_d0_nondif, alpha_L_d0_passed_nondif,
      n_L_d0_nondif, alpha_R_passed_d0_nondif, n_R_d0_nondif);

  RooAddPdf d0_model_passed_mc_dif(
      "d0_model_passed_mc_dif", "d0_model_passed_mc_dif",
      d0_gauss_passed_mc_dif, d0_cb_passed_mc_dif, f_d0_passed_dif);
  RooAddPdf d0_model_failed_mc_dif(
      "d0_model_failed_mc_dif", "d0_model_failed_mc_dif",
      d0_gauss_failed_mc_dif, d0_cb_failed_mc_dif, f_d0_failed_dif);

  RooAddPdf d0_model_passed_mc_nondif(
      "d0_model_passed_mc_nondif", "d0_model_passed_mc_nondif",
      d0_gauss_passed_mc_nondif, d0_cb_passed_mc_nondif, f_d0_passed_nondif);
  RooAddPdf d0_model_failed_mc_nondif(
      "d0_model_failed_mc_nondif", "d0_model_failed_mc_nondif",
      d0_gauss_failed_mc_nondif, d0_cb_failed_mc_nondif, f_d0_failed_nondif);

  // Data PDFs with floating shift
  RooGaussian d0_gauss_passed_calib_dif(
      "d0_gauss_passed_calib_dif", "d0_gauss_passed_calib_dif", d0_m_var,
      mean_shifted_d0_passed_dif, width_d0_scaled_passed_dif);
  RooGaussian d0_gauss_passed_calib_nondif(
      "d0_gauss_passed_calib_nondif", "d0_gauss_passed_calib_nondif", d0_m_var,
      mean_shifted_d0_passed_nondif, width_d0_scaled_passed_nondif);

  RooGaussian d0_gauss_failed_calib_dif(
      "d0_gauss_failed_calib_dif", "d0_gauss_failed_calib_dif", d0_m_var,
      mean_shifted_d0_failed_dif, width_d0_scaled_failed_dif);
  RooGaussian d0_gauss_failed_calib_nondif(
      "d0_gauss_failed_calib_nondif", "d0_gauss_failed_calib_nondif", d0_m_var,
      mean_shifted_d0_failed_nondif, width_d0_scaled_failed_nondif);

  RooCrystalBall d0_cb_passed_calib_dif(
      "d0_cb_passed_calib_dif", "d0_cb_passed_calib_dif", d0_m_var,
      mean_shifted_d0_passed_dif, width_cb_d0_scaled_passed_dif,
      alpha_L_d0_passed_dif, n_L_d0_dif, alpha_R_passed_d0_dif, n_R_d0_dif);
  RooCrystalBall d0_cb_passed_calib_nondif(
      "d0_cb_passed_calib_nondif", "d0_cb_passed_calib_nondif", d0_m_var,
      mean_shifted_d0_passed_nondif, width_cb_d0_scaled_passed_nondif,
      alpha_L_d0_passed_nondif, n_L_d0_nondif, alpha_R_passed_d0_nondif,
      n_R_d0_nondif);

  RooCrystalBall d0_cb_failed_calib_dif(
      "d0_cb_failed_calib_dif", "d0_cb_failed_calib_dif", d0_m_var,
      mean_shifted_d0_failed_dif, width_cb_d0_scaled_failed_dif,
      alpha_L_d0_failed_dif, n_L_d0_dif, alpha_R_failed_d0_dif, n_R_d0_dif);
  RooCrystalBall d0_cb_failed_calib_nondif(
      "d0_cb_failed_calib_nondif", "d0_cb_failed_calib_nondif", d0_m_var,
      mean_shifted_d0_failed_nondif, width_cb_d0_scaled_failed_nondif,
      alpha_L_d0_failed_nondif, n_L_d0_nondif, alpha_R_failed_d0_nondif,
      n_R_d0_nondif);

  RooAddPdf d0_model_passed_calib_dif(
      "d0_model_passed_calib_dif", "d0_model_passed_calib_dif",
      d0_gauss_passed_calib_dif, d0_cb_passed_calib_dif, f_d0_passed_dif);
  RooAddPdf d0_model_passed_calib_nondif(
      "d0_model_passed_calib_nondif", "d0_model_passed_calib_nondif",
      d0_gauss_passed_calib_nondif, d0_cb_passed_calib_nondif,
      f_d0_passed_nondif);

  RooAddPdf d0_model_failed_calib_dif(
      "d0_model_failed_calib_dif", "d0_model_failed_calib_dif",
      d0_gauss_failed_calib_dif, d0_cb_failed_calib_dif, f_d0_failed_dif);
  RooAddPdf d0_model_failed_calib_nondif(
      "d0_model_failed_calib_nondif", "d0_model_failed_calib_nondif",
      d0_gauss_failed_calib_nondif, d0_cb_failed_calib_nondif,
      f_d0_failed_nondif);

  // dm signal: CrystalBall + Gaussian

  // Core
  RooRealVar mean_dm_failed_dif("mean_dm_failed_dif", "mean_dm_failed_dif",
                                0.145, 0.146, "GeV");
  RooRealVar mean_dm_passed_dif("mean_dm_passed_dif", "mean_dm_passed_dif",
                                0.145, 0.146, "GeV");

  RooRealVar mean_dm_failed_nondif(
      "mean_dm_failed_nondif", "mean_dm_failed_nondif", 0.145, 0.146, "GeV");
  RooRealVar mean_dm_passed_nondif(
      "mean_dm_passed_nondif", "mean_dm_passed_nondif", 0.145, 0.146, "GeV");

  RooRealVar width_dm_dif("width_dm_dif", "width_dm_dif", 0.00021, 0.002,
                          "GeV");
  RooRealVar width_dm_nondif("width_dm_nondif", "width_dm_nondif", 0.00021,
                             0.002, "GeV");

  RooRealVar width_cb_dm_dif("width_cb_dm_dif", "width_cb_dm_dif", 0.0001,
                             0.003, "GeV");
  RooRealVar width_cb_dm_nondif("width_cb_dm_nondif", "width_cb_dm_nondif",
                                0.0001, 0.003, "GeV");

  RooRealVar m_shift_dm("m_shift_dm", "m_shift_dm", -0.0001, 0.0001, "GeV");
  RooRealVar w_scale_dm("w_scale_dm", "w_scale_dm", 0.8, 1.2, "");

  RooFormulaVar mean_shifted_dm_failed_dif(
      "mean_shifted_dm_failed_dif", "mean_shifted_dm_failed_dif", "x[0] + x[1]",
      RooArgList(mean_dm_failed_dif, m_shift_dm));
  RooFormulaVar mean_shifted_dm_failed_nondif(
      "mean_shifted_dm_failed_nondif", "mean_shifted_dm_failed_nondif",
      "x[0] + x[1]", RooArgList(mean_dm_failed_nondif, m_shift_dm));
  RooFormulaVar mean_shifted_dm_passed_dif(
      "mean_shifted_dm_passed_dif", "mean_shifted_dm_passed_dif", "x[0] + x[1]",
      RooArgList(mean_dm_passed_dif, m_shift_dm));
  RooFormulaVar mean_shifted_dm_passed_nondif(
      "mean_shifted_dm_passed_nondif", "mean_shifted_dm_passed_nondif",
      "x[0] + x[1]", RooArgList(mean_dm_passed_nondif, m_shift_dm));

  RooFormulaVar width_dm_scaled_failed_dif(
      "width_dm_scaled_failed_dif", "width_dm_scaled_failed_dif", "x[0] * x[1]",
      RooArgList(width_dm_dif, w_scale_dm));
  RooFormulaVar width_dm_scaled_failed_nondif(
      "width_dm_scaled_failed_nondif", "width_dm_scaled_failed_nondif",
      "x[0] * x[1]", RooArgList(width_dm_nondif, w_scale_dm));
  RooFormulaVar width_dm_scaled_passed_dif(
      "width_dm_scaled_passed_dif", "width_dm_scaled_passed_dif", "x[0] * x[1]",
      RooArgList(width_dm_dif, w_scale_dm));
  RooFormulaVar width_dm_scaled_passed_nondif(
      "width_dm_scaled_passed_nondif", "width_dm_scaled_passed_nondif",
      "x[0] * x[1]", RooArgList(width_dm_nondif, w_scale_dm));

  RooFormulaVar width_cb_dm_scaled_failed_dif(
      "width_cb_dm_scaled_failed_dif", "width_cb_dm_scaled_failed_dif",
      "x[0] * x[1]", RooArgList(width_cb_dm_dif, w_scale_dm));
  RooFormulaVar width_cb_dm_scaled_failed_nondif(
      "width_cb_dm_scaled_failed_nondif", "width_cb_dm_scaled_failed_nondif",
      "x[0] * x[1]", RooArgList(width_cb_dm_nondif, w_scale_dm));
  RooFormulaVar width_cb_dm_scaled_passed_dif(
      "width_cb_dm_scaled_passed_dif", "width_cb_dm_scaled_passed_dif",
      "x[0] * x[1]", RooArgList(width_cb_dm_dif, w_scale_dm));
  RooFormulaVar width_cb_dm_scaled_passed_nondif(
      "width_cb_dm_scaled_passed_nondif", "width_cb_dm_scaled_passed_nondif",
      "x[0] * x[1]", RooArgList(width_cb_dm_nondif, w_scale_dm));

  // Left tail
  RooRealVar alpha_L_dm_dif("alpha_L_dm_dif", "alpha_L_dm_dif", 1e-1, 3., "");
  RooRealVar alpha_L_dm_nondif("alpha_L_dm_nondif", "alpha_L_dm_nondif", 1e-1,
                               3., "");
  RooRealVar n_L_dm_dif("n_L_dm_dif", "n_L_dm_dif", 0.1, 20., "");
  RooRealVar n_L_dm_nondif("n_L_dm_nondif", "n_L_dm_nondif", 0.1, 20., "");

  // Right tail
  RooRealVar alpha_R_dm_dif("alpha_R_dm_dif", "alpha_R_dm_dif", 1e-1, 3., "");
  RooRealVar alpha_R_dm_nondif("alpha_R_dm_nondif", "alpha_R_dm_nondif", 1e-1,
                               3., "");
  RooRealVar n_R_dm_dif("n_R_dm_dif", "n_R_dm_dif", 0.1, 20., "");
  RooRealVar n_R_dm_nondif("n_R_dm_nondif", "n_R_dm_nondif", 0.1, 20., "");

  // Gaussian fraction
  RooRealVar f_dm_passed_dif("f_dm_passed_dif", "f_dm_passed_dif", 0., 1.);
  RooRealVar f_dm_failed_dif("f_dm_failed_dif", "f_dm_failed_dif", 0., 1.);

  RooRealVar f_dm_passed_nondif("f_dm_passed_nondif", "f_dm_passed_nondif", 0.,
                                1.);
  RooRealVar f_dm_failed_nondif("f_dm_failed_nondif", "f_dm_failed_nondif", 0.,
                                1.);

  // MC PDFs
  RooGaussian dm_gauss_mc_failed_dif("dm_gauss_mc_failed_dif",
                                     "dm_gauss_mc_failed_dif", dm_var,
                                     mean_dm_failed_dif, width_dm_dif);
  RooGaussian dm_gauss_mc_passed_dif("dm_gauss_mc_passed_dif",
                                     "dm_gauss_mc_passed_dif", dm_var,
                                     mean_dm_passed_dif, width_dm_dif);

  RooGaussian dm_gauss_mc_failed_nondif("dm_gauss_mc_failed_nondif",
                                        "dm_gauss_mc_failed_nondif", dm_var,
                                        mean_dm_failed_nondif, width_dm_nondif);
  RooGaussian dm_gauss_mc_passed_nondif("dm_gauss_mc_passed_nondif",
                                        "dm_gauss_mc_passed_nondif", dm_var,
                                        mean_dm_passed_nondif, width_dm_nondif);

  RooCrystalBall dm_cb_mc_failed_dif(
      "dm_cb_mc_failed_dif", "dm_cb_mc_failed_dif", dm_var, mean_dm_failed_dif,
      width_cb_dm_dif, alpha_L_dm_dif, n_L_dm_dif, alpha_R_dm_dif, n_R_dm_dif);
  RooCrystalBall dm_cb_mc_passed_dif(
      "dm_cb_mc_passed_dif", "dm_cb_mc_passed_dif", dm_var, mean_dm_passed_dif,
      width_cb_dm_dif, alpha_L_dm_dif, n_L_dm_dif, alpha_R_dm_dif, n_R_dm_dif);

  RooCrystalBall dm_cb_mc_failed_nondif(
      "dm_cb_mc_failed_nondif", "dm_cb_mc_failed_nondif", dm_var,
      mean_dm_failed_nondif, width_cb_dm_nondif, alpha_L_dm_nondif,
      n_L_dm_nondif, alpha_R_dm_nondif, n_R_dm_nondif);
  RooCrystalBall dm_cb_mc_passed_nondif(
      "dm_cb_mc_passed_nondif", "dm_cb_mc_passed_nondif", dm_var,
      mean_dm_passed_nondif, width_cb_dm_nondif, alpha_L_dm_nondif,
      n_L_dm_nondif, alpha_R_dm_nondif, n_R_dm_nondif);

  RooAddPdf dm_model_passed_mc_dif(
      "dm_model_passed_mc_dif", "dm_model_passed_mc_dif",
      dm_gauss_mc_passed_dif, dm_cb_mc_passed_dif, f_dm_passed_dif);
  RooAddPdf dm_model_failed_mc_dif(
      "dm_model_failed_mc_dif", "dm_model_failed_mc_dif",
      dm_gauss_mc_failed_dif, dm_cb_mc_failed_dif, f_dm_failed_dif);

  RooAddPdf dm_model_passed_mc_nondif(
      "dm_model_passed_mc_nondif", "dm_model_passed_mc_nondif",
      dm_gauss_mc_passed_nondif, dm_cb_mc_passed_nondif, f_dm_passed_nondif);
  RooAddPdf dm_model_failed_mc_nondif(
      "dm_model_failed_mc_nondif", "dm_model_failed_mc_nondif",
      dm_gauss_mc_failed_nondif, dm_cb_mc_failed_nondif, f_dm_failed_nondif);

  // Data PDFs with floating shift
  RooGaussian dm_gauss_passed_calib_dif(
      "dm_gauss_passed_calib_dif", "dm_gauss_passed_calib_dif", dm_var,
      mean_shifted_dm_passed_dif, width_dm_scaled_passed_dif);
  RooGaussian dm_gauss_passed_calib_nondif(
      "dm_gauss_passed_calib_nondif", "dm_gauss_passed_calib_nondif", dm_var,
      mean_shifted_dm_passed_nondif, width_dm_scaled_passed_nondif);

  RooGaussian dm_gauss_failed_calib_dif(
      "dm_gauss_failed_calib_dif", "dm_gauss_failed_calib_dif", dm_var,
      mean_shifted_dm_failed_dif, width_dm_scaled_failed_dif);
  RooGaussian dm_gauss_failed_calib_nondif(
      "dm_gauss_failed_calib_nondif", "dm_gauss_failed_calib_nondif", dm_var,
      mean_shifted_dm_failed_nondif, width_dm_scaled_failed_nondif);

  RooCrystalBall dm_cb_passed_calib_dif(
      "dm_cb_passed_calib_dif", "dm_cb_passed_calib_dif", dm_var,
      mean_shifted_dm_passed_dif, width_cb_dm_scaled_passed_dif, alpha_L_dm_dif,
      n_L_dm_dif, alpha_R_dm_dif, n_R_dm_dif);
  RooCrystalBall dm_cb_passed_calib_nondif(
      "dm_cb_passed_calib_nondif", "dm_cb_passed_calib_nondif", dm_var,
      mean_shifted_dm_passed_nondif, width_cb_dm_scaled_passed_nondif,
      alpha_L_dm_nondif, n_L_dm_nondif, alpha_R_dm_nondif, n_R_dm_nondif);

  RooCrystalBall dm_cb_failed_calib_dif(
      "dm_cb_failed_calib_dif", "dm_cb_failed_calib_dif", dm_var,
      mean_shifted_dm_failed_dif, width_cb_dm_scaled_failed_dif, alpha_L_dm_dif,
      n_L_dm_dif, alpha_R_dm_dif, n_R_dm_dif);
  RooCrystalBall dm_cb_failed_calib_nondif(
      "dm_cb_failed_calib_nondif", "dm_cb_failed_calib_nondif", dm_var,
      mean_shifted_dm_failed_nondif, width_cb_dm_scaled_failed_nondif,
      alpha_L_dm_nondif, n_L_dm_nondif, alpha_R_dm_nondif, n_R_dm_nondif);

  RooAddPdf dm_model_passed_calib_dif(
      "dm_model_passed_calib_dif", "dm_model_passed_calib_dif",
      dm_gauss_passed_calib_dif, dm_cb_passed_calib_dif, f_dm_passed_dif);
  RooAddPdf dm_model_passed_calib_nondif(
      "dm_model_passed_calib_nondif", "dm_model_passed_calib_nondif",
      dm_gauss_passed_calib_nondif, dm_cb_passed_calib_nondif,
      f_dm_passed_nondif);

  RooAddPdf dm_model_failed_calib_dif(
      "dm_model_failed_calib_dif", "dm_model_failed_calib_dif",
      dm_gauss_failed_calib_dif, dm_cb_failed_calib_dif, f_dm_failed_dif);
  RooAddPdf dm_model_failed_calib_nondif(
      "dm_model_failed_calib_nondif", "dm_model_failed_calib_nondif",
      dm_gauss_failed_calib_nondif, dm_cb_failed_calib_nondif,
      f_dm_failed_nondif);

  // dm comb model: threshold function
  RooConstVar dm0("dm0", "dm0", PI_M);

  RooRealVar c_comb_spi("c_comb_spi", "c_comb_spi", 0., 1., "");
  RooRealVar c_comb_all("c_comb_all", "c_comb_all", 0., 1., "");

  RooPowerLaw dm_comb_spi_failed("dm_comb_spi_failed", "dm_comb_spi_failed",
                                 dm_var, dm0, c_comb_spi);
  RooPowerLaw dm_comb_spi_passed("dm_comb_spi_passed", "dm_comb_spi_passed",
                                 dm_var, dm0, c_comb_spi);

  RooPowerLaw dm_comb_all_failed("dm_comb_all_failed", "dm_comb_all_failed",
                                 dm_var, dm0, c_comb_all);
  RooPowerLaw dm_comb_all_passed("dm_comb_all_passed", "dm_comb_all_passed",
                                 dm_var, dm0, c_comb_all);

  // D0 comb model: exponential distribution
  RooRealVar k_part_reco_failed("k_part_reco_failed", "k_part_reco_failed",
                                -20., 0., "");
  RooRealVar k_part_reco_passed("k_part_reco_passed", "k_part_reco_passed",
                                -20., 0., "");

  RooRealVar k_comb_all("k_comb_all", "k_comb_all", -20., 0., "");

  RooExponential d0_part_reco_failed("d0_part_reco_failed",
                                     "d0_part_reco_failed", d0_m_var,
                                     k_part_reco_failed);
  RooExponential d0_part_reco_passed("d0_part_reco_passed",
                                     "d0_part_reco_passed", d0_m_var,
                                     k_part_reco_passed);

  RooExponential d0_comb_all_failed("d0_comb_all_failed", "d0_comb_all_failed",
                                    d0_m_var, k_comb_all);
  RooExponential d0_comb_all_passed("d0_comb_all_passed", "d0_comb_all_passed",
                                    d0_m_var, k_comb_all);

  ///////////////////////////
  // PDFs for MC-only fits //
  ///////////////////////////

  // D0 mass
  RooSimultaneous d0_model_mc_dif("d0_model_mc_dif", "d0_model_mc_dif",
                                  {{"mc_passed_dif", &d0_model_passed_mc_dif},
                                   {"mc_failed_dif", &d0_model_failed_mc_dif}},
                                  sample);

  RooSimultaneous d0_model_mc_nondif(
      "d0_model_mc_nondif", "d0_model_mc_nondif",
      {{"mc_passed_nondif", &d0_model_passed_mc_nondif},
       {"mc_failed_nondif", &d0_model_failed_mc_nondif}},
      sample);

  // dm
  RooSimultaneous dm_model_mc_dif("dm_model_mc_dif", "dm_model_mc_dif",
                                  {{"mc_passed_dif", &dm_model_passed_mc_dif},
                                   {"mc_failed_dif", &dm_model_failed_mc_dif}},
                                  sample);

  RooSimultaneous dm_model_mc_nondif(
      "dm_model_mc_nondif", "dm_model_mc_nondif",
      {{"mc_passed_nondif", &dm_model_passed_mc_nondif},
       {"mc_failed_nondif", &dm_model_failed_mc_nondif}},
      sample);

  // Full fit
  RooRealVar n_inmw_passed_dif("n_inmw_passed_dif", "n_inmw_passed_dif", 0.,
                               "");
  RooRealVar n_inmw_passed_nondif("n_inmw_passed_nondif",
                                  "n_inmw_passed_nondif", 0., "");
  RooRealVar n_inmw_failed_dif("n_inmw_failed_dif", "n_inmw_failed_dif", 0.,
                               "");
  RooRealVar n_inmw_failed_nondif("n_inmw_failed_nondif",
                                  "n_inmw_failed_nondif", 0., "");
  RooRealVar n_total_passed_dif("n_total_passed_dif", "n_total_passed_dif", 0.,
                                "");
  RooRealVar n_total_passed_nondif("n_total_passed_nondif",
                                   "n_total_passed_nondif", 0., "");
  RooRealVar n_total_failed_dif("n_total_failed_dif", "n_total_failed_dif", 0.,
                                "");
  RooRealVar n_total_failed_nondif("n_total_failed_nondif",
                                   "n_total_failed_nondif", 0., "");

  RooFormulaVar f_dif_mc_passed(
      "f_dif_mc_passed", "f_dif_mc_passed", "x[0]/(x[0]+x[1])",
      RooArgList(n_inmw_passed_dif, n_inmw_passed_nondif));
  RooFormulaVar f_dif_mc_failed(
      "f_dif_mc_failed", "f_dif_mc_failed", "x[0]/(x[0]+x[1])",
      RooArgList(n_inmw_failed_dif, n_inmw_failed_nondif));

  // Scale controlling increase/decrease of mis-identified non-dif tracks.
  // 0 means no non-dif misid; 1 means all tracks are misid'ed.
  RooRealVar scale_nondif_calib("scale_nondif_calib", "scale_nondif_calib", 0.,
                                1., "");
  // Number of non-dif events shifted from failed to passed
  RooFormulaVar nondif_yield_add_passed(
      "nondif_yield_add_passed", "nondif_yield_add_passed",
      "(x[2]+x[1])*x[0]-x[1]",
      RooArgList(scale_nondif_calib, n_inmw_passed_nondif,
                 n_inmw_failed_nondif));

  RooFormulaVar f_dif_calib_passed(
      "f_dif_calib_passed", "f_dif_calib_passed", "x[0]/(x[0]+x[1]+x[2])",
      RooArgList(n_inmw_passed_dif, n_inmw_passed_nondif,
                 nondif_yield_add_passed));
  RooFormulaVar f_dif_calib_failed(
      "f_dif_calib_failed", "f_dif_calib_failed", "x[0]/(x[0]+x[1]-x[2])",
      RooArgList(n_inmw_failed_dif, n_inmw_failed_nondif,
                 nondif_yield_add_passed));

  RooAddPdf d0_model_passed_mc("d0_model_passed_mc", "d0_model_passed_mc",
                               d0_model_passed_mc_dif,
                               d0_model_passed_mc_nondif, f_dif_mc_passed);
  RooAddPdf dm_model_passed_mc("dm_model_passed_mc", "dm_model_passed_mc",
                               dm_model_passed_mc_dif,
                               dm_model_passed_mc_nondif, f_dif_mc_passed);
  RooAddPdf d0_model_failed_mc("d0_model_failed_mc", "d0_model_failed_mc",
                               d0_model_failed_mc_dif,
                               d0_model_failed_mc_nondif, f_dif_mc_failed);
  RooAddPdf dm_model_failed_mc("dm_model_failed_mc", "dm_model_failed_mc",
                               dm_model_failed_mc_dif,
                               dm_model_failed_mc_nondif, f_dif_mc_failed);

  RooAddPdf d0_model_passed_calib(
      "d0_model_passed_calib", "d0_model_passed_calib",
      d0_model_passed_calib_dif, d0_model_passed_calib_nondif,
      f_dif_calib_passed);
  RooAddPdf dm_model_passed_calib(
      "dm_model_passed_calib", "dm_model_passed_calib",
      dm_model_passed_calib_dif, dm_model_passed_calib_nondif,
      f_dif_calib_passed);
  RooAddPdf d0_model_failed_calib(
      "d0_model_failed_calib", "d0_model_failed_calib",
      d0_model_failed_calib_dif, d0_model_failed_calib_nondif,
      f_dif_calib_failed);
  RooAddPdf dm_model_failed_calib(
      "dm_model_failed_calib", "dm_model_failed_calib",
      dm_model_failed_calib_dif, dm_model_failed_calib_nondif,
      f_dif_calib_failed);

  RooProdPdf sig_passed_mc("sig_passed_mc", "sig_passed_mc", d0_model_passed_mc,
                           dm_model_passed_mc);
  RooProdPdf sig_failed_mc("sig_failed_mc", "sig_failed_mc", d0_model_failed_mc,
                           dm_model_failed_mc);
  RooProdPdf sig_passed_calib("sig_passed_calib", "sig_passed_calib",
                              d0_model_passed_calib, dm_model_passed_calib);
  RooProdPdf sig_failed_calib("sig_failed_calib", "sig_failed_calib",
                              d0_model_failed_calib, dm_model_failed_calib);

  RooProdPdf dm_comb_passed("dm_comb_passed", "dm_comb_passed",
                            d0_model_passed_calib, dm_comb_spi_passed);
  RooProdPdf dm_comb_failed("dm_comb_failed", "dm_comb_failed",
                            d0_model_failed_calib, dm_comb_spi_failed);

  RooProdPdf d0_comb_passed("d0_comb_passed", "d0_comb_passed",
                            d0_part_reco_passed, dm_model_passed_calib);
  RooProdPdf d0_comb_failed("d0_comb_failed", "d0_comb_failed",
                            d0_part_reco_failed, dm_model_failed_calib);

  RooProdPdf dmd0_comb_passed("dmd0_comb_passed", "dmd0_comb_passed",
                              d0_comb_all_passed, dm_comb_all_passed);
  RooProdPdf dmd0_comb_failed("dmd0_comb_failed", "dmd0_comb_failed",
                              d0_comb_all_failed, dm_comb_all_failed);

  RooRealVar f1_bkg_passed("f1_bkg_passed", "f1_bkg_passed", 0., 1., "");
  RooRealVar f2_bkg_passed("f2_bkg_passed", "f2_bkg_passed", 0., 1., "");
  RooRealVar f1_bkg_failed("f1_bkg_failed", "f1_bkg_failed", 0., 1., "");
  RooRealVar f2_bkg_failed("f2_bkg_failed", "f2_bkg_failed", 0., 1., "");

  RooAddPdf bkg_passed(
      "bkg_passed", "bkg_passed",
      RooArgList(dm_comb_passed, d0_comb_passed, dmd0_comb_passed),
      RooArgList(f1_bkg_passed, f2_bkg_passed), true);
  RooAddPdf bkg_failed(
      "bkg_failed", "bkg_failed",
      RooArgList(dm_comb_failed, d0_comb_failed, dmd0_comb_failed),
      RooArgList(f1_bkg_failed, f2_bkg_failed), true);

  RooRealVar n_mc_passed("n_mc_passed", "n_mc_passed", 0., "");
  RooRealVar n_mc_failed("n_mc_failed", "n_mc_failed", 0., "");

  RooExtendPdf model_mc_passed("model_mc_passed", "model_mc_passed",
                               sig_passed_mc, n_mc_passed);
  RooExtendPdf model_mc_failed("model_mc_failed", "model_mc_failed",
                               sig_failed_mc, n_mc_failed);

  // PID efficiency
  RooRealVar eff("eff", "eff", 0., 1., "");
  if (fake_mu)
    eff.setMin(0.70);
  else
    eff.setMax(0.05);

  // PIDCalib mass window efficiencies
  RooFormulaVar eff_mw_passed(
      "eff_mw_passed", "eff_mw_passed",
      "(x[0] + x[1] + x[4]) / (x[2] + x[3] + x[4] * x[3]/x[1])",
      RooArgList(n_inmw_passed_dif, n_inmw_passed_nondif, n_total_passed_dif,
                 n_total_passed_nondif, nondif_yield_add_passed));
  RooFormulaVar eff_mw_failed(
      "eff_mw_failed", "eff_mw_failed",
      "(x[0] + x[1] - x[4]) / (x[2] + x[3] - x[4] * x[3]/x[1])",
      RooArgList(n_inmw_failed_dif, n_inmw_failed_nondif, n_total_failed_dif,
                 n_total_failed_nondif, nondif_yield_add_passed));

  RooFormulaVar eff_corrected("eff_corrected", "eff_corrected",
                              "(x[0]/x[1]) / ((x[0]/x[1]) + ((1-x[0])/x[2]))",
                              RooArgList(eff, eff_mw_passed, eff_mw_failed));

  // Normalizations
  RooRealVar n_sig("n_sig", "n_sig", 0, 1, "");
  RooRealVar n_bkg_passed("n_bkg_passed", "n_bkg_passed", 0., 1., "");
  RooRealVar n_bkg_failed("n_bkg_failed", "n_bkg_failed", 0., 1., "");

  //   RooFormulaVar n_sig_passed("n_sig_passed", "n_sig_passed",
  //   "x[0]*x[1]*x[2]",
  //                              RooArgList(n_sig, eff, eff_mw_passed));
  //   RooFormulaVar n_sig_failed("n_sig_failed", "n_sig_failed",
  //                              "x[0] * (1. - x[1]) * x[2]",
  //                              RooArgList(n_sig, eff, eff_mw_failed));

  // The definitions above would allow the mass window efficiencies to be
  // incorporated in the fit and automatically adjusted by the floating amount
  // of non-dif events. HOWEVER, the fit becomes much slower. So we fit the
  // uncorrected efficiency, and correct it afterwards.
  RooFormulaVar n_sig_passed("n_sig_passed", "n_sig_passed", "x[0]*x[1]",
                             RooArgList(n_sig, eff));
  RooFormulaVar n_sig_failed("n_sig_failed", "n_sig_failed",
                             "x[0] * (1. - x[1])", RooArgList(n_sig, eff));

  RooAddPdf model_passed("model_passed", "model_passed",
                         RooArgList(sig_passed_calib, bkg_passed),
                         RooArgList(n_sig_passed, n_bkg_passed));
  RooAddPdf model_failed("model_failed", "model_failed",
                         RooArgList(sig_failed_calib, bkg_failed),
                         RooArgList(n_sig_failed, n_bkg_failed));

  // PDF for simultaneous calib fit
  RooSimultaneous model_calib(
      "model_calib", "model_calib",
      {{"calib_passed", &model_passed}, {"calib_failed", &model_failed}},
      sample);

  // PDF for simultaneous MC + calib fit
  RooSimultaneous model("model", "model",
                        {{"mc_passed_dif", &model_mc_passed},
                         {"mc_passed_nondif", &model_mc_passed},
                         {"mc_failed_dif", &model_mc_failed},
                         {"mc_failed_nondif", &model_mc_failed},
                         {"calib_passed", &model_passed},
                         {"calib_failed", &model_failed}},
                        sample);

  auto set_parameters_d0_dif = [&](const string &probe) {
    if (probe == "k") {
      mean_d0_failed_dif.setVal(D0_M);
      mean_d0_passed_dif.setVal(D0_M);
      width_d0_dif.setVal(0.008);
      width_cb_d0_dif.setVal(0.008);
      alpha_L_d0_failed_dif.setVal(0.7);
      alpha_L_d0_passed_dif.setVal(0.9);
      n_L_d0_dif.setVal(3.5);
      alpha_R_failed_d0_dif.setVal(0.9);
      alpha_R_passed_d0_dif.setVal(0.9);
      n_R_d0_dif.setVal(7.);
      f_d0_passed_dif.setVal(0.2);
      f_d0_failed_dif.setVal(0.2);
    } else if (probe == "pi") {
      mean_d0_failed_dif.setVal(D0_M);
      mean_d0_passed_dif.setVal(D0_M);
      width_d0_dif.setVal(0.007);
      width_cb_d0_dif.setVal(0.010);
      alpha_L_d0_failed_dif.setVal(1.4);
      alpha_L_d0_passed_dif.setVal(1.4);
      n_L_d0_dif.setVal(5.5);
      alpha_R_failed_d0_dif.setVal(1.7);
      alpha_R_passed_d0_dif.setVal(1.7);
      n_R_d0_dif.setVal(7.);
      f_d0_passed_dif.setVal(0.2);
      f_d0_failed_dif.setVal(0.2);
    }
  };

  auto set_parameters_d0_nondif = [&](const string &probe) {
    if (probe == "k") {
      mean_d0_failed_nondif.setVal(D0_M);
      mean_d0_passed_nondif.setVal(D0_M);
      width_d0_nondif.setVal(0.008);
      width_cb_d0_nondif.setVal(0.008);
      alpha_L_d0_failed_nondif.setVal(1.1);
      alpha_L_d0_passed_nondif.setVal(1.2);
      n_L_d0_nondif.setVal(5);
      alpha_R_failed_d0_nondif.setVal(1.1);
      alpha_R_passed_d0_nondif.setVal(1.1);
      n_R_d0_nondif.setVal(6.5);
      f_d0_passed_nondif.setVal(0.2);
      f_d0_failed_nondif.setVal(0.2);
    } else if (probe == "pi") {
      mean_d0_failed_nondif.setVal(D0_M);
      mean_d0_passed_nondif.setVal(D0_M);
      width_d0_nondif.setVal(0.007);
      width_cb_d0_nondif.setVal(0.010);
      alpha_L_d0_failed_nondif.setVal(0.6);
      alpha_L_d0_passed_nondif.setVal(0.6);
      n_L_d0_nondif.setVal(5);
      alpha_R_failed_d0_nondif.setVal(1.8);
      alpha_R_passed_d0_nondif.setVal(1.8);
      n_R_d0_nondif.setVal(6.5);
      f_d0_passed_nondif.setVal(0.2);
      f_d0_failed_nondif.setVal(0.2);
    }
  };

  auto set_parameters_dm_dif = [&](const string &probe) {
    if (probe == "k") {
      mean_dm_failed_dif.setVal(DM);
      mean_dm_passed_dif.setVal(DM);
      width_dm_dif.setVal(0.0006);
      width_cb_dm_dif.setVal(0.0007);
      alpha_L_dm_dif.setVal(1.4);
      n_L_dm_dif.setVal(5.9);
      alpha_R_dm_dif.setVal(1.2);
      n_R_dm_dif.setVal(5.1);
      f_dm_passed_dif.setVal(0.2);
      f_dm_failed_dif.setVal(0.2);
    } else if (probe == "pi") {
      mean_dm_failed_dif.setVal(DM);
      mean_dm_passed_dif.setVal(DM);
      width_dm_dif.setVal(0.0005);
      width_cb_dm_dif.setVal(0.0008);
      alpha_L_dm_dif.setVal(1.6);
      n_L_dm_dif.setVal(6.);
      alpha_R_dm_dif.setVal(1.3);
      n_R_dm_dif.setVal(5.);
      f_dm_passed_dif.setVal(0.2);
      f_dm_failed_dif.setVal(0.2);
    }
  };

  auto set_parameters_dm_nondif = [&](const string &probe) {
    if (probe == "k") {
      mean_dm_failed_nondif.setVal(DM);
      mean_dm_passed_nondif.setVal(DM);
      width_dm_nondif.setVal(0.0006);
      width_cb_dm_nondif.setVal(0.0007);
      alpha_L_dm_nondif.setVal(1.2);
      n_L_dm_nondif.setVal(6.5);
      alpha_R_dm_nondif.setVal(1.0);
      n_R_dm_nondif.setVal(5.5);
      f_dm_failed_nondif.setVal(0.2);
      f_dm_passed_nondif.setVal(0.2);
    } else if (probe == "pi") {
      mean_dm_failed_nondif.setVal(DM);
      mean_dm_passed_nondif.setVal(DM);
      width_dm_nondif.setVal(0.0005);
      width_cb_dm_nondif.setVal(0.0008);
      alpha_L_dm_nondif.setVal(1.5);
      n_L_dm_nondif.setVal(6.);
      alpha_R_dm_nondif.setVal(1.2);
      n_R_dm_nondif.setVal(5.5);
      f_dm_failed_nondif.setVal(0.2);
      f_dm_passed_nondif.setVal(0.2);
    }
  };

  auto set_parameters_calib = [&](const string &probe) {
    m_shift_d0_failed.setVal(0.);
    m_shift_d0_passed.setVal(0.);
    w_scale_d0.setVal(1.);
    m_shift_dm.setVal(0.);
    w_scale_dm.setVal(1.);
    scale_nondif_calib.setVal(
        n_inmw_passed_nondif.getVal() /
        (n_inmw_passed_nondif.getVal() + n_inmw_failed_nondif.getVal()));
    if (probe == "k") {
      c_comb_spi.setVal(0.50);
      c_comb_all.setVal(0.32);
    } else if (probe == "pi") {
      c_comb_spi.setVal(0.55);
      c_comb_all.setVal(0.33);
    }
  };

  RooArgSet fit_vars(d0_m_var, dm_var);

  RooArgSet argset_minos(eff, scale_nondif_calib);

  RooArgSet params_d0_dif(
      mean_d0_failed_dif, mean_d0_passed_dif, width_d0_dif, width_cb_d0_dif,
      width_cb_d0_nondif, alpha_L_d0_failed_dif, alpha_L_d0_passed_dif,
      n_L_d0_dif, alpha_R_failed_d0_dif, alpha_R_passed_d0_dif, n_R_d0_dif,
      f_d0_passed_dif, f_d0_failed_dif);

  RooArgSet params_d0_nondif(
      mean_d0_failed_nondif, mean_d0_passed_nondif, width_d0_nondif,
      width_cb_d0_nondif, alpha_L_d0_failed_nondif, alpha_L_d0_passed_nondif,
      n_L_d0_nondif, alpha_R_failed_d0_nondif, alpha_R_passed_d0_nondif,
      n_R_d0_nondif, f_d0_passed_nondif, f_d0_failed_nondif);

  RooArgSet params_dm_dif(mean_dm_failed_dif, mean_dm_passed_dif, width_dm_dif,
                          width_cb_dm_dif, alpha_L_dm_dif, n_L_dm_dif,
                          alpha_R_dm_dif, n_R_dm_dif, f_dm_passed_dif,
                          f_dm_failed_dif);

  RooArgSet params_dm_nondif(
      mean_dm_failed_nondif, mean_dm_passed_nondif, width_dm_nondif,
      width_cb_dm_nondif, alpha_L_dm_nondif, n_L_dm_nondif, alpha_R_dm_nondif,
      n_R_dm_nondif, f_dm_failed_nondif, f_dm_passed_nondif);

  RooFormulaVar m_shift_d0_delta(
      "m_shift_d0_delta", "m_shift_d0_delta", "x[0]-x[1]",
      RooArgList(m_shift_d0_passed, m_shift_d0_failed));
  RooRealVar m_shift_d0_pull("m_shift_d0_pull", "m_shift_d0_pull", -5., 5., "");

  RooArgSet params_calib(m_shift_d0_failed, m_shift_d0_passed, w_scale_d0,
                         m_shift_dm, w_scale_dm, c_comb_spi, c_comb_all,
                         k_part_reco_failed, k_part_reco_passed, k_comb_all,
                         eff, f1_bkg_passed, f2_bkg_passed, f1_bkg_failed,
                         f2_bkg_failed, m_shift_d0_pull, scale_nondif_calib);

  for (auto probe : particles) {
    cout << "INFO Selecting " << probe << endl;

    TString tag;
    if (probe == "pi")
      tag = "k";
    else
      tag = "pi";

    // Define MC datasets to fit
    // With the dif and non-dif events separated, there is only a small
    // difference in the D0 mass shape within each population for different
    // multiplicity bins. Therefore, for the dif events (which have lower
    // stats), the multipliciply bins are combined. For the non-dif population,
    // stats are much higher and merging is not benefitial. In fact, due to the
    // higher precision, the difference actually becomes noticeable (although
    // small).
    RooDataSet *datasets_mc_passed_dif[N_BINS_ETA][N_BINS_P] = {{nullptr}};
    RooDataSet *datasets_mc_failed_dif[N_BINS_ETA][N_BINS_P] = {{nullptr}};
    RooDataSet *datasets_mc_passed_nondif[N_BINS_NTRACKS][N_BINS_ETA]
                                         [N_BINS_P] = {{{nullptr}}};
    RooDataSet *datasets_mc_failed_nondif[N_BINS_NTRACKS][N_BINS_ETA]
                                         [N_BINS_P] = {{{nullptr}}};

    // Initialize MC datasets
    cout << "INFO Initializing MC datasets " << endl;
    for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
      for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
        TString suffix =
            TString::Format("%s_%d_%d", probe.c_str(), eta_idx, p_idx);
        datasets_mc_passed_dif[eta_idx][p_idx] =
            new RooDataSet("dataset_mc_passed_dif_" + suffix,
                           "dataset_mc_passed_dif_" + suffix, fit_vars);
        datasets_mc_failed_dif[eta_idx][p_idx] =
            new RooDataSet("dataset_mc_failed_dif_" + suffix,
                           "dataset_mc_failed_dif_" + suffix, fit_vars);
        for (int ntrks_idx = 0; ntrks_idx < N_BINS_NTRACKS; ntrks_idx++) {
          suffix.Form("%s_%d_%d_%d", probe.c_str(), ntrks_idx, eta_idx, p_idx);
          datasets_mc_passed_nondif[ntrks_idx][eta_idx][p_idx] =
              new RooDataSet("dataset_mc_passed_nondif_" + suffix,
                             "dataset_mc_passed_nondif_" + suffix, fit_vars);
          datasets_mc_failed_nondif[ntrks_idx][eta_idx][p_idx] =
              new RooDataSet("dataset_mc_failed_nondif_" + suffix,
                             "dataset_mc_failed_nondif_" + suffix, fit_vars);
        }
      }
    }

    int eta_bin = 0, p_bin = 0, ntrks_bin = 0;

    // Loop over MC files first to build combined MC sample,
    // but calculating separate mass window efficiencies
    for (auto year : years) {
      // Open and loop over MC files
      const auto mc_path = ymlConfig["mc_ntps"][year][probe].as<string>();
      cout << "INFO Opening MC files: " << mc_path << endl;
      TChain ch_mc("tree");
      ch_mc.Add(mc_path.c_str());
      cout << "INFO Opened MC files:" << endl;
      print_files(ch_mc);

      cout << "INFO Setting input branches " << endl;
      ch_mc.SetBranchStatus("*", false);
      ch_mc.SetBranchAddress("dst_M", &dst_m);
      ch_mc.SetBranchAddress("dst_TRUEID", &dst_id);
      ch_mc.SetBranchAddress("d0_M", &d0_m);
      ch_mc.SetBranchAddress("d0_PT", &d0_pt);
      ch_mc.SetBranchAddress("d0_TRUEID", &d0_id);
      ch_mc.SetBranchAddress("spi_P", &spi_p);
      ch_mc.SetBranchAddress("spi_PT", &spi_pt);
      ch_mc.SetBranchAddress(tag + "_P", &tag_p);
      ch_mc.SetBranchAddress(tag + "_PT", &tag_pt);
      ch_mc.SetBranchAddress((probe + "_P").c_str(), &probe_p);
      ch_mc.SetBranchAddress((probe + "_PT").c_str(), &probe_pt);
      ch_mc.SetBranchAddress((probe + "_PZ").c_str(), &probe_pz);
      ch_mc.SetBranchAddress((probe + "_isMuon").c_str(), &probe_ismuon);
      ch_mc.SetBranchAddress((probe + "_hasMuon").c_str(), &probe_hasmuon);
      ch_mc.SetBranchAddress((probe + "_PIDmu").c_str(), &probe_dllmu);
      ch_mc.SetBranchAddress((probe + "_PIDe").c_str(), &probe_dlle);
      ch_mc.SetBranchAddress((probe + "_bdt_mu").c_str(), &probe_mu_ubdt);
      ch_mc.SetBranchAddress((probe + "_TRUEID").c_str(), &probe_trueid);
      ch_mc.SetBranchAddress((probe + "_DAUGHTER0_ID").c_str(),
                             &probe_daughter0_trueid);
      ch_mc.SetBranchAddress("k_TRACK_CHI2NDOF", &k_track_chi2ndof);
      ch_mc.SetBranchAddress("pi_TRACK_CHI2NDOF", &pi_track_chi2ndof);
      ch_mc.SetBranchAddress("spi_TRACK_CHI2NDOF", &spi_track_chi2ndof);
      ch_mc.SetBranchAddress("nTracks", &ntracks);

      int count_tm = 0, count_cond = 0, count_calib_sel = 0, count_range = 0,
          count_mw = 0;

      int    kin_bin        = -1;
      bool   in_mass_window = false, pid_ok = false, dif = false;
      double dm = 0.;

      const int entries_mc = ch_mc.GetEntries();
      cout << "INFO Starting MC event loop over " << entries_mc << " entries"
           << endl;

      auto start_mc_loop = high_resolution_clock::now();
      for (int evt = 0; evt < entries_mc; evt++) {
        ch_mc.GetEntry(evt);

        // Offline truth-matching
        if (abs(d0_id) != D0_ID || abs(dst_id) != Dst_ID) continue;
        count_tm++;

        // Conditional cuts
        if (!probe_hasmuon) continue;
        count_cond++;

        // Calib sample cuts
        // See tables 36 and 48 in LHCb-PUB-2016-005
        // See also
        // https://lhcb-pid-wgp-plots.web.cern.ch/lhcb-pid-wgp-plots/Run2/ and
        // https://gitlab.cern.ch/lhcb-datapkg/WG/SemilepConfig/-/blob/master/options/Filter_Dstar2D02KpiNoPID_2016MC.py?ref_type=heads
        if (tag_p < 2000. || probe_p < 2000. || spi_p < 1000.) continue;
        // if (tag_p < 2000. || tag_pt < 250.) continue;
        // if (probe_p < 2000. || probe_pt < 250.) continue;
        // if (spi_p < 1000. || spi_pt < 100.) continue;
        if (k_track_chi2ndof > 3. || pi_track_chi2ndof > 3. ||
            spi_track_chi2ndof > 3.)
          continue;
        // if (std::max(probe_pt, tag_pt) < 1000. || d0_pt < 1500.) continue;
        count_calib_sel++;

        // MuonUnbiased equivalent to L0+HLT1 TIS. See MuonUnbiased defition at
        // https://gitlab.cern.ch/lhcb-datapkg/WG/PIDCalib/-/blob/master/scriptsR2/makeTuples_pp_2016_reprocessing.py#L71
        // and https://mattermost.web.cern.ch/lhcb/pl/893yre484jggigooti5u3gqb8w

        // Not applying MuonUnbiased and GHOSTPROB conditional cuts in MC since
        // some kinematical bins have very low statistics

        if (probe_p < 3000. || probe_p >= 100000. || ntracks >= 600) continue;

        probe_eta = 0.5 * log((probe_p + probe_pz) / (probe_p - probe_pz));
        if (probe_eta < 1.7 || probe_eta >= 5.0) continue;
        count_range++;

        // Determine kinematical bin
        kin_bin = histo_binnning.FindBin(probe_p, probe_eta, ntracks);
        histo_binnning.GetBinXYZ(kin_bin, p_bin, eta_bin, ntrks_bin);

        // Fill histograms
        dm   = (dst_m - d0_m) * 0.001;
        d0_m = d0_m * 0.001;

        in_mass_window = in_range(d0_m_var, d0_m) && in_range(dm_var, dm);
        if (in_mass_window) count_mw++;

        // Check for decay in flight of probe hadron
        dif = (abs(probe_trueid) == MU_ID) ||
              (abs(probe_daughter0_trueid) == MU_ID);

        // PID
        if (fake_mu) {
          pid_ok = !probe_ismuon;
        } else if (vmu) {
          pid_ok = probe_ismuon && probe_dllmu > 2.0 && probe_dlle < 1.0 &&
                   probe_mu_ubdt < 0.25;
        } else {
          pid_ok = probe_ismuon && probe_dllmu > 2.0 && probe_dlle < 1.0 &&
                   probe_mu_ubdt > 0.25;
        }

        if (pid_ok) {
          // Fill "passed" sample
          if (dif) {
            count_total_passed_dif[year_idx.at(year)][ntrks_bin - 1]
                                  [eta_bin - 1][p_bin - 1]++;
            if (in_mass_window) {
              count_in_mw_passed_dif[year_idx.at(year)][ntrks_bin - 1]
                                    [eta_bin - 1][p_bin - 1]++;
              d0_m_var.setVal(d0_m);
              dm_var.setVal(dm);
              datasets_mc_passed_dif[eta_bin - 1][p_bin - 1]->addFast(fit_vars);
            }
          } else {
            count_total_passed_nondif[year_idx.at(year)][ntrks_bin - 1]
                                     [eta_bin - 1][p_bin - 1]++;
            if (in_mass_window) {
              count_in_mw_passed_nondif[year_idx.at(year)][ntrks_bin - 1]
                                       [eta_bin - 1][p_bin - 1]++;
              d0_m_var.setVal(d0_m);
              dm_var.setVal(dm);
              datasets_mc_passed_nondif[ntrks_bin - 1][eta_bin - 1][p_bin - 1]
                  ->addFast(fit_vars);
            }
          }
        } else {
          // Fill "failed" sample
          if (dif) {
            count_total_failed_dif[year_idx.at(year)][ntrks_bin - 1]
                                  [eta_bin - 1][p_bin - 1]++;
            if (in_mass_window) {
              count_in_mw_failed_dif[year_idx.at(year)][ntrks_bin - 1]
                                    [eta_bin - 1][p_bin - 1]++;
              d0_m_var.setVal(d0_m);
              dm_var.setVal(dm);
              datasets_mc_failed_dif[eta_bin - 1][p_bin - 1]->addFast(fit_vars);
            }
          } else {
            count_total_failed_nondif[year_idx.at(year)][ntrks_bin - 1]
                                     [eta_bin - 1][p_bin - 1]++;
            if (in_mass_window) {
              count_in_mw_failed_nondif[year_idx.at(year)][ntrks_bin - 1]
                                       [eta_bin - 1][p_bin - 1]++;
              d0_m_var.setVal(d0_m);
              dm_var.setVal(dm);
              datasets_mc_failed_nondif[ntrks_bin - 1][eta_bin - 1][p_bin - 1]
                  ->addFast(fit_vars);
            }
          }
        }
      }

      auto stop_mc_loop = high_resolution_clock::now();

      const int duration_mc_loop =
          duration_cast<seconds>(stop_mc_loop - start_mc_loop).count();
      cout << "INFO MC loop took " << format_time(duration_mc_loop) << endl;

      cout << "INFO MC cutflow:" << endl;
      cout << " - Total MC events:              " << entries_mc << endl;
      cout << " - After offline truth-matching: " << count_tm << "("
           << count_tm * 100. / entries_mc << "%)" << endl;
      cout << " - After conditional cuts:       " << count_cond << "("
           << count_cond * 100. / count_tm << "%)" << endl;
      cout << " - After PIDCalib cuts:          " << count_calib_sel << "("
           << count_calib_sel * 100. / count_cond << "%)" << endl;
      cout << " - After range cuts:             " << count_range << "("
           << count_range * 100. / count_calib_sel << "%)" << endl;
      cout << " - After mass window cuts:       " << count_mw << "("
           << count_mw * 100. / count_range << "%)" << endl;
      cout << " Overall: " << "(" << count_mw * 100. / entries_mc << "%)\n"
           << endl;
    }

    // Now loop over calib samples and make fits
    for (auto year : years) {
      // TODO Run all years
      if (year != "2016") continue;

      // Define calib datasets
      RooDataSet *datasets_calib_passed[N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] =
          {{{nullptr}}};
      RooDataSet *datasets_calib_failed[N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] =
          {{{nullptr}}};

      // Initialize calib datasets
      cout << "INFO Initializing calib datasets " << endl;
      for (int ntrks_idx = 0; ntrks_idx < N_BINS_NTRACKS; ntrks_idx++) {
        for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
          for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
            TString suffix =
                TString::Format("%s_%s_%d_%d_%d", year.c_str(), probe.c_str(),
                                ntrks_idx, eta_idx, p_idx);
            datasets_calib_passed[ntrks_idx][eta_idx][p_idx] =
                new RooDataSet("dataset_calib_passed_" + suffix,
                               "dataset_calib_passed_" + suffix, fit_vars);
            datasets_calib_failed[ntrks_idx][eta_idx][p_idx] =
                new RooDataSet("dataset_calib_failed_" + suffix,
                               "dataset_calib_failed_" + suffix, fit_vars);
          }
        }
      }

      // Open and loop over PIDCalib files
      const auto tree_names =
          ymlConfig["pid_calib"]["trees"][probe].as<vector<string>>();
      const auto calib_paths =
          ymlConfig["pid_calib"]["calib_samples"][year].as<vector<string>>();
      const auto friends_paths =
          ymlConfig["pid_calib"]["friends"][year].as<vector<string>>();

      for (int p = 0; p < calib_paths.size(); p++) {
        auto calib_path   = calib_paths[p];
        auto friends_path = friends_paths[p];
        cout << "INFO Opening PIDCalib files: " << calib_path << endl;

        for (auto tree_name : tree_names) {
          TChain ch_calib(tree_name.c_str());

          ch_calib.Add(calib_path.c_str());
          cout << "INFO Opened PIDCalib files:" << endl;
          print_files(ch_calib);

          cout << "INFO Opening friend ntuples with uBDT: " << friends_path
               << endl;
          TChain ch_friends(tree_name.c_str());
          ch_friends.Add(friends_path.c_str());
          cout << "INFO Opened friend input files:" << endl;
          print_files(ch_friends);

          // TODO Make sure order is correct
          ch_calib.AddFriend(&ch_friends);

          cout << "INFO Setting input branches " << endl;
          ch_calib.SetBranchStatus("*", 0);
          ch_calib.SetBranchAddress("Dst_M", &dst_m);
          ch_calib.SetBranchAddress("Dz_M", &d0_m);
          ch_calib.SetBranchAddress("probe_Brunel_P", &probe_p);
          ch_calib.SetBranchAddress("probe_Brunel_ETA", &probe_eta);
          ch_calib.SetBranchAddress("probe_Brunel_isMuon", &probe_ismuon);
          ch_calib.SetBranchAddress("probe_Brunel_hasMuon", &probe_hasmuon);
          ch_calib.SetBranchAddress("probe_Brunel_PIDmu", &probe_dllmu);
          ch_calib.SetBranchAddress("probe_Brunel_PIDe", &probe_dlle);
          ch_calib.SetBranchAddress("probe_Brunel_TRACK_GHOSTPROB",
                                    &probe_ghostprob);
          ch_calib.SetBranchAddress("probe_Brunel_MuonUnbiased",
                                    &probe_mu_unbiased);
          ch_calib.SetBranchAddress("probe_UBDT", &probe_mu_ubdt);
          ch_calib.SetBranchAddress("nTracks_Brunel", &ntracks_calib);
          // PID cut for companion particle
          //   if (probe == "pi")
          //     ch_calib.SetBranchAddress("k_Brunel_MC15TuneV1_ProbNNk",
          //                               &tag_probnn);
          //   else
          //     ch_calib.SetBranchAddress("pi_Brunel_MC15TuneV1_ProbNNpi",
          //                               &tag_probnn);

          int count_cond = 0, count_range = 0, count_companion = 0;

          bool pid_ok  = false;
          int  kin_bin = -1;

          const int entries_calib = (dry_run || mc_only)
                                        ? ch_calib.GetEntries() * 0.1
                                        : ch_calib.GetEntries();
          cout << "INFO Starting PIDCalib event loop over " << entries_calib
               << " entries in tree " << tree_name << endl;
          auto start_calib_loop = high_resolution_clock::now();
          for (int evt = 0; evt < entries_calib; evt++) {
            ch_calib.GetEntry(evt);

            // Conditional cuts

            if (!probe_hasmuon || probe_ghostprob > 0.5 || !probe_mu_unbiased)
              continue;
            count_cond++;

            if (probe_p < 3000. || probe_p >= 100000. || ntracks_calib >= 600)
              continue;

            if (probe_eta < 1.7 || probe_eta >= 5.0) continue;
            count_range++;

            // if (tag_probnn > 0.5) continue;
            count_companion++;

            // Fill histograms

            d0_m_var.setVal(d0_m * 0.001);
            dm_var.setVal((dst_m - d0_m) * 0.001);

            // Determine kinematical bin
            kin_bin = histo_binnning.FindBin(probe_p, probe_eta, ntracks_calib);

            histo_binnning.GetBinXYZ(kin_bin, p_bin, eta_bin, ntrks_bin);

            if (fake_mu) {
              pid_ok = !probe_ismuon;
            } else if (vmu) {
              pid_ok = probe_ismuon && probe_dllmu > 2.0 && probe_dlle < 1.0 &&
                       probe_mu_ubdt < 0.25;
            } else {
              pid_ok = probe_ismuon && probe_dllmu > 2.0 && probe_dlle < 1.0 &&
                       probe_mu_ubdt > 0.25;
            }

            if (pid_ok) {
              datasets_calib_passed[ntrks_bin - 1][eta_bin - 1][p_bin - 1]
                  ->addFast(fit_vars);
            } else {
              datasets_calib_failed[ntrks_bin - 1][eta_bin - 1][p_bin - 1]
                  ->addFast(fit_vars);
            }
          }

          auto stop_calib_loop = high_resolution_clock::now();

          const int duration_calib_loop =
              duration_cast<seconds>(stop_calib_loop - start_calib_loop)
                  .count();
          cout << "INFO Calib loop took " << format_time(duration_calib_loop)
               << endl;

          cout << "INFO Calib cutflow:" << endl;
          cout << " - Total MC events:              " << entries_calib << endl;
          cout << " - After conditional cuts:       " << count_cond << "("
               << count_cond * 100. / entries_calib << "%)" << endl;
          cout << " - After range cuts:             " << count_range << "("
               << count_range * 100. / count_cond << "%)" << endl;
          cout << " - After companion cuts:             " << count_companion
               << "(" << count_companion * 100. / count_range << "%)" << endl;
          cout << " Overall: " << "(" << count_companion * 100. / entries_calib
               << "%)\n"
               << endl;
        }
      }

      // Create histograms to store efficiency of PIDCalib mass windows
      cout << "INFO Initializing efficiency histograms " << endl;
      TString suffix = TString::Format("%s_%s", year.c_str(), probe.c_str());

      TH3D effs_passed(histo_binnning);
      effs_passed.SetName("effs_passed_" + suffix);

      TH3D effs_passed_unc_hi(histo_binnning);
      effs_passed_unc_hi.SetName("effs_passed_" + suffix + "_unc_hi");

      TH3D effs_passed_unc_lo(histo_binnning);
      effs_passed_unc_lo.SetName("effs_passed_" + suffix + "_unc_lo");

      TH3D effs_failed(histo_binnning);
      effs_failed.SetName("effs_failed_" + suffix);

      TH3D effs_failed_unc_hi(histo_binnning);
      effs_failed_unc_hi.SetName("effs_failed_" + suffix + "_unc_hi");

      TH3D effs_failed_unc_lo(histo_binnning);
      effs_failed_unc_lo.SetName("effs_failed_" + suffix + "_unc_lo");

      // Check distribution of fitted parameters
      RooDataSet ds_params_d0_dif("ds_params_d0_dif", "ds_params_d0_dif",
                                  params_d0_dif);
      RooDataSet ds_params_d0_nondif("ds_params_d0_nondif",
                                     "ds_params_d0_nondif", params_d0_nondif);
      RooDataSet ds_params_dm_dif("ds_params_dm_dif", "ds_params_dm_dif",
                                  params_dm_dif);
      RooDataSet ds_params_dm_nondif("ds_params_dm_nondif",
                                     "ds_params_dm_nondif", params_dm_nondif);
      RooDataSet ds_params_calib("ds_params_calib", "ds_params_calib",
                                 params_calib);

      TH3D fit_status_calib(histo_binnning);
      fit_status_calib.SetName("fit_status_calib_" + suffix);
      fit_status_calib.SetMinimum(-0.1);
      fit_status_calib.SetMaximum(12.);

      TH3D fit_cov_qual_calib(histo_binnning);
      fit_cov_qual_calib.SetName("fit_cov_qual_calib_" + suffix);
      fit_cov_qual_calib.SetMinimum(-0.1);
      fit_cov_qual_calib.SetMaximum(6.);

      TH3D fit_status_d0_mc_dif(histo_binnning);
      fit_status_d0_mc_dif.SetName("fit_status_d0_mc_dif_" + suffix);
      fit_status_d0_mc_dif.SetMinimum(-0.1);
      fit_status_d0_mc_dif.SetMaximum(12.);

      TH3D fit_status_d0_mc_nondif(histo_binnning);
      fit_status_d0_mc_nondif.SetName("fit_status_d0_mc_nondif_" + suffix);
      fit_status_d0_mc_nondif.SetMinimum(-0.1);
      fit_status_d0_mc_nondif.SetMaximum(12.);

      TH3D fit_cov_qual_d0_mc_dif(histo_binnning);
      fit_cov_qual_d0_mc_dif.SetName("fit_cov_qual_d0_mc_dif_" + suffix);
      fit_cov_qual_d0_mc_dif.SetMinimum(-0.1);
      fit_cov_qual_d0_mc_dif.SetMaximum(6.);

      TH3D fit_cov_qual_d0_mc_nondif(histo_binnning);
      fit_cov_qual_d0_mc_nondif.SetName("fit_cov_qual_d0_mc_nondif_" + suffix);
      fit_cov_qual_d0_mc_nondif.SetMinimum(-0.1);
      fit_cov_qual_d0_mc_nondif.SetMaximum(6.);

      TH3D fit_status_dm_mc_dif(histo_binnning);
      fit_status_dm_mc_dif.SetName("fit_status_dm_mc_dif_" + suffix);
      fit_status_dm_mc_dif.SetMinimum(-0.1);
      fit_status_dm_mc_dif.SetMaximum(12.);

      TH3D fit_status_dm_mc_nondif(histo_binnning);
      fit_status_dm_mc_nondif.SetName("fit_status_dm_mc_nondif_" + suffix);
      fit_status_dm_mc_nondif.SetMinimum(-0.1);
      fit_status_dm_mc_nondif.SetMaximum(12.);

      TH3D fit_cov_qual_dm_mc_dif(histo_binnning);
      fit_cov_qual_dm_mc_dif.SetName("fit_cov_qual_dm_mc_dif_" + suffix);
      fit_cov_qual_dm_mc_dif.SetMinimum(-0.1);
      fit_cov_qual_dm_mc_dif.SetMaximum(6.);

      TH3D fit_cov_qual_dm_mc_nondif(histo_binnning);
      fit_cov_qual_dm_mc_nondif.SetName("fit_cov_qual_dm_mc_nondif_" + suffix);
      fit_cov_qual_dm_mc_nondif.SetMinimum(-0.1);
      fit_cov_qual_dm_mc_nondif.SetMaximum(6.);

      TH3D fit_status(histo_binnning);
      fit_status.SetName("fit_status");
      fit_status.SetMinimum(-0.1);
      fit_status.SetMaximum(12.);

      TH3D fit_cov_qual(histo_binnning);
      fit_cov_qual.SetName("fit_cov_qual");
      fit_cov_qual.SetMinimum(-0.1);
      fit_cov_qual.SetMaximum(6.);

      TH1D d0_retries_dif("d0_retries_dif", "d0_m MC fit #retries;#;;",
                          max_fix_reattempts + 1, 0., max_fix_reattempts + 1.);
      TH1D d0_retries_nondif("d0_retries_nondif", "d0_m MC fit #retries;#;;",
                             max_fix_reattempts + 1, 0.,
                             max_fix_reattempts + 1.);
      TH1D dm_retries_dif("dm_retries_dif", "dm MC fit #retries;#;;",
                          max_fix_reattempts + 1, 0., max_fix_reattempts + 1.);
      TH1D dm_retries_nondif("dm_retries_nondif", "dm MC fit #retries;#;;",
                             max_fix_reattempts + 1, 0.,
                             max_fix_reattempts + 1.);
      TH1D calib_retries("calib_retries", "calib fit #retries;#;;",
                         max_fix_reattempts + 1, 0., max_fix_reattempts + 1.);

      ofile.cd();

      cout << "INFO Calculating efficiencies " << endl;

      // Create histogram to hold PID efficiency
      // Use "eff" as name to match PIDCalib output
      TH3D histo_pid(histo_binnning);
      histo_pid.SetName("eff");

      TH3D histo_pid_raw(histo_binnning);
      histo_pid_raw.SetName("eff_raw");

      // Create histogram to store fitted DiF fractions
      TH3D histo_f_dif(histo_binnning);
      histo_f_dif.SetName("f_dif_" + TString(probe));

      for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
        for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
          for (int ntrks_idx = 0; ntrks_idx < N_BINS_NTRACKS; ntrks_idx++) {
            // Since we combine ntrack bins in the DiF samples, iterating over
            // ntrks_idx in the inner loop allows the DiF fits to be reused

            suffix.Form("%s_%s_%d_%d_%d", year.c_str(), probe.c_str(),
                        ntrks_idx, eta_idx, p_idx);
            const TString tag = TString::Format(
                "(%.0f < nTracks < %.0f, %.1f < #eta < %.1f, %.0f < p < %.0f)",
                BINS_NTRACKS[ntrks_idx], BINS_NTRACKS[ntrks_idx + 1],
                BINS_ETA[eta_idx], BINS_ETA[eta_idx + 1], BINS_P[p_idx],
                BINS_P[p_idx + 1]);

            // Get p, eta bin
            const int kin_bin =
                effs_passed.GetBin(p_idx + 1, eta_idx + 1, ntrks_idx + 1);

            //////////////////////////////////////////////////////////
            // Simultaneous fits to passed and failed distributions //
            //////////////////////////////////////////////////////////

            mean_d0_failed_dif.setConstant(false);
            mean_d0_passed_dif.setConstant(false);
            mean_d0_failed_nondif.setConstant(false);
            mean_d0_passed_nondif.setConstant(false);
            width_d0_dif.setConstant(false);
            width_d0_nondif.setConstant(false);
            width_cb_d0_dif.setConstant(false);
            width_cb_d0_nondif.setConstant(false);
            alpha_L_d0_failed_dif.setConstant(false);
            alpha_L_d0_passed_dif.setConstant(false);
            alpha_L_d0_failed_nondif.setConstant(false);
            alpha_L_d0_passed_nondif.setConstant(false);
            n_L_d0_dif.setConstant(false);
            n_L_d0_nondif.setConstant(false);
            alpha_R_passed_d0_dif.setConstant(false);
            alpha_R_failed_d0_dif.setConstant(false);
            alpha_R_passed_d0_nondif.setConstant(false);
            alpha_R_failed_d0_nondif.setConstant(false);
            n_R_d0_dif.setConstant(false);
            n_R_d0_nondif.setConstant(false);
            f_d0_passed_dif.setConstant(false);
            f_d0_failed_dif.setConstant(false);
            f_d0_passed_nondif.setConstant(false);
            f_d0_failed_nondif.setConstant(false);
            mean_dm_failed_dif.setConstant(false);
            mean_dm_passed_dif.setConstant(false);
            mean_dm_failed_nondif.setConstant(false);
            mean_dm_passed_nondif.setConstant(false);
            width_dm_dif.setConstant(false);
            width_dm_nondif.setConstant(false);
            width_cb_dm_dif.setConstant(false);
            width_cb_dm_nondif.setConstant(false);
            alpha_L_dm_dif.setConstant(false);
            alpha_L_dm_nondif.setConstant(false);
            n_L_dm_dif.setConstant(false);
            n_L_dm_nondif.setConstant(false);
            alpha_R_dm_dif.setConstant(false);
            alpha_R_dm_nondif.setConstant(false);
            n_R_dm_dif.setConstant(false);
            n_R_dm_nondif.setConstant(false);
            f_dm_passed_dif.setConstant(false);
            f_dm_failed_dif.setConstant(false);
            f_dm_passed_nondif.setConstant(false);
            f_dm_failed_nondif.setConstant(false);

            // Build combined dataset
            auto &dataset_mc_passed_dif =
                datasets_mc_passed_dif[eta_idx][p_idx];
            auto &dataset_mc_passed_nondif =
                datasets_mc_passed_nondif[ntrks_idx][eta_idx][p_idx];
            auto &dataset_mc_failed_dif =
                datasets_mc_failed_dif[eta_idx][p_idx];
            auto &dataset_mc_failed_nondif =
                datasets_mc_failed_nondif[ntrks_idx][eta_idx][p_idx];
            auto &dataset_calib_passed =
                datasets_calib_passed[ntrks_idx][eta_idx][p_idx];
            auto &dataset_calib_failed =
                datasets_calib_failed[ntrks_idx][eta_idx][p_idx];

            double       n_calib_passed = dataset_calib_passed->numEntries();
            double       n_calib_failed = dataset_calib_failed->numEntries();
            const double n_mc_guess_passed_dif =
                dataset_mc_passed_dif->numEntries();
            const double n_mc_guess_passed_nondif =
                dataset_mc_passed_nondif->numEntries();
            const double n_mc_guess_failed_dif =
                dataset_mc_failed_dif->numEntries();
            const double n_mc_guess_failed_nondif =
                dataset_mc_failed_nondif->numEntries();

            // Plot 2D distributions in calib sample
            unique_ptr<TH2F> th2_calib_passed(dataset_calib_passed->createHistogram(
                d0_m_var, dm_var, 40, 40, "", "th2_calib_passed"));
            unique_ptr<TH2F> th2_calib_failed(dataset_calib_failed->createHistogram(
                d0_m_var, dm_var, 100, 100, "", "th2_calib_failed"));
            c_single.cd();
            th2_calib_passed->Draw("COLZ");
            c_single.SaveAs(opath + "/figs/calib/th2_" + suffix +
                            "_passed.pdf");
            th2_calib_failed->Draw("COLZ");
            c_single.SaveAs(opath + "/figs/calib/th2_" + suffix +
                            "_failed.pdf");

            cout << "\nINFO Building combined dataset with:" << endl;
            cout << "  - " << n_calib_passed << " passed calib events" << endl;
            cout << "  - " << n_calib_failed << " failed calib events" << endl;
            cout << "  - " << n_mc_guess_passed_dif << " passed MC dif events"
                 << endl;
            cout << "  - " << n_mc_guess_passed_nondif
                 << " passed MC non-dif events" << endl;
            cout << "  - " << n_mc_guess_failed_dif << " failed MC dif events"
                 << endl;
            cout << "  - " << n_mc_guess_failed_nondif
                 << " failed MC non-dif events" << endl;

            RooDataSet dataset(
                "dataset", "combined calib + MC data", fit_vars, Index(sample),
                Import("mc_passed_dif", *dataset_mc_passed_dif),
                Import("mc_passed_nondif", *dataset_mc_passed_nondif),
                Import("mc_failed_dif", *dataset_mc_failed_dif),
                Import("mc_failed_nondif", *dataset_mc_failed_nondif),
                Import("calib_passed", *dataset_calib_passed),
                Import("calib_failed", *dataset_calib_failed));

            cout << "INFO Built combined dataset with " << dataset.numEntries()
                 << " calib + MC events" << endl;

            // Create binnined dataset to calculate fit chi2ndof
            // FIXME ndof partially hard-coded
            d0_m_var.setBins(50);
            dm_var.setBins(50);
            unique_ptr<RooDataHist> datahist_mc_passed_dif(
                dataset_mc_passed_dif->binnedClone(nullptr, nullptr));
            unique_ptr<RooDataHist> datahist_mc_passed_nondif(
                dataset_mc_passed_nondif->binnedClone(nullptr, nullptr));
            unique_ptr<RooDataHist> datahist_mc_failed_dif(
                dataset_mc_failed_dif->binnedClone(nullptr, nullptr));
            unique_ptr<RooDataHist> datahist_mc_failed_nondif(
                dataset_mc_failed_nondif->binnedClone(nullptr, nullptr));
            unique_ptr<RooDataHist> datahist_calib_passed(
                dataset_calib_passed->binnedClone(nullptr, nullptr));
            unique_ptr<RooDataHist> datahist_calib_failed(
                dataset_calib_failed->binnedClone(nullptr, nullptr));

            RooDataHist datahist_mc_dif_d0_m(
                "datahist_mc_dif_d0_m", "MC DiF datahist", d0_m_var,
                Index(sample), Import("mc_passed_dif", *datahist_mc_passed_dif),
                Import("mc_failed_dif", *datahist_mc_failed_dif));

            RooDataHist datahist_mc_dif_dm(
                "datahist_mc_dif_dm", "MC DiF datahist", dm_var, Index(sample),
                Import("mc_passed_dif", *datahist_mc_passed_dif),
                Import("mc_failed_dif", *datahist_mc_failed_dif));

            RooDataHist datahist_mc_nondif_d0_m(
                "datahist_mc_nondif_d0_m", "MC non-DiF datahist", d0_m_var,
                Index(sample),
                Import("mc_passed_nondif", *datahist_mc_passed_nondif),
                Import("mc_failed_nondif", *datahist_mc_failed_nondif));

            RooDataHist datahist_mc_nondif_dm(
                "datahist_mc_nondif_dm", "MC non-DiF datahist", dm_var,
                Index(sample),
                Import("mc_passed_nondif", *datahist_mc_passed_nondif),
                Import("mc_failed_nondif", *datahist_mc_failed_nondif));

            RooDataHist datahist_calib(
                "datahist_calib", "calib datahist", fit_vars, Index(sample),
                Import("calib_passed", *datahist_calib_passed),
                Import("calib_failed", *datahist_calib_failed));

            // d0_m fit without hadron decay in flight
            // With "migradimproved", if MIGRAD finds a minimum, IMPROVE is
            // called to try to find a better one.
            cout << "\nINFO Starting D0 MC fit without dif " << suffix << endl;

            if (!dry_run) {
              auto start_mc_d0_nondif = high_resolution_clock::now();

              set_parameters_d0_nondif(probe);
              RooFitResult *fit_res_d0_mc_nondif = d0_model_mc_nondif.fitTo(
                  dataset, Strategy(2), NumCPU(8), PrintLevel(0), Save(),
                  Offset(), Minimizer("Minuit", "migradimproved"), Timer());

              cout << "\nINFO Fit status: " << fit_res_d0_mc_nondif->status()
                   << endl;
              cout << "INFO Covariance matrix status: "
                   << fit_res_d0_mc_nondif->covQual() << "\n"
                   << endl;

              int d0_nondif_fit_reattempts = 0;
              // If the fit does not converge, mimic the logic of "minimize" and
              // keep alternating between "simplex" and "migradimproved" until
              // it converges or hits the limit of retries. Here we hope that
              // SIMPLEX can at least find a better starting point for MIGRAD.
              while (!fit_ok(fit_res_d0_mc_nondif) &&
                     (d0_nondif_fit_reattempts < max_fix_reattempts)) {
                d0_nondif_fit_reattempts++;
                cout << "INFO Retrying D0 MC non-dif fit " << suffix
                     << " (retry #" << d0_nondif_fit_reattempts << ")" << endl;

                delete fit_res_d0_mc_nondif;

                d0_model_mc_nondif.fitTo(
                    dataset, NumCPU(8), PrintLevel(0), Offset(),
                    Minimizer("Minuit", "simplex"), Hesse(false), Timer());
                fit_res_d0_mc_nondif = d0_model_mc_nondif.fitTo(
                    dataset, Strategy(2), NumCPU(8), PrintLevel(0), Save(),
                    Offset(), Minimizer("Minuit", "migradimproved"), Timer());

                cout << "\nINFO Fit status: " << fit_res_d0_mc_nondif->status()
                     << endl;
                cout << "INFO Covariance matrix status: "
                     << fit_res_d0_mc_nondif->covQual() << "\n"
                     << endl;
              }

              auto stop_mc_d0_nondif = high_resolution_clock::now();

              const int duration_mc_d0_nondif =
                  duration_cast<seconds>(stop_mc_d0_nondif - start_mc_d0_nondif)
                      .count();

              cout << "INFO D0 non-dif MC fits took "
                   << format_time(duration_mc_d0_nondif) << endl;
              cout << "INFO Needed " << d0_nondif_fit_reattempts
                   << " retries with "
                   << format_time(duration_mc_d0_nondif /
                                  (d0_nondif_fit_reattempts + 1))
                   << " on average\n"
                   << endl;

              // Get reduced chi2
              RooChi2Var chi2_d0_mc_nondif(
                  "chi2_d0_mc_nondif", "chi2_d0_mc_nondif", d0_model_mc_nondif,
                  datahist_mc_nondif_d0_m, Extended(true),
                  DataError(RooAbsData::Poisson));
              // FIXME ndof depends on hard-coded number of bins in roodatahist
              // FIXME chi2ndof not correct
              const int ndof_d0_mc_nondif =
                  2 * 50 - fit_res_d0_mc_nondif->floatParsFinal().size();
              const double chi2ndof_d0_mc_nondif =
                  chi2_d0_mc_nondif.getVal() / ndof_d0_mc_nondif;

              cout << "INFO Finished D0 non-dif MC fit " << suffix
                   << " with chi2ndof = " << chi2ndof_d0_mc_nondif << "\n"
                   << endl;

              d0_retries_nondif.Fill(d0_nondif_fit_reattempts);

              if (fit_ok(fit_res_d0_mc_nondif)) {
                ds_params_d0_nondif.addFast(params_d0_nondif);
              }

              fit_status_d0_mc_nondif.SetBinContent(
                  kin_bin, fit_res_d0_mc_nondif->status());
              fit_cov_qual_d0_mc_nondif.SetBinContent(
                  kin_bin, fit_res_d0_mc_nondif->covQual());

              delete fit_res_d0_mc_nondif;
            }

            // Plot fit results
            unique_ptr<RooPlot> frame_d0_passed_nondif(
                d0_m_var.frame(Title("D0 M Passed " + tag)));
            unique_ptr<RooPlot> frame_d0_failed_nondif(
                d0_m_var.frame(Title("D0 M Failed " + tag)));
            dataset.plotOn(frame_d0_passed_nondif.get(),
                           Binning(bins_d0_m_passed),
                           Cut("sample==sample::mc_passed_nondif"));
            dataset.plotOn(frame_d0_failed_nondif.get(), Binning(bins_d0_m),
                           Cut("sample==sample::mc_failed_nondif"));
            d0_model_mc_nondif.plotOn(frame_d0_passed_nondif.get(),
                                      Slice(sample, "mc_passed_nondif"),
                                      ProjWData(sample, dataset));
            d0_model_mc_nondif.plotOn(
                frame_d0_passed_nondif.get(), Slice(sample, "mc_passed_nondif"),
                ProjWData(sample, dataset),
                Components(d0_gauss_passed_mc_nondif), LineColor(kRed));
            d0_model_mc_nondif.plotOn(
                frame_d0_passed_nondif.get(), Slice(sample, "mc_passed_nondif"),
                ProjWData(sample, dataset), Components(d0_cb_passed_mc_nondif),
                LineColor(kMagenta));
            d0_model_mc_nondif.plotOn(frame_d0_failed_nondif.get(),
                                      Slice(sample, "mc_failed_nondif"),
                                      ProjWData(sample, dataset));
            d0_model_mc_nondif.plotOn(
                frame_d0_failed_nondif.get(), Slice(sample, "mc_failed_nondif"),
                ProjWData(sample, dataset),
                Components(d0_gauss_failed_mc_nondif), LineColor(kRed));
            d0_model_mc_nondif.plotOn(
                frame_d0_failed_nondif.get(), Slice(sample, "mc_failed_nondif"),
                ProjWData(sample, dataset), Components(d0_cb_failed_mc_nondif),
                LineColor(kMagenta));

            c_double.cd(1);
            frame_d0_passed_nondif->Draw();
            c_double.cd(2);
            frame_d0_failed_nondif->Draw();
            c_double.SaveAs(opath + "/figs/fits/d0_m_" + suffix +
                            "_nondif.pdf");

            // d0_m fit with hadron decay in flight
            cout << "\nINFO Starting D0 MC fit with dif " << suffix << endl;

            if (!dry_run && ntrks_idx == 0) {
              auto start_mc_d0_dif = high_resolution_clock::now();

              set_parameters_d0_dif(probe);
              // Tuning for dif D0 MC fit
              if (probe == "pi" && eta_idx == 0 && p_idx == 0) {
                f_d0_passed_dif.setConstant();
                f_d0_passed_dif.setVal(0.);
              }
              if (probe == "pi" && eta_idx == 1 && p_idx == 1) {
                f_d0_passed_dif.setConstant();
                f_d0_passed_dif.setVal(0.);
              }

              RooFitResult *fit_res_d0_mc_dif = d0_model_mc_dif.fitTo(
                  dataset, Strategy(2), NumCPU(8), PrintLevel(0), Save(),
                  Offset(), Minimizer("Minuit", "migradimproved"), Timer());

              cout << "\nINFO Fit status: " << fit_res_d0_mc_dif->status()
                   << endl;
              cout << "INFO Covariance matrix status: "
                   << fit_res_d0_mc_dif->covQual() << "\n"
                   << endl;

              int d0_dif_fit_reattempts = 0;
              while (!fit_ok(fit_res_d0_mc_dif) &&
                     (d0_dif_fit_reattempts < max_fix_reattempts)) {
                d0_dif_fit_reattempts++;
                cout << "INFO Retrying D0 MC dif fit " << suffix << " (retry #"
                     << d0_dif_fit_reattempts << ")" << endl;

                delete fit_res_d0_mc_dif;

                d0_model_mc_dif.fitTo(dataset, NumCPU(8), PrintLevel(0),
                                      Offset(), Minimizer("Minuit", "simplex"),
                                      Hesse(false), Timer());
                fit_res_d0_mc_dif = d0_model_mc_dif.fitTo(
                    dataset, Strategy(2), NumCPU(8), PrintLevel(0), Save(),
                    Offset(), Minimizer("Minuit", "migradimproved"), Timer());

                cout << "\nINFO Fit status: " << fit_res_d0_mc_dif->status()
                     << endl;
                cout << "INFO Covariance matrix status: "
                     << fit_res_d0_mc_dif->covQual() << "\n"
                     << endl;
              }

              auto stop_mc_d0_dif = high_resolution_clock::now();

              const int duration_mc_d0_dif =
                  duration_cast<seconds>(stop_mc_d0_dif - start_mc_d0_dif)
                      .count();

              cout << "INFO D0 dif MC fits took "
                   << format_time(duration_mc_d0_dif) << endl;
              cout << "INFO Needed " << d0_dif_fit_reattempts
                   << " retries with "
                   << format_time(duration_mc_d0_dif /
                                  (d0_dif_fit_reattempts + 1))
                   << " on average\n"
                   << endl;

              // Get reduced chi2
              RooChi2Var chi2_d0_mc_dif("chi2_d0_mc_dif", "chi2_d0_mc_dif",
                                        d0_model_mc_dif, datahist_mc_dif_d0_m,
                                        Extended(true),
                                        DataError(RooAbsData::Poisson));
              // FIXME ndof depends on hard-coded number of bins in roodatahist
              const int ndof_d0_mc_dif =
                  2 * 50 - fit_res_d0_mc_dif->floatParsFinal().size();
              const double chi2ndof_d0_mc_dif =
                  chi2_d0_mc_dif.getVal() / ndof_d0_mc_dif;

              cout << "INFO Finished D0 dif MC fit " << suffix
                   << " with chi2ndof = " << chi2ndof_d0_mc_dif << "\n"
                   << endl;

              d0_retries_dif.Fill(d0_dif_fit_reattempts);

              if (fit_ok(fit_res_d0_mc_dif)) {
                ds_params_d0_dif.addFast(params_d0_dif);
              }

              fit_status_d0_mc_dif.SetBinContent(kin_bin,
                                                 fit_res_d0_mc_dif->status());
              fit_cov_qual_d0_mc_dif.SetBinContent(
                  kin_bin, fit_res_d0_mc_dif->covQual());

              delete fit_res_d0_mc_dif;
            }

            // Plot fit results
            unique_ptr<RooPlot> frame_d0_passed_dif(
                d0_m_var.frame(Title("D0 M Passed " + tag)));
            unique_ptr<RooPlot> frame_d0_failed_dif(
                d0_m_var.frame(Title("D0 M Failed " + tag)));
            dataset.plotOn(frame_d0_passed_dif.get(), Binning(bins_d0_m_passed),
                           Cut("sample==sample::mc_passed_dif"));
            dataset.plotOn(frame_d0_failed_dif.get(), Binning(bins_d0_m),
                           Cut("sample==sample::mc_failed_dif"));
            d0_model_mc_dif.plotOn(frame_d0_passed_dif.get(),
                                   Slice(sample, "mc_passed_dif"),
                                   ProjWData(sample, dataset));
            d0_model_mc_dif.plotOn(
                frame_d0_passed_dif.get(), Slice(sample, "mc_passed_dif"),
                ProjWData(sample, dataset), Components(d0_gauss_passed_mc_dif),
                LineColor(kRed));
            d0_model_mc_dif.plotOn(
                frame_d0_passed_dif.get(), Slice(sample, "mc_passed_dif"),
                ProjWData(sample, dataset), Components(d0_cb_passed_mc_dif),
                LineColor(kMagenta));
            d0_model_mc_dif.plotOn(frame_d0_failed_dif.get(),
                                   Slice(sample, "mc_failed_dif"),
                                   ProjWData(sample, dataset));
            d0_model_mc_dif.plotOn(
                frame_d0_failed_dif.get(), Slice(sample, "mc_failed_dif"),
                ProjWData(sample, dataset), Components(d0_gauss_failed_mc_dif),
                LineColor(kRed));
            d0_model_mc_dif.plotOn(
                frame_d0_failed_dif.get(), Slice(sample, "mc_failed_dif"),
                ProjWData(sample, dataset), Components(d0_cb_failed_mc_dif),
                LineColor(kMagenta));

            c_double.cd(1);
            frame_d0_passed_dif->Draw();
            c_double.cd(2);
            frame_d0_failed_dif->Draw();
            c_double.SaveAs(opath + "/figs/fits/d0_m_" + suffix + "_dif.pdf");

            // dm fit without hadron dif
            cout << "\nINFO Starting dm MC fit without dif " << suffix << endl;

            if (!dry_run) {
              auto start_mc_dm_nondif = high_resolution_clock::now();

              set_parameters_dm_nondif(probe);
              // Tuning for nondif dm MC fit
              if (probe == "pi" && eta_idx == 0 && p_idx == 0) {
                f_d0_passed_dif.setConstant();
                f_d0_passed_dif.setVal(0.);
              }
              RooFitResult *fit_res_dm_mc_nondif = dm_model_mc_nondif.fitTo(
                  dataset, Strategy(2), NumCPU(8), Save(), PrintLevel(0),
                  Offset(), Minimizer("Minuit", "migradimproved"), Timer());

              cout << "\nINFO Fit status: " << fit_res_dm_mc_nondif->status()
                   << endl;
              cout << "INFO Covariance matrix status: "
                   << fit_res_dm_mc_nondif->covQual() << "\n"
                   << endl;

              int dm_nondif_fit_reattempts = 0;
              while (!fit_ok(fit_res_dm_mc_nondif) &&
                     (dm_nondif_fit_reattempts < max_fix_reattempts)) {
                dm_nondif_fit_reattempts++;
                cout << "INFO Retrying dm MC non-dif fit " << suffix
                     << " (retry #" << dm_nondif_fit_reattempts << ")" << endl;

                delete fit_res_dm_mc_nondif;

                dm_model_mc_nondif.fitTo(
                    dataset, NumCPU(8), PrintLevel(0), Offset(),
                    Minimizer("Minuit", "simplex"), Hesse(false), Timer());
                fit_res_dm_mc_nondif = dm_model_mc_nondif.fitTo(
                    dataset, Strategy(2), NumCPU(8), PrintLevel(0), Save(),
                    Offset(), Minimizer("Minuit", "migradimproved"), Timer());

                cout << "\nINFO Fit status: " << fit_res_dm_mc_nondif->status()
                     << endl;
                cout << "INFO Covariance matrix status: "
                     << fit_res_dm_mc_nondif->covQual() << "\n"
                     << endl;
              }

              auto stop_mc_dm_nondif = high_resolution_clock::now();

              const int duration_mc_dm_nondif =
                  duration_cast<seconds>(stop_mc_dm_nondif - start_mc_dm_nondif)
                      .count();

              cout << "INFO dm non-dif MC fits took "
                   << format_time(duration_mc_dm_nondif) << endl;
              cout << "INFO Needed " << dm_nondif_fit_reattempts
                   << " retries with "
                   << format_time(duration_mc_dm_nondif /
                                  (dm_nondif_fit_reattempts + 1))
                   << " on average\n"
                   << endl;

              // Get reduced chi2
              RooChi2Var chi2_dm_mc_nondif(
                  "chi2_dm_mc_nondif", "chi2_dm_mc_nondif", dm_model_mc_nondif,
                  datahist_mc_nondif_dm, Extended(true),
                  DataError(RooAbsData::Poisson));
              // FIXME ndof depends on hard-coded number of bins in roodatahist
              const int ndof_dm_mc_nondif =
                  2 * 50 - fit_res_dm_mc_nondif->floatParsFinal().size();
              const double chi2ndof_dm_mc_nondif =
                  chi2_dm_mc_nondif.getVal() / ndof_dm_mc_nondif;

              cout << "INFO Finished dm non-dif MC fit " << suffix
                   << " with chi2ndof = " << chi2ndof_dm_mc_nondif << "\n"
                   << endl;

              dm_retries_nondif.Fill(dm_nondif_fit_reattempts);

              if (fit_ok(fit_res_dm_mc_nondif)) {
                ds_params_dm_nondif.addFast(params_dm_nondif);
              }

              fit_status_dm_mc_nondif.SetBinContent(
                  kin_bin, fit_res_dm_mc_nondif->status());
              fit_cov_qual_dm_mc_nondif.SetBinContent(
                  kin_bin, fit_res_dm_mc_nondif->covQual());

              delete fit_res_dm_mc_nondif;
            }

            // Plot fit results
            unique_ptr<RooPlot> frame_dm_passed_nondif(
                dm_var.frame(Title("dm Passed " + tag)));
            unique_ptr<RooPlot> frame_dm_failed_nondif(
                dm_var.frame(Title("dm Failed " + tag)));
            dataset.plotOn(frame_dm_passed_nondif.get(),
                           Binning(bins_dm_passed),
                           Cut("sample==sample::mc_passed_nondif"));
            dataset.plotOn(frame_dm_failed_nondif.get(), Binning(bins_dm),
                           Cut("sample==sample::mc_failed_nondif"));
            dm_model_mc_nondif.plotOn(frame_dm_passed_nondif.get(),
                                      Slice(sample, "mc_passed_nondif"),
                                      ProjWData(sample, dataset));
            dm_model_mc_nondif.plotOn(
                frame_dm_passed_nondif.get(), Slice(sample, "mc_passed_nondif"),
                ProjWData(sample, dataset),
                Components(dm_gauss_mc_passed_nondif), LineColor(kRed));
            dm_model_mc_nondif.plotOn(
                frame_dm_passed_nondif.get(), Slice(sample, "mc_passed_nondif"),
                ProjWData(sample, dataset), Components(dm_cb_mc_passed_nondif),
                LineColor(kMagenta));
            dm_model_mc_nondif.plotOn(frame_dm_failed_nondif.get(),
                                      Slice(sample, "mc_failed_nondif"),
                                      ProjWData(sample, dataset));
            dm_model_mc_nondif.plotOn(
                frame_dm_failed_nondif.get(), Slice(sample, "mc_failed_nondif"),
                ProjWData(sample, dataset),
                Components(dm_gauss_mc_failed_nondif), LineColor(kRed));
            dm_model_mc_nondif.plotOn(
                frame_dm_failed_nondif.get(), Slice(sample, "mc_failed_nondif"),
                ProjWData(sample, dataset), Components(dm_cb_mc_failed_nondif),
                LineColor(kMagenta));

            c_double.cd(1);
            frame_dm_passed_nondif->Draw();
            c_double.cd(2);
            frame_dm_failed_nondif->Draw();
            c_double.SaveAs(opath + "/figs/fits/dm_" + suffix + "_nondif.pdf");

            // dm fit with hadron dif
            cout << "\nINFO Starting dm MC fit with dif " << suffix << endl;

            if (!dry_run && ntrks_idx == 0) {
              auto start_mc_dm_dif = high_resolution_clock::now();

              set_parameters_dm_dif(probe);
              RooFitResult *fit_res_dm_mc_dif = dm_model_mc_dif.fitTo(
                  dataset, Strategy(2), NumCPU(8), Save(), PrintLevel(0),
                  Offset(), Minimizer("Minuit", "migradimproved"), Timer());

              cout << "\nINFO Fit status: " << fit_res_dm_mc_dif->status()
                   << endl;
              cout << "INFO Covariance matrix status: "
                   << fit_res_dm_mc_dif->covQual() << "\n"
                   << endl;

              int dm_dif_fit_reattempts = 0;
              while (!fit_ok(fit_res_dm_mc_dif) &&
                     (dm_dif_fit_reattempts < max_fix_reattempts)) {
                dm_dif_fit_reattempts++;
                cout << "INFO Retrying dm MC dif fit " << suffix << " (retry #"
                     << dm_dif_fit_reattempts << ")" << endl;

                delete fit_res_dm_mc_dif;

                dm_model_mc_dif.fitTo(dataset, NumCPU(8), PrintLevel(0),
                                      Offset(), Minimizer("Minuit", "simplex"),
                                      Hesse(false), Timer());
                fit_res_dm_mc_dif = dm_model_mc_dif.fitTo(
                    dataset, Strategy(2), NumCPU(8), PrintLevel(0), Save(),
                    Offset(), Minimizer("Minuit", "migradimproved"), Timer());

                cout << "\nINFO Fit status: " << fit_res_dm_mc_dif->status()
                     << endl;
                cout << "INFO Covariance matrix status: "
                     << fit_res_dm_mc_dif->covQual() << "\n"
                     << endl;
              }

              auto stop_mc_dm_dif = high_resolution_clock::now();

              const int duration_mc_dm_dif =
                  duration_cast<seconds>(stop_mc_dm_dif - start_mc_dm_dif)
                      .count();

              cout << "INFO dm dif MC fits took "
                   << format_time(duration_mc_dm_dif) << endl;
              cout << "INFO Needed " << dm_dif_fit_reattempts
                   << " retries with "
                   << format_time(duration_mc_dm_dif /
                                  (dm_dif_fit_reattempts + 1))
                   << " on average\n"
                   << endl;

              // Get reduced chi2
              RooChi2Var chi2_dm_mc_dif(
                  "chi2_dm_mc_dif", "chi2_dm_mc_dif", dm_model_mc_dif,
                  datahist_mc_dif_dm, Extended(true),
                  DataError(RooAbsData::Poisson), IntegrateBins(1e-6));
              // FIXME ndof depends on hard-coded number of bins in roodatahist
              const int ndof_dm_mc_dif =
                  2 * 50 - fit_res_dm_mc_dif->floatParsFinal().size();
              const double chi2ndof_dm_mc_dif =
                  chi2_dm_mc_dif.getVal() / ndof_dm_mc_dif;

              cout << "INFO Finished dm dif MC fit " << suffix
                   << " with chi2ndof = " << chi2ndof_dm_mc_dif << "\n"
                   << endl;

              dm_retries_dif.Fill(dm_dif_fit_reattempts);

              if (fit_ok(fit_res_dm_mc_dif)) {
                ds_params_dm_dif.addFast(params_dm_dif);
              }

              fit_status_dm_mc_dif.SetBinContent(kin_bin,
                                                 fit_res_dm_mc_dif->status());
              fit_cov_qual_dm_mc_dif.SetBinContent(
                  kin_bin, fit_res_dm_mc_dif->covQual());

              delete fit_res_dm_mc_dif;
            }

            // Plot fit results
            unique_ptr<RooPlot> frame_dm_passed_dif(
                dm_var.frame(Title("dm Passed " + tag)));
            unique_ptr<RooPlot> frame_dm_failed_dif(
                dm_var.frame(Title("dm Failed " + tag)));
            dataset.plotOn(frame_dm_passed_dif.get(), Binning(bins_dm_passed),
                           Cut("sample==sample::mc_passed_dif"));
            dataset.plotOn(frame_dm_failed_dif.get(), Binning(bins_dm),
                           Cut("sample==sample::mc_failed_dif"));
            dm_model_mc_dif.plotOn(frame_dm_passed_dif.get(),
                                   Slice(sample, "mc_passed_dif"),
                                   ProjWData(sample, dataset));
            dm_model_mc_dif.plotOn(
                frame_dm_passed_dif.get(), Slice(sample, "mc_passed_dif"),
                ProjWData(sample, dataset), Components(dm_gauss_mc_passed_dif),
                LineColor(kRed));
            dm_model_mc_dif.plotOn(
                frame_dm_passed_dif.get(), Slice(sample, "mc_passed_dif"),
                ProjWData(sample, dataset), Components(dm_cb_mc_passed_dif),
                LineColor(kMagenta));
            dm_model_mc_dif.plotOn(frame_dm_failed_dif.get(),
                                   Slice(sample, "mc_failed_dif"),
                                   ProjWData(sample, dataset));
            dm_model_mc_dif.plotOn(
                frame_dm_failed_dif.get(), Slice(sample, "mc_failed_dif"),
                ProjWData(sample, dataset), Components(dm_gauss_mc_failed_dif),
                LineColor(kRed));
            dm_model_mc_dif.plotOn(
                frame_dm_failed_dif.get(), Slice(sample, "mc_failed_dif"),
                ProjWData(sample, dataset), Components(dm_cb_mc_failed_dif),
                LineColor(kMagenta));

            c_double.cd(1);
            frame_dm_passed_dif->Draw();
            c_double.cd(2);
            frame_dm_failed_dif->Draw();
            c_double.SaveAs(opath + "/figs/fits/dm_" + suffix + "_dif.pdf");

            mean_d0_failed_dif.setConstant();
            mean_d0_passed_dif.setConstant();
            mean_d0_failed_nondif.setConstant();
            mean_d0_passed_nondif.setConstant();
            width_d0_dif.setConstant();
            width_d0_nondif.setConstant();
            width_cb_d0_dif.setConstant();
            width_cb_d0_nondif.setConstant();
            alpha_L_d0_failed_dif.setConstant();
            alpha_L_d0_passed_dif.setConstant();
            alpha_L_d0_failed_nondif.setConstant();
            alpha_L_d0_passed_nondif.setConstant();
            n_L_d0_dif.setConstant();
            n_L_d0_nondif.setConstant();
            alpha_R_passed_d0_dif.setConstant();
            alpha_R_failed_d0_dif.setConstant();
            alpha_R_passed_d0_nondif.setConstant();
            alpha_R_failed_d0_nondif.setConstant();
            n_R_d0_dif.setConstant();
            n_R_d0_nondif.setConstant();
            f_d0_passed_dif.setConstant();
            f_d0_failed_dif.setConstant();
            f_d0_passed_nondif.setConstant();
            f_d0_failed_nondif.setConstant();
            mean_dm_failed_dif.setConstant();
            mean_dm_passed_dif.setConstant();
            mean_dm_failed_nondif.setConstant();
            mean_dm_passed_nondif.setConstant();
            width_dm_dif.setConstant();
            width_dm_nondif.setConstant();
            width_cb_dm_dif.setConstant();
            width_cb_dm_nondif.setConstant();
            alpha_L_dm_dif.setConstant();
            alpha_L_dm_nondif.setConstant();
            n_L_dm_dif.setConstant();
            n_L_dm_nondif.setConstant();
            alpha_R_dm_dif.setConstant();
            alpha_R_dm_nondif.setConstant();
            n_R_dm_dif.setConstant();
            n_R_dm_nondif.setConstant();
            f_dm_passed_dif.setConstant();
            f_dm_failed_dif.setConstant();
            f_dm_passed_nondif.setConstant();
            f_dm_failed_nondif.setConstant();

            // Fit only calib sample

            // Set mass-window and total yields in MC for this kinematic bin.
            // Together with the fitted alpha parameter, this allows the
            // mass-window efficiencies to be calculated.
            n_inmw_passed_dif.setVal(count_in_mw_passed_dif[year_idx.at(
                year)][ntrks_idx][eta_idx][p_idx]);
            n_inmw_passed_nondif.setVal(count_in_mw_passed_nondif[year_idx.at(
                year)][ntrks_idx][eta_idx][p_idx]);
            n_inmw_failed_dif.setVal(count_in_mw_failed_dif[year_idx.at(
                year)][ntrks_idx][eta_idx][p_idx]);
            n_inmw_failed_nondif.setVal(count_in_mw_failed_nondif[year_idx.at(
                year)][ntrks_idx][eta_idx][p_idx]);
            n_total_passed_dif.setVal(count_total_passed_dif[year_idx.at(
                year)][ntrks_idx][eta_idx][p_idx]);
            n_total_passed_nondif.setVal(count_total_passed_nondif[year_idx.at(
                year)][ntrks_idx][eta_idx][p_idx]);
            n_total_failed_dif.setVal(count_total_failed_dif[year_idx.at(
                year)][ntrks_idx][eta_idx][p_idx]);
            n_total_failed_nondif.setVal(count_total_failed_nondif[year_idx.at(
                year)][ntrks_idx][eta_idx][p_idx]);

            // Roughly estimate amount of background
            const TString sb_cut =
                "(dm_var > 0.149) || (dm_var < 0.143) || (d0_m_var > 1.890) || "
                "(d0_m_var < 1.840)";
            constexpr double A_sig = (1.890 - 1.840) * (0.149 - 0.143);
            constexpr double A_fit = (1.910 - 1.825) * (0.153 - 0.141);
            constexpr double A_sb  = A_fit - A_sig;
            const double     n_sb_calib_passed =
                dataset_calib_passed->sumEntries(sb_cut);
            const double n_sb_calib_failed =
                dataset_calib_failed->sumEntries(sb_cut);
            const double n_bkg_calib_guess_passed =
                n_sb_calib_passed * A_fit / A_sb;
            const double n_bkg_calib_guess_failed =
                n_sb_calib_failed * A_fit / A_sb;
            const double f_guess_passed =
                1. - n_bkg_calib_guess_passed / n_calib_passed;
            const double f_guess_failed =
                1. - n_bkg_calib_guess_failed / n_calib_failed;
            const TString sb1_cut =
                "((dm_var > 0.149) || (dm_var < 0.143)) && ((d0_m_var < 1.890) "
                "&& "
                "(d0_m_var > 1.840))";
            const TString sb2_cut =
                "((dm_var < 0.149) && (dm_var > 0.143)) && ((d0_m_var > 1.890) "
                "|| "
                "(d0_m_var < 1.840))";
            const TString sb3_cut =
                "((dm_var > 0.149) || (dm_var < 0.143)) && ((d0_m_var > 1.890) "
                "|| "
                "(d0_m_var < 1.840))";
            const double n_sb1_calib_passed =
                dataset_calib_passed->sumEntries(sb1_cut);
            const double n_sb1_calib_failed =
                dataset_calib_failed->sumEntries(sb1_cut);
            const double n_sb2_calib_passed =
                dataset_calib_passed->sumEntries(sb2_cut);
            const double n_sb2_calib_failed =
                dataset_calib_failed->sumEntries(sb2_cut);
            const double n_sb3_calib_passed =
                dataset_calib_passed->sumEntries(sb3_cut);
            const double n_sb3_calib_failed =
                dataset_calib_failed->sumEntries(sb3_cut);

            constexpr double A_sb3 = (0.153 - 0.149) * (1.840 - 1.825) +
                                     (0.153 - 0.149) * (1.910 - 1.890) +
                                     (0.143 - 0.141) * (1.840 - 1.825) +
                                     (0.143 - 0.141) * (1.910 - 1.890);
            constexpr double A_sb2 = (0.149 - 0.143) * (1.840 - 1.825) +
                                     (0.149 - 0.143) * (1.910 - 1.890);
            constexpr double A_sb1 = (0.153 - 0.149) * (1.890 - 1.840) +
                                     (0.143 - 0.141) * (1.890 - 1.840);

            // FIXME f_sb{1,2} area actually overestimated. For a proper
            // estimate, I would need to subtract the contribution of bkg 3 from
            // sidebands 1 and 2. However, the fit converges as it is so it's
            // not a priority.
            const double f_sb1_calib_passed =
                (n_sb1_calib_passed - A_sb1 / A_sb3 * n_sb3_calib_passed) /
                n_sb_calib_passed;
            const double f_sb1_calib_failed =
                (n_sb1_calib_failed - A_sb1 / A_sb3 * n_sb3_calib_failed) /
                n_sb_calib_failed;
            const double f_sb2_calib_passed =
                (n_sb2_calib_passed - A_sb2 / A_sb3 * n_sb3_calib_passed) /
                n_sb_calib_passed;
            const double f_sb2_calib_failed =
                (n_sb2_calib_failed - A_sb2 / A_sb3 * n_sb3_calib_failed) /
                n_sb_calib_failed;
            const double f_sb3_calib_passed =
                n_sb3_calib_passed / n_sb_calib_passed;
            const double f_sb3_calib_failed =
                n_sb3_calib_failed / n_sb_calib_failed;

            const double c1_guess_passed = f_sb1_calib_passed;
            const double c1_guess_failed = f_sb1_calib_failed;
            const double c2_guess_passed =
                f_sb2_calib_passed / (1 - f_sb1_calib_passed);
            const double c2_guess_failed =
                f_sb2_calib_failed / (1 - f_sb1_calib_failed);

            f1_bkg_failed.setVal(c1_guess_failed);
            f1_bkg_passed.setVal(c1_guess_passed);
            f2_bkg_failed.setVal(c2_guess_failed);
            f2_bkg_passed.setVal(c2_guess_passed);

            // Guess normalizations based on dm sideband
            n_bkg_passed.setRange(0.4 * n_bkg_calib_guess_passed,
                                  2.0 * n_bkg_calib_guess_passed);
            n_bkg_passed.setVal(n_bkg_calib_guess_passed);
            n_bkg_failed.setRange(0.4 * n_bkg_calib_guess_failed,
                                  2.0 * n_bkg_calib_guess_failed);
            n_bkg_failed.setVal(n_bkg_calib_guess_failed);
            // Calib signal
            const double n_sig_guess =
                (n_calib_passed + n_calib_failed) -
                (n_bkg_calib_guess_passed + n_bkg_calib_guess_failed);
            n_sig.setRange(0.5 * n_sig_guess, 2.0 * n_sig_guess);
            n_sig.setVal(n_sig_guess);
            const double n_sig_guess_passed = f_guess_passed * n_calib_passed;
            const double eff_guess          = n_sig_guess_passed / n_sig_guess;
            eff.setVal(eff_guess);
            // MC
            // These normalizations are known so they are set as constant at
            // construction
            n_mc_passed.setVal(n_mc_guess_passed_dif +
                               n_mc_guess_passed_nondif);
            n_mc_failed.setVal(n_mc_guess_failed_dif +
                               n_mc_guess_failed_nondif);

            // Estimate exponential coefficients
            constexpr double dx = 1.8975 - 1.8275;

            const double y1_passed = dataset_calib_passed->sumEntries(
                "(d0_m_var > 1.825) && (d0_m_var < 1.830)");
            const double y2_passed = dataset_calib_passed->sumEntries(
                "(d0_m_var > 1.895) && (d0_m_var < 1.900)");
            const double k_passed_guess = log(y2_passed / y1_passed) / dx;
            if (k_passed_guess < 0.) {
              k_part_reco_passed.setVal(k_passed_guess);
            } else {
              k_part_reco_passed.setVal(-1e-4);
            }

            const double y1_failed = dataset_calib_failed->sumEntries(
                "(d0_m_var > 1.825) && (d0_m_var < 1.830)");
            const double y2_failed = dataset_calib_failed->sumEntries(
                "(d0_m_var > 1.895) && (d0_m_var < 1.900)");
            const double k_failed_guess = log(y2_failed / y1_failed) / dx;
            if (k_failed_guess < 0.) {
              k_comb_all.setVal(k_failed_guess);
              k_part_reco_failed.setVal(k_failed_guess);
            } else {
              k_comb_all.setVal(-1e-4);
              k_part_reco_failed.setVal(-1e-4);
            }

            // Fine tuning
            if ((probe == "pi") && (p_idx <= 1)) {
              // Use reduced fit range excluding higher D0 mass region to avoid
              // mass threshold in certain kinematic bins
              d0_m_var.setRange("fit", 1.825, 1.900);
              dm_var.setRange("fit", 0.141, 0.153);
            } else {
              // Use full fit range
              d0_m_var.setRange("fit", 1.825, 1.910);
              dm_var.setRange("fit", 0.141, 0.153);
            }

            //

            RooCmdArg fit_strat  = Strategy(2);
            RooCmdArg numcpu     = NumCPU(6);
            RooCmdArg offset     = Offset();
            RooCmdArg extended   = Extended();
            RooCmdArg printlevel = PrintLevel(0);
            RooCmdArg minos      = Minos(argset_minos);
            RooCmdArg timer      = Timer();
            RooCmdArg save       = Save();
            RooCmdArg range      = Range("fit");
            RooCmdArg minimize   = Minimizer("Minuit", "minimize");

            RooLinkedList fit_args;
            fit_args.Add(&fit_strat);
            fit_args.Add(&numcpu);
            fit_args.Add(&offset);
            fit_args.Add(&extended);
            fit_args.Add(&printlevel);
            fit_args.Add(&timer);
            fit_args.Add(&save);
            fit_args.Add(&range);
            fit_args.Add(&minimize);
            if (use_minos) fit_args.Add(&minos);

            cout << "\nINFO Starting calib only fit " << suffix << endl;

            set_parameters_calib(probe);

            if (!vmu && !fake_mu) {
              if ((probe == "pi") && (p_idx < 1)) {
                m_shift_d0_passed.setConstant();
                m_shift_d0_passed.setVal(0.);
              } else {
                m_shift_d0_passed.setConstant(false);
              }
            }

            // if (probe == "k") {
            //   //   f2_bkg_failed.setConstant();
            //   f2_bkg_passed.setConstant();
            //   //   k_part_reco_failed.setConstant();
            //   k_part_reco_passed.setConstant();
            //   //   f2_bkg_failed.setVal(0);
            //   f2_bkg_passed.setVal(0);
            // } else {
            //   //   f2_bkg_failed.setConstant(false);
            //   f2_bkg_passed.setConstant(false);
            //   //   k_part_reco_failed.setConstant(false);
            //   k_part_reco_passed.setConstant(false);
            // }

            // Create RooFitResult pointer out of if block to track lastest data
            // fit (either the standalone calib fit, or the full simultaneous
            // fit if it is performed)
            std::unique_ptr<RooFitResult> fit_result(nullptr);

            if (!dry_run && !mc_only) {
              auto start_calib = high_resolution_clock::now();

              fit_result.reset(model_calib.fitTo(dataset, fit_args));

              cout << "\nINFO Fit status: " << fit_result->status() << endl;
              cout << "INFO Covariance matrix status: " << fit_result->covQual()
                   << "\n"
                   << endl;

              int calib_fit_reattempts = 0;
              while (!fit_ok(fit_result.get()) &&
                     (calib_fit_reattempts < max_fix_reattempts)) {
                calib_fit_reattempts++;
                cout << "INFO Retrying calib only fit " << suffix << " (retry #"
                     << calib_fit_reattempts << ")" << endl;

                model_calib.fitTo(dataset, Extended(), NumCPU(6), PrintLevel(0),
                                  Offset(), Range("fit"),
                                  Minimizer("Minuit", "simplex"), Hesse(false),
                                  Timer());
                fit_result.reset(model_calib.fitTo(dataset, fit_args));

                cout << "\nINFO Fit status: " << fit_result->status() << endl;
                cout << "INFO Covariance matrix status: "
                     << fit_result->covQual() << "\n"
                     << endl;
              }

              auto stop_calib = high_resolution_clock::now();

              const int duration_calib =
                  duration_cast<seconds>(stop_calib - start_calib).count();

              cout << "INFO Calib fits took " << format_time(duration_calib)
                   << endl;
              cout << "INFO Needed " << calib_fit_reattempts << " retries with "
                   << format_time(duration_calib / (calib_fit_reattempts + 1))
                   << " on average\n"
                   << endl;

              // Get reduced chi2
              RooChi2Var chi2_calib(
                  "chi2_calib", "chi2_calib", model_calib, datahist_calib,
                  Extended(true), Range("fit"), DataError(RooAbsData::Poisson));
              // FIXME ndof depends on hard-coded number of bins in roodatahist
              const int ndof_calib =
                  2 * 250 - fit_result->floatParsFinal().size();
              const double chi2ndof_calib = chi2_calib.getVal() / ndof_calib;

              cout << "INFO Finished calib only fit " << suffix
                   << " with chi2ndof = " << chi2ndof_calib << "\n"
                   << endl;

              calib_retries.Fill(calib_fit_reattempts);

              if (fit_ok(fit_result.get())) {
                m_shift_d0_pull =
                    m_shift_d0_delta.getVal() /
                    m_shift_d0_delta.getPropagatedError(*fit_result);

                ds_params_calib.addFast(params_calib);
              }

              fit_status_calib.SetBinContent(kin_bin, fit_result->status());
              fit_cov_qual_calib.SetBinContent(kin_bin, fit_result->covQual());
            }

            // Plot fit results
            unique_ptr<RooPlot> frame_d0_calib_only_passed(
                d0_m_var.frame(Title("D0 M Calib Passed " + tag)));
            unique_ptr<RooPlot> frame_d0_calib_only_failed(
                d0_m_var.frame(Title("D0 M Calib Failed " + tag)));
            unique_ptr<RooPlot> frame_dm_calib_only_passed(
                dm_var.frame(Title("dm Calib Passed " + tag)));
            unique_ptr<RooPlot> frame_dm_calib_only_failed(
                dm_var.frame(Title("dm Calib Failed " + tag)));

            // TODO
            // https://root-forum.cern.ch/t/simultaneous-fit-normalization-issue/33965
            dataset.plotOn(frame_d0_calib_only_passed.get(),
                           Binning(bins_d0_m_passed),
                           Cut("sample==sample::calib_passed"));
            dataset.plotOn(frame_d0_calib_only_failed.get(), Binning(bins_d0_m),
                           Cut("sample==sample::calib_failed"));
            dataset.plotOn(
                frame_dm_calib_only_passed.get(), Binning(bins_dm_passed),
                Cut("sample==sample::calib_passed"), CutRange("fit"));
            dataset.plotOn(frame_dm_calib_only_failed.get(), Binning(bins_dm),
                           Cut("sample==sample::calib_failed"),
                           CutRange("fit"));

            model_calib.plotOn(frame_d0_calib_only_passed.get(),
                               Slice(sample, "calib_passed"),
                               ProjWData(sample, dataset), LineWidth(1));
            model_calib.plotOn(frame_d0_calib_only_passed.get(),
                               Slice(sample, "calib_passed"),
                               ProjWData(sample, dataset), LineWidth(1),
                               Components(d0_comb_passed), LineColor(kMagenta));
            model_calib.plotOn(frame_d0_calib_only_passed.get(),
                               Slice(sample, "calib_passed"),
                               ProjWData(sample, dataset), LineWidth(1),
                               Components(dm_comb_passed), LineColor(kGreen));
            model_calib.plotOn(frame_d0_calib_only_passed.get(),
                               Slice(sample, "calib_passed"),
                               ProjWData(sample, dataset), LineWidth(1),
                               Components(dmd0_comb_passed), LineColor(kRed));

            model_calib.plotOn(frame_d0_calib_only_failed.get(),
                               Slice(sample, "calib_failed"),
                               ProjWData(sample, dataset), LineWidth(1));
            model_calib.plotOn(frame_d0_calib_only_failed.get(),
                               Slice(sample, "calib_failed"),
                               ProjWData(sample, dataset), LineWidth(1),
                               Components(d0_comb_failed), LineColor(kMagenta));
            model_calib.plotOn(frame_d0_calib_only_failed.get(),
                               Slice(sample, "calib_failed"),
                               ProjWData(sample, dataset), LineWidth(1),
                               Components(dm_comb_failed), LineColor(kGreen));
            model_calib.plotOn(frame_d0_calib_only_failed.get(),
                               Slice(sample, "calib_failed"),
                               ProjWData(sample, dataset), LineWidth(1),
                               Components(dmd0_comb_failed), LineColor(kRed));

            model_calib.plotOn(frame_dm_calib_only_passed.get(),
                               Slice(sample, "calib_passed"),
                               ProjWData(sample, dataset),
                               ProjectionRange("fit"), LineWidth(1));
            model_calib.plotOn(
                frame_dm_calib_only_passed.get(), Slice(sample, "calib_passed"),
                ProjWData(sample, dataset), ProjectionRange("fit"),
                LineWidth(1), Components(d0_comb_passed), LineColor(kMagenta));
            model_calib.plotOn(
                frame_dm_calib_only_passed.get(), Slice(sample, "calib_passed"),
                ProjWData(sample, dataset), ProjectionRange("fit"),
                LineWidth(1), Components(dm_comb_passed), LineColor(kGreen));
            model_calib.plotOn(
                frame_dm_calib_only_passed.get(), Slice(sample, "calib_passed"),
                ProjWData(sample, dataset), ProjectionRange("fit"),
                LineWidth(1), Components(dmd0_comb_passed), LineColor(kRed));

            model_calib.plotOn(frame_dm_calib_only_failed.get(),
                               Slice(sample, "calib_failed"),
                               ProjWData(sample, dataset),
                               ProjectionRange("fit"), LineWidth(1));
            model_calib.plotOn(
                frame_dm_calib_only_failed.get(), Slice(sample, "calib_failed"),
                ProjWData(sample, dataset), ProjectionRange("fit"),
                LineWidth(1), Components(d0_comb_failed), LineColor(kMagenta));
            model_calib.plotOn(
                frame_dm_calib_only_failed.get(), Slice(sample, "calib_failed"),
                ProjWData(sample, dataset), ProjectionRange("fit"),
                LineWidth(1), Components(dm_comb_failed), LineColor(kGreen));
            model_calib.plotOn(
                frame_dm_calib_only_failed.get(), Slice(sample, "calib_failed"),
                ProjWData(sample, dataset), ProjectionRange("fit"),
                LineWidth(1), Components(dmd0_comb_failed), LineColor(kRed));

            c_four.cd(1);
            frame_d0_calib_only_passed->Draw();
            c_four.cd(2);
            frame_dm_calib_only_passed->Draw();
            c_four.cd(3);
            frame_d0_calib_only_failed->Draw();
            c_four.cd(4);
            frame_dm_calib_only_failed->Draw();

            c_four.SaveAs(opath + "/figs/fits/fit_calib_" + suffix + ".pdf");

            if (run_simultaneous_mc_calib) {
              cout << "\nINFO Starting full simultaneous fit " << suffix
                   << endl;

              mean_d0_failed_dif.setConstant(false);
              mean_d0_passed_dif.setConstant(false);
              mean_d0_failed_nondif.setConstant(false);
              mean_d0_passed_nondif.setConstant(false);
              width_d0_dif.setConstant(false);
              width_d0_nondif.setConstant(false);
              width_cb_d0_dif.setConstant(false);
              width_cb_d0_nondif.setConstant(false);
              alpha_L_d0_failed_dif.setConstant(false);
              alpha_L_d0_passed_dif.setConstant(false);
              alpha_L_d0_failed_nondif.setConstant(false);
              alpha_L_d0_passed_nondif.setConstant(false);
              n_L_d0_dif.setConstant(false);
              n_L_d0_nondif.setConstant(false);
              alpha_R_passed_d0_dif.setConstant(false);
              alpha_R_failed_d0_dif.setConstant(false);
              alpha_R_passed_d0_nondif.setConstant(false);
              alpha_R_failed_d0_nondif.setConstant(false);
              n_R_d0_dif.setConstant(false);
              n_R_d0_nondif.setConstant(false);
              f_d0_passed_dif.setConstant(false);
              f_d0_failed_dif.setConstant(false);
              f_d0_passed_nondif.setConstant(false);
              f_d0_failed_nondif.setConstant(false);
              mean_dm_failed_dif.setConstant(false);
              mean_dm_passed_dif.setConstant(false);
              mean_dm_failed_nondif.setConstant(false);
              mean_dm_passed_nondif.setConstant(false);
              width_dm_dif.setConstant(false);
              width_dm_nondif.setConstant(false);
              width_cb_dm_dif.setConstant(false);
              width_cb_dm_nondif.setConstant(false);
              alpha_L_dm_dif.setConstant(false);
              alpha_L_dm_nondif.setConstant(false);
              n_L_dm_dif.setConstant(false);
              n_L_dm_nondif.setConstant(false);
              alpha_R_dm_dif.setConstant(false);
              alpha_R_dm_nondif.setConstant(false);
              n_R_dm_dif.setConstant(false);
              n_R_dm_nondif.setConstant(false);
              f_dm_passed_dif.setConstant(false);
              f_dm_failed_dif.setConstant(false);
              f_dm_passed_nondif.setConstant(false);
              f_dm_failed_nondif.setConstant(false);

              if (!dry_run && !mc_only) {
                auto start_full = high_resolution_clock::now();

                fit_result.reset(model.fitTo(dataset, fit_args));

                cout << "\nINFO Fit status: " << fit_result.get()->status()
                     << endl;
                cout << "INFO Covariance matrix status: "
                     << fit_result.get()->covQual() << "\n"
                     << endl;

                int full_fit_reattempts = 0;
                while (!fit_ok(fit_result.get()) &&
                       (full_fit_reattempts < max_fix_reattempts)) {
                  full_fit_reattempts++;
                  cout << "INFO Retrying full fit " << suffix << " (retry #"
                       << full_fit_reattempts << ")" << endl;

                  model.fitTo(dataset, Extended(), NumCPU(8), PrintLevel(0),
                              Offset(), Range("fit"),
                              Minimizer("Minuit", "scan"), Timer(),
                              Hesse(false));
                  fit_result.reset(model_calib.fitTo(dataset, fit_args));

                  cout << "\nINFO Fit status: " << fit_result.get()->status()
                       << endl;
                  cout << "INFO Covariance matrix status: "
                       << fit_result.get()->covQual() << "\n"
                       << endl;
                }

                auto stop_full = high_resolution_clock::now();

                const int duration_full =
                    duration_cast<seconds>(stop_full - start_full).count();

                cout << "INFO Full fits took " << format_time(duration_full)
                     << endl;
                cout << "INFO Needed " << full_fit_reattempts
                     << " retries with "
                     << format_time(duration_full / (full_fit_reattempts + 1))
                     << " on average\n"
                     << endl;

                fit_status.SetBinContent(kin_bin, fit_result.get()->status());
                fit_cov_qual.SetBinContent(kin_bin,
                                           fit_result.get()->covQual());
              }
            }

            // Save fitted PID efficiency
            histo_pid_raw.SetBinContent(kin_bin, eff.getVal());
            histo_pid.SetBinContent(kin_bin, eff_corrected.getVal());
            // Save fitted DiF fractions
            histo_f_dif.SetBinContent(kin_bin, f_dif_calib_passed.getVal());
            // Save errors if fit was performed
            if (fit_result.get()) {
              histo_pid_raw.SetBinError(kin_bin, eff.getError());
              histo_pid.SetBinError(
                  kin_bin, eff_corrected.getPropagatedError(*fit_result));
              histo_f_dif.SetBinError(
                  kin_bin, f_dif_calib_passed.getPropagatedError(*fit_result));
            }

            // Store final mass-window efficiencies in histogram
            effs_passed.SetBinContent(kin_bin, eff_mw_passed.getVal());
            // effs_passed_unc_hi.SetBinContent(kin_bin, eff_passed_unc_hi);
            // effs_passed_unc_lo.SetBinContent(kin_bin, eff_passed_unc_lo);
            effs_failed.SetBinContent(kin_bin, eff_mw_failed.getVal());
            // effs_failed_unc_hi.SetBinContent(kin_bin, eff_failed_unc_hi);
            // effs_failed_unc_lo.SetBinContent(kin_bin, eff_failed_unc_lo);

            // // Print mass-window efficiencies
            // cout << "\nINFO Mass-window efficiency for passed sample: "
            //      << in_mw_passed << " / " << total_passed << " = " <<
            //      eff_passed
            //      << endl;
            // cout << "\nINFO Mass-window efficiency for failed sample: "
            //      << in_mw_failed << " / " << total_failed << " = " <<
            //      eff_failed
            //      << endl;

            // Check normalization estimates
            cout << "\nINFO Fitted vs estimated:" << endl;
            cout << " - Fitted n_sig = " << n_sig.getVal() << " vs estimated "
                 << n_sig_guess << " (" << n_sig.getVal() / n_sig_guess << ")"
                 << endl;
            cout << " - Fitted n_bkg_failed = " << n_bkg_failed.getVal()
                 << " vs estimated " << n_bkg_calib_guess_failed << " ("
                 << n_bkg_failed.getVal() / n_bkg_calib_guess_failed << ")"
                 << endl;
            cout << " - Fitted n_bkg_passed = " << n_bkg_passed.getVal()
                 << " vs estimated " << n_bkg_calib_guess_passed << " ("
                 << n_bkg_passed.getVal() / n_bkg_calib_guess_passed << ")\n"
                 << endl;

            cout << " - Fitted k_part_reco_passed = "
                 << k_part_reco_passed.getVal() << " vs estimated "
                 << k_passed_guess << " ("
                 << k_part_reco_passed.getVal() / k_passed_guess << ")" << endl;
            cout << " - Fitted k_part_reco_failed = "
                 << k_part_reco_failed.getVal() << " vs estimated "
                 << k_failed_guess << " ("
                 << k_part_reco_failed.getVal() / k_failed_guess << ")" << endl;
            cout << " - Fitted k_comb_all = " << k_comb_all.getVal()
                 << " vs estimated " << k_failed_guess << " ("
                 << k_comb_all.getVal() / k_failed_guess << ")\n"
                 << endl;

            // Plot fit results
            unique_ptr<RooPlot> frame_d0_mc_passed_nondif(
                d0_m_var.frame(Title("D0 M MC Passed " + tag)));
            unique_ptr<RooPlot> frame_d0_calib_passed(
                d0_m_var.frame(Title("D0 M Calib Passed " + tag)));
            unique_ptr<RooPlot> frame_d0_mc_failed_nondif(
                d0_m_var.frame(Title("D0 M MC Failed " + tag)));
            unique_ptr<RooPlot> frame_d0_calib_failed(
                d0_m_var.frame(Title("D0 M Calib Failed " + tag)));
            unique_ptr<RooPlot> frame_dm_mc_passed_nondif(
                dm_var.frame(Title("dm MC Passed " + tag)));
            unique_ptr<RooPlot> frame_dm_calib_passed(
                dm_var.frame(Title("dm Calib Passed " + tag)));
            unique_ptr<RooPlot> frame_dm_mc_failed_nondif(
                dm_var.frame(Title("dm MC Failed " + tag)));
            unique_ptr<RooPlot> frame_dm_calib_failed(
                dm_var.frame(Title("dm Calib Failed " + tag)));

            dataset.plotOn(frame_d0_mc_passed_nondif.get(),
                           Binning(bins_d0_m_passed),
                           Cut("sample==sample::mc_passed_nondif"));
            dataset.plotOn(frame_d0_mc_failed_nondif.get(), Binning(bins_d0_m),
                           Cut("sample==sample::mc_failed_nondif"));
            dataset.plotOn(frame_d0_calib_passed.get(),
                           Binning(bins_d0_m_passed),
                           Cut("sample==sample::calib_passed"));
            dataset.plotOn(frame_d0_calib_failed.get(), Binning(bins_d0_m),
                           Cut("sample==sample::calib_failed"));
            dataset.plotOn(frame_dm_mc_passed_nondif.get(),
                           Binning(bins_dm_passed),
                           Cut("sample==sample::mc_passed_nondif"));
            dataset.plotOn(frame_dm_mc_failed_nondif.get(), Binning(bins_dm),
                           Cut("sample==sample::mc_failed_nondif"));
            dataset.plotOn(frame_dm_calib_passed.get(), Binning(bins_dm_passed),
                           Cut("sample==sample::calib_passed"),
                           CutRange("fit"));
            dataset.plotOn(frame_dm_calib_failed.get(), Binning(bins_dm),
                           Cut("sample==sample::calib_failed"),
                           CutRange("fit"));

            d0_model_mc_nondif.plotOn(frame_d0_mc_passed_nondif.get(),
                                      Slice(sample, "mc_passed_nondif"),
                                      ProjWData(sample, dataset), LineWidth(1));
            d0_model_mc_nondif.plotOn(frame_d0_mc_failed_nondif.get(),
                                      Slice(sample, "mc_failed_nondif"),
                                      ProjWData(sample, dataset), LineWidth(1));
            model_calib.plotOn(frame_d0_calib_passed.get(),
                               Slice(sample, "calib_passed"),
                               ProjWData(sample, dataset), LineWidth(1));
            model_calib.plotOn(frame_d0_calib_failed.get(),
                               Slice(sample, "calib_failed"),
                               ProjWData(sample, dataset), LineWidth(1));
            dm_model_mc_nondif.plotOn(frame_dm_mc_passed_nondif.get(),
                                      Slice(sample, "mc_passed_nondif"),
                                      ProjWData(sample, dataset), LineWidth(1));
            dm_model_mc_nondif.plotOn(frame_dm_mc_failed_nondif.get(),
                                      Slice(sample, "mc_failed_nondif"),
                                      ProjWData(sample, dataset), LineWidth(1));
            model_calib.plotOn(frame_dm_calib_passed.get(),
                               Slice(sample, "calib_passed"),
                               ProjWData(sample, dataset),
                               ProjectionRange("fit"), LineWidth(1));
            model_calib.plotOn(frame_dm_calib_failed.get(),
                               Slice(sample, "calib_failed"),
                               ProjWData(sample, dataset),
                               ProjectionRange("fit"), LineWidth(1));

            c_mult.cd(1);
            frame_d0_mc_passed_nondif->Draw();
            c_mult.cd(2);
            frame_dm_mc_passed_nondif->Draw();
            c_mult.cd(3);
            frame_d0_calib_passed->Draw();
            c_mult.cd(4);
            frame_dm_calib_passed->Draw();
            c_mult.cd(5);
            frame_d0_mc_failed_nondif->Draw();
            c_mult.cd(6);
            frame_dm_mc_failed_nondif->Draw();
            c_mult.cd(7);
            frame_d0_calib_failed->Draw();
            c_mult.cd(8);
            frame_dm_calib_failed->Draw();

            c_mult.SaveAs(opath + "/figs/fits/fit_" + suffix + ".pdf");
          }
        }
      }

      // Define output file por PID efficiency
      // File name matches the output of PIDCalib
      const TString fname_suffix = fake_mu ? "_denom" : "_nom";
      const TString opath_full_probe =
          opath + "/" + probe + "TrueToMuTag" + fname_suffix + ".root";
      cout << "INFO Creating output file: " << opath_full_probe << endl;
      TFile ofile_probe(opath_full_probe, "RECREATE");

      // Save PID efficiency
      ofile_probe.cd();
      histo_pid_raw.Write();
      histo_pid.Write();
      ofile_probe.Close();

      // Define output file for fitted dif fractions
      const TString opath_f_dif = opath + "/dif_fractions_" + probe + ".root";
      cout << "INFO Creating output file: " << opath_f_dif << endl;
      TFile ofile_f_dif(opath_f_dif, "RECREATE");

      // Save fitted DiF fractions
      ofile_f_dif.cd();
      histo_f_dif.Write();
      ofile_f_dif.Close();

      // Now save other distributions in common output file
      ofile.cd();
      fit_status_d0_mc_dif.Write();
      fit_cov_qual_d0_mc_dif.Write();
      fit_status_d0_mc_nondif.Write();
      fit_cov_qual_d0_mc_nondif.Write();
      fit_status_dm_mc_dif.Write();
      fit_cov_qual_dm_mc_dif.Write();
      fit_status_dm_mc_nondif.Write();
      fit_cov_qual_dm_mc_nondif.Write();
      fit_status_calib.Write();
      fit_cov_qual_calib.Write();
      if (run_simultaneous_mc_calib) {
        fit_status.Write();
        fit_cov_qual.Write();
      }

      // Save fit status
      suffix.Form("%s_%s", year.c_str(), probe.c_str());
      c_single.cd();
      fit_status_d0_mc_dif.Draw("BOX2Z");
      c_single.SaveAs(opath + "/figs/fit_status_d0_mc_dif_" + suffix + ".pdf");
      fit_cov_qual_d0_mc_dif.Draw("BOX2Z");
      c_single.SaveAs(opath + "/figs/fit_cov_qual_d0_mc_dif_" + suffix +
                      ".pdf");
      fit_status_d0_mc_nondif.Draw("BOX2Z");
      c_single.SaveAs(opath + "/figs/fit_status_d0_mc_nondif_" + suffix +
                      ".pdf");
      fit_cov_qual_d0_mc_nondif.Draw("BOX2Z");
      c_single.SaveAs(opath + "/figs/fit_cov_qual_d0_mc_nondif_" + suffix +
                      ".pdf");
      fit_status_dm_mc_dif.Draw("BOX2Z");
      c_single.SaveAs(opath + "/figs/fit_status_dm_mc_dif_" + suffix + ".pdf");
      fit_cov_qual_dm_mc_dif.Draw("BOX2Z");
      c_single.SaveAs(opath + "/figs/fit_cov_qual_dm_mc_dif_" + suffix +
                      ".pdf");
      fit_status_dm_mc_nondif.Draw("BOX2Z");
      c_single.SaveAs(opath + "/figs/fit_status_dm_mc_nondif_" + suffix +
                      ".pdf");
      fit_cov_qual_dm_mc_nondif.Draw("BOX2Z");
      c_single.SaveAs(opath + "/figs/fit_cov_qual_dm_mc_nondif_" + suffix +
                      ".pdf");
      fit_status_calib.Draw("BOX2Z");
      c_single.SaveAs(opath + "/figs/fit_status_calib_" + suffix + ".pdf");
      fit_cov_qual_calib.Draw("BOX2Z");
      c_single.SaveAs(opath + "/figs/fit_cov_qual_calib_" + suffix + ".pdf");
      if (run_simultaneous_mc_calib) {
        fit_status.Draw("BOX2Z");
        c_single.SaveAs(opath + "/figs/fit_status_" + suffix + ".pdf");
        fit_cov_qual.Draw("BOX2Z");
        c_single.SaveAs(opath + "/figs/fit_cov_qual_" + suffix + ".pdf");
      }

      // Save number of fit retries
      d0_retries_dif.Draw();
      c_single.SaveAs(opath + "/figs/fit_retries_d0_mc_" + suffix + "_dif.pdf");
      d0_retries_nondif.Draw();
      c_single.SaveAs(opath + "/figs/fit_retries_d0_mc_" + suffix +
                      "_nondif.pdf");
      dm_retries_dif.Draw();
      c_single.SaveAs(opath + "/figs/fit_retries_dm_mc_" + suffix + "_dif.pdf");
      dm_retries_nondif.Draw();
      c_single.SaveAs(opath + "/figs/fit_retries_dm_mc_" + suffix +
                      "_nondif.pdf");
      calib_retries.Draw();
      c_single.SaveAs(opath + "/figs/fit_retries_calib_" + suffix + ".pdf");

      // Plot distribution of fitted variables
      if (!dry_run) {
        plot_dataset(ds_params_d0_dif, opath + "/figs/params/" + suffix + "_");
        plot_dataset(ds_params_d0_nondif,
                     opath + "/figs/params/" + suffix + "_");
        plot_dataset(ds_params_dm_dif, opath + "/figs/params/" + suffix + "_");
        plot_dataset(ds_params_dm_nondif,
                     opath + "/figs/params/" + suffix + "_");
        if (!mc_only) {
          plot_dataset(ds_params_calib, opath + "/figs/params/" + suffix + "_");
        } else {
          cout << "INFO MC fits only: will not plot distributions of calib fit "
                  "parameters"
               << endl;
        }
      } else {
        cout << "INFO Dry run: will not plot distributions of fit parameters"
             << endl;
      }

      // Save efficiencies
      cout << "INFO Saving efficiencies " << endl;
      effs_passed.Write();
      effs_passed_unc_hi.Write();
      effs_passed_unc_lo.Write();
      effs_failed.Write();
      effs_failed_unc_hi.Write();
      effs_failed_unc_lo.Write();

      // Delete calib datasets
      cout << "INFO Deleting MC datasets " << endl;
      for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
        for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
          for (int ntrks_idx = 0; ntrks_idx < N_BINS_NTRACKS; ntrks_idx++) {
            delete datasets_calib_passed[ntrks_idx][eta_idx][p_idx];
            delete datasets_calib_failed[ntrks_idx][eta_idx][p_idx];
          }
        }
      }
    }

    // Delete MC datasets
    cout << "INFO Deleting MC datasets " << endl;
    for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
      for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
        delete datasets_mc_passed_dif[eta_idx][p_idx];
        delete datasets_mc_failed_dif[eta_idx][p_idx];
        for (int ntrks_idx = 0; ntrks_idx < N_BINS_NTRACKS; ntrks_idx++) {
          delete datasets_mc_passed_nondif[ntrks_idx][eta_idx][p_idx];
          delete datasets_mc_failed_nondif[ntrks_idx][eta_idx][p_idx];
        }
      }
    }
  }

  ofile.Close();

  auto stop = high_resolution_clock::now();

  const int duration = duration_cast<seconds>(stop - start).count();
  cout << "INFO Finished processing in " << format_time(duration) << endl;
}
