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
#include "RooArgSet.h"
#include "RooBinning.h"
#include "RooCategory.h"
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
constexpr double DM_min    = 141.0;
constexpr double DM_max    = 153.0;
constexpr int    DM_nbins  = 60;
constexpr int    D0_ID     = 421;
constexpr double D0_M_min  = 1825.0;
constexpr double D0_M_max  = 1910.0;
constexpr double D0_M      = 1864.84;
constexpr int    D0_nbins  = 40;
constexpr int    PI_ID     = 211;
constexpr double PI_M      = 139.57039;
constexpr int    K_ID      = 321;
constexpr double ONE_SIGMA = 0.682689492137;

constexpr double BINS_P[]       = {3.e+3,  6.e+3,  10.e+3, 15.6e+3,
                                   27.e+3, 60.e+3, 100.e+3};
constexpr double BINS_ETA[]     = {1.7, 3.6, 5.0};
constexpr double BINS_NTRACKS[] = {0, 200, 600};

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
  auto    vars        = ds.get(0);
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
    ("h,help", "print help")
    ("d,debug", "enable debug mode",
     cxxopts::value<bool>()->default_value("false"))
    ("p,particles", "specify probed particle",
     cxxopts::value<vector<string>>()->default_value("pi,k"))
    ("c,config", "specify input YAML config file",
     cxxopts::value<string>())
    ("years", "specify input YAML config file",
     cxxopts::value<vector<string>>()->default_value("2016,2017,2018"))
    ("o,output", "specify output folder",
     cxxopts::value<string>()->default_value("gen/"))
    ("vmu", "flag misid validation PID cuts",
     cxxopts::value<bool>()->default_value("false"))
    ("fake_mu", "flag fake muon sample PID cuts",
     cxxopts::value<bool>()->default_value("false"))
    ;
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  // Read arguments
  const auto    ymlFile   = parsedArgs["config"].as<string>();
  const auto    years     = parsedArgs["years"].as<vector<string>>();
  const auto    particles = parsedArgs["particles"].as<vector<string>>();
  const TString opath     = parsedArgs["output"].as<string>();
  const auto    vmu       = parsedArgs["vmu"].as<bool>();
  const auto    fake_mu   = parsedArgs["fake_mu"].as<bool>();

  const bool run_simultaneous_mc_calib = false;

  if (vmu) cout << "INFO Using VMU pid cuts" << endl;

  if (fake_mu) cout << "INFO Using FAKE_MU pid cuts" << endl;

  if (vmu && fake_mu) {
    cout << "WARNING Both VMU and FAKE_MU flags set. The VMU flag is redundant "
            "as the FAKE MU flag implies no muBDT cut"
         << endl;
  }

  if (run_simultaneous_mc_calib) {
    cout << "WARNING Configured to run full simultaneous fit, which can take "
            "2 days to finish!"
         << endl;
  }

  constexpr int max_fix_reattempts = 10;

  // Parse YAML config
  const auto ymlConfig = YAML::LoadFile(ymlFile)["misid_corrections"];

  // Define histogram to easily determine kinematical bins
  TH3D histo_binnning("histo_binnning", ";p;#eta;nTracks", N_BINS_P, BINS_P,
                      N_BINS_ETA, BINS_ETA, N_BINS_NTRACKS, BINS_NTRACKS);

  // Define counters for mass window efficiencies and initialize all to 0
  int count_in_mw_passed[3][N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] = {{{{0}}}};
  int count_total_passed[3][N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] = {{{{0}}}};
  int count_in_mw_failed[3][N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] = {{{{0}}}};
  int count_total_failed[3][N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] = {{{{0}}}};

  // Define variable to access input ntuples
  int    dst_id, d0_id, ntracks;
  double dst_m, d0_m, d0_pt, probe_p, probe_pz, probe_pt, probe_eta,
      probe_dllmu, probe_dlle, probe_ghostprob, probe_mu_unbiased,
      ntracks_calib, spi_p, spi_pt, tag_p, tag_pt;
  float probe_mu_ubdt, probe_track_chi2ndof;
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
  RooRealVar d0_m_var("d0_m_var", "d0_m_var", 1825, 1910, "MeV/c^{2}");
  RooRealVar dm_var("dm_var", "dm_var", 141, 153, "MeV/c^{2}");
  RooBinning bins_d0_m(85, 1825, 1910, "bins_d0_m");
  RooBinning bins_dm(80, 141, 153, "bins_dm");
  RooBinning bins_d0_m_passed(50, 1825, 1910, "bins_mc_d0_passed");
  RooBinning bins_dm_passed(50, 141, 153, "bins_dm_passed");

  // Define categories
  RooCategory sample("sample", "sample");
  sample.defineType("mc_passed");
  sample.defineType("mc_failed");
  sample.defineType("calib_passed");
  sample.defineType("calib_failed");

  ///////////////////////////////
  // Individual fit components //
  ///////////////////////////////

  // D0 signal: CrystalBall + Gaussian

  // Core
  RooRealVar mean_d0_failed("mean_d0_failed", "mean_d0_failed", 1855, 1870,
                            "MeV");
  RooRealVar mean_d0_passed("mean_d0_passed", "mean_d0_passed", 1855, 1870,
                            "MeV");
  RooRealVar width_d0("width_d0", "width_d0", 1, 15, "MeV");
  RooRealVar width_cb_d0("width_cb_d0", "width_cb_d0", 1, 20, "MeV");

  RooRealVar m_shift_d0_failed("m_shift_d0_failed", "m_shift_d0_failed", -4, 4,
                               "MeV");
  RooRealVar m_shift_d0_passed("m_shift_d0_passed", "m_shift_d0_passed", -4, 4,
                               "MeV");
  RooRealVar w_scale_d0_failed("w_scale_d0_failed", "w_scale_d0_failed", 0.5,
                               1.5, "");
  RooRealVar w_scale_d0_passed("w_scale_d0_passed", "w_scale_d0_passed", 0.5,
                               1.5, "");

  RooFormulaVar mean_shifted_d0_failed(
      "mean_shifted_d0_failed", "mean_shifted_d0_failed", "x[0] + x[1]",
      RooArgList(mean_d0_failed, m_shift_d0_failed));
  RooFormulaVar mean_shifted_d0_passed(
      "mean_shifted_d0_passed", "mean_shifted_d0_passed", "x[0] + x[1]",
      RooArgList(mean_d0_passed, m_shift_d0_passed));
  RooFormulaVar width_d0_scaled_failed("width_d0_scaled_failed",
                                       "width_d0_scaled_failed", "x[0] * x[1]",
                                       RooArgList(width_d0, w_scale_d0_failed));
  RooFormulaVar width_d0_scaled_passed("width_d0_scaled_passed",
                                       "width_d0_scaled_passed", "x[0] * x[1]",
                                       RooArgList(width_d0, w_scale_d0_passed));
  RooFormulaVar width_cb_d0_scaled_failed(
      "width_cb_d0_scaled_failed", "width_cb_d0_scaled_failed", "x[0] * x[1]",
      RooArgList(width_cb_d0, w_scale_d0_failed));
  RooFormulaVar width_cb_d0_scaled_passed(
      "width_cb_d0_scaled_passed", "width_cb_d0_scaled_passed", "x[0] * x[1]",
      RooArgList(width_cb_d0, w_scale_d0_passed));

  // Left tail
  RooRealVar alpha_L_d0_failed("alpha_L_d0_failed", "alpha_L_d0_failed", 1e-3,
                               3, "");
  RooRealVar alpha_L_d0_passed("alpha_L_d0_passed", "alpha_L_d0_passed", 1e-3,
                               3, "");
  RooRealVar n_L_d0_failed("n_L_d0_failed", "n_L_d0_failed", 0.01, 50, "");
  RooRealVar n_L_d0_passed("n_L_d0_passed", "n_L_d0_passed", 0.01, 50, "");

  // Right tail
  RooRealVar alpha_R_d0_failed("alpha_R_d0_failed", "alpha_R_d0_failed", 1e-3,
                               3, "");
  RooRealVar alpha_R_d0_passed("alpha_R_d0_passed", "alpha_R_d0_passed", 1e-3,
                               3, "");
  RooRealVar n_R_d0_failed("n_R_d0_failed", "n_R_d0_failed", 0.01, 50, "");
  RooRealVar n_R_d0_passed("n_R_d0_passed", "n_R_d0_passed", 0.01, 50, "");

  // Gaussian fraction
  RooRealVar f_d0_passed("f_d0_passed", "f_d0_passed", 0, 1);
  RooRealVar f_d0_failed("f_d0_failed", "f_d0_failed", 0, 1);

  // MC PDFs
  RooGaussian d0_gauss_passed_mc("d0_gauss_passed_mc", "d0_gauss_passed_mc",
                                 d0_m_var, mean_d0_passed, width_d0);
  RooGaussian d0_gauss_failed_mc("d0_gauss_failed_mc", "d0_gauss_failed_mc",
                                 d0_m_var, mean_d0_failed, width_d0);

  RooCrystalBall d0_cb_passed_mc("d0_cb_passed_mc", "d0_cb_passed_mc", d0_m_var,
                                 mean_d0_passed, width_cb_d0, alpha_L_d0_passed,
                                 n_L_d0_passed, alpha_R_d0_passed,
                                 n_R_d0_passed);
  RooCrystalBall d0_cb_failed_mc("d0_cb_failed_mc", "d0_cb_failed_mc", d0_m_var,
                                 mean_d0_failed, width_cb_d0, alpha_L_d0_failed,
                                 n_L_d0_failed, alpha_R_d0_failed,
                                 n_R_d0_failed);

  RooAddPdf d0_model_passed_mc("d0_model_passed_mc", "d0_model_passed_mc",
                               d0_gauss_passed_mc, d0_cb_passed_mc,
                               f_d0_passed);
  RooAddPdf d0_model_failed_mc("d0_model_failed_mc", "d0_model_failed_mc",
                               d0_gauss_failed_mc, d0_cb_failed_mc,
                               f_d0_failed);
  // Data PDFs with floating shift
  RooGaussian d0_gauss_passed_calib(
      "d0_gauss_passed_calib", "d0_gauss_passed_calib", d0_m_var,
      mean_shifted_d0_passed, width_d0_scaled_passed);
  RooGaussian d0_gauss_failed_calib(
      "d0_gauss_failed_calib", "d0_gauss_failed_calib", d0_m_var,
      mean_shifted_d0_failed, width_d0_scaled_failed);

  RooCrystalBall d0_cb_passed_calib(
      "d0_cb_passed_calib", "d0_cb_passed_calib", d0_m_var,
      mean_shifted_d0_passed, width_cb_d0_scaled_passed, alpha_L_d0_passed,
      n_L_d0_passed, alpha_R_d0_passed, n_R_d0_passed);
  RooCrystalBall d0_cb_failed_calib(
      "d0_cb_failed_calib", "d0_cb_failed_calib", d0_m_var,
      mean_shifted_d0_failed, width_cb_d0_scaled_failed, alpha_L_d0_failed,
      n_L_d0_failed, alpha_R_d0_failed, n_R_d0_failed);

  RooAddPdf d0_model_passed_calib(
      "d0_model_passed_calib", "d0_model_passed_calib", d0_gauss_passed_calib,
      d0_cb_passed_calib, f_d0_passed);
  RooAddPdf d0_model_failed_calib(
      "d0_model_failed_calib", "d0_model_failed_calib", d0_gauss_failed_calib,
      d0_cb_failed_calib, f_d0_failed);

  // dm signal: CrystalBall + Gaussian

  // Core
  RooRealVar mean_dm_failed("mean_dm_failed", "mean_dm_failed", 145, 146,
                            "MeV");
  RooRealVar mean_dm_passed("mean_dm_passed", "mean_dm_passed", 145, 146,
                            "MeV");
  RooRealVar width_dm("width_dm", "width_dm", 0.21, 2, "MeV");
  RooRealVar width_cb_dm("width_cb_dm", "width_cb_dm", 0.1, 3, "MeV");

  RooRealVar m_shift_dm_failed("m_shift_dm_failed", "m_shift_dm_failed", -0.4,
                               0.4, "MeV");
  RooRealVar m_shift_dm_passed("m_shift_dm_passed", "m_shift_dm_passed", -0.4,
                               0.4, "MeV");
  RooRealVar w_scale_dm_failed("w_scale_dm_failed", "w_scale_dm_failed", 0.5,
                               1.5, "");
  RooRealVar w_scale_dm_passed("w_scale_dm_passed", "w_scale_dm_passed", 0.5,
                               1.5, "");

  RooFormulaVar mean_shifted_dm_failed(
      "mean_shifted_dm_failed", "mean_shifted_dm_failed", "x[0] + x[1]",
      RooArgList(mean_dm_failed, m_shift_dm_failed));
  RooFormulaVar mean_shifted_dm_passed(
      "mean_shifted_dm_passed", "mean_shifted_dm_passed", "x[0] + x[1]",
      RooArgList(mean_dm_passed, m_shift_dm_passed));
  RooFormulaVar width_dm_scaled_failed("width_dm_scaled_failed",
                                       "width_dm_scaled_failed", "x[0] * x[1]",
                                       RooArgList(width_dm, w_scale_dm_failed));
  RooFormulaVar width_dm_scaled_passed("width_dm_scaled_passed",
                                       "width_dm_scaled_passed", "x[0] * x[1]",
                                       RooArgList(width_dm, w_scale_dm_passed));
  RooFormulaVar width_cb_dm_scaled_failed(
      "width_cb_dm_scaled_failed", "width_cb_dm_scaled_failed", "x[0] * x[1]",
      RooArgList(width_cb_dm, w_scale_dm_failed));
  RooFormulaVar width_cb_dm_scaled_passed(
      "width_cb_dm_scaled_passed", "width_cb_dm_scaled_passed", "x[0] * x[1]",
      RooArgList(width_cb_dm, w_scale_dm_passed));

  // Left tail
  RooRealVar alpha_L_dm_failed("alpha_L_dm_failed", "alpha_L_dm_failed", 1e-3,
                               3, "");
  RooRealVar alpha_L_dm_passed("alpha_L_dm_passed", "alpha_L_dm_passed", 1e-3,
                               3, "");
  RooRealVar n_L_dm_failed("n_L_dm_failed", "n_L_dm_failed", 0.01, 50, "");
  RooRealVar n_L_dm_passed("n_L_dm_passed", "n_L_dm_passed", 0.01, 50, "");

  // Right tail
  RooRealVar alpha_R_dm_failed("alpha_R_dm_failed", "alpha_R_dm_failed", 1e-3,
                               3, "");
  RooRealVar alpha_R_dm_passed("alpha_R_dm_passed", "alpha_R_dm_passed", 1e-3,
                               3, "");
  RooRealVar n_R_dm_failed("n_R_dm_failed", "n_R_dm_failed", 0.01, 50, "");
  RooRealVar n_R_dm_passed("n_R_dm_passed", "n_R_dm_passed", 0.01, 50, "");

  // Gaussian fraction
  RooRealVar f_dm_passed("f_dm_passed", "f_dm_passed", 0, 1);
  RooRealVar f_dm_failed("f_dm_failed", "f_dm_failed", 0, 1);

  // MC PDFs
  RooGaussian dm_gauss_passed_mc("dm_gauss_passed_mc", "dm_gauss_passed_mc",
                                 dm_var, mean_dm_passed, width_dm);
  RooGaussian dm_gauss_failed_mc("dm_gauss_failed_mc", "dm_gauss_failed_mc",
                                 dm_var, mean_dm_failed, width_dm);

  RooCrystalBall dm_cb_passed_mc(
      "dm_cb_passed_mc", "dm_cb_passed_mc", dm_var, mean_dm_passed, width_cb_dm,
      alpha_L_dm_passed, n_L_dm_passed, alpha_R_dm_passed, n_R_dm_passed);
  RooCrystalBall dm_cb_failed_mc(
      "dm_cb_failed_mc", "dm_cb_failed_mc", dm_var, mean_dm_failed, width_cb_dm,
      alpha_L_dm_failed, n_L_dm_failed, alpha_R_dm_failed, n_R_dm_failed);

  RooAddPdf dm_model_passed_mc("dm_model_passed_mc", "dm_model_passed_mc",
                               dm_gauss_passed_mc, dm_cb_passed_mc,
                               f_dm_passed);
  RooAddPdf dm_model_failed_mc("dm_model_failed_mc", "dm_model_failed_mc",
                               dm_gauss_failed_mc, dm_cb_failed_mc,
                               f_dm_failed);

  // Data PDFs with floating shift and scaled width
  RooGaussian dm_gauss_passed_calib(
      "dm_gauss_passed_calib", "dm_gauss_passed_calib", dm_var,
      mean_shifted_dm_passed, width_dm_scaled_passed);
  RooGaussian dm_gauss_failed_calib(
      "dm_gauss_failed_calib", "dm_gauss_failed_calib", dm_var,
      mean_shifted_dm_failed, width_dm_scaled_failed);

  RooCrystalBall dm_cb_passed_calib(
      "dm_cb_passed_calib", "dm_cb_passed_calib", dm_var,
      mean_shifted_dm_passed, width_cb_dm_scaled_passed, alpha_L_dm_passed,
      n_L_dm_passed, alpha_R_dm_passed, n_R_dm_passed);
  RooCrystalBall dm_cb_failed_calib(
      "dm_cb_failed_calib", "dm_cb_failed_calib", dm_var,
      mean_shifted_dm_failed, width_cb_dm_scaled_failed, alpha_L_dm_failed,
      n_L_dm_failed, alpha_R_dm_failed, n_R_dm_failed);

  RooAddPdf dm_model_passed_calib(
      "dm_model_passed_calib", "dm_model_passed_calib", dm_gauss_passed_calib,
      dm_cb_passed_calib, f_dm_passed);
  RooAddPdf dm_model_failed_calib(
      "dm_model_failed_calib", "dm_model_failed_calib", dm_gauss_failed_calib,
      dm_cb_failed_calib, f_dm_failed);

  // dm comb model: threshold function
  RooConstVar dm0("dm0", "dm0", PI_M);

  RooRealVar c_comb_spi_failed("c_comb_spi_failed", "c_comb_spi_failed", 0, 1,
                               "");
  RooRealVar c_comb_spi_passed("c_comb_spi_passed", "c_comb_spi_passed", 0, 1,
                               "");

  RooRealVar c_comb_all_failed("c_comb_all_failed", "c_comb_all_failed", 0, 1,
                               "");
  RooRealVar c_comb_all_passed("c_comb_all_passed", "c_comb_all_passed", 0, 1,
                               "");

  RooPowerLaw dm_comb_spi_failed("dm_comb_spi_failed", "dm_comb_spi_failed",
                                 dm_var, dm0, c_comb_spi_failed);
  RooPowerLaw dm_comb_spi_passed("dm_comb_spi_passed", "dm_comb_spi_passed",
                                 dm_var, dm0, c_comb_spi_passed);

  RooPowerLaw dm_comb_all_failed("dm_comb_all_failed", "dm_comb_all_failed",
                                 dm_var, dm0, c_comb_all_failed);
  RooPowerLaw dm_comb_all_passed("dm_comb_all_passed", "dm_comb_all_passed",
                                 dm_var, dm0, c_comb_all_passed);

  // D0 comb model: exponential distribution
  RooRealVar k_part_reco_failed("k_part_reco_failed", "k_part_reco_failed", -1,
                                0, "");
  RooRealVar k_part_reco_passed("k_part_reco_passed", "k_part_reco_passed", -1,
                                0, "");

  RooRealVar k_comb_all_failed("k_comb_all_failed", "k_comb_all_failed", -0.10,
                               0, "");
  RooRealVar k_comb_all_passed("k_comb_all_passed", "k_comb_all_passed", -0.10,
                               0, "");

  // Use shifted d0_m for exponentials to avoid hitting numerical precision
  // limits. Amounts to a simply normalizations change which is cancelled out by
  // Roofit.
  RooConstVar   d0_m_pdg("d0_m_pdg", "d0_m_pdg", D0_M);
  RooFormulaVar d0_m_var_shifted("d0_m_var_shifted", "d0_m_var_shifted",
                                 "x[0] - x[1]", RooArgList(d0_m_var, d0_m_pdg));

  RooExponential d0_part_reco_failed("d0_part_reco_failed",
                                     "d0_part_reco_failed", d0_m_var_shifted,
                                     k_part_reco_failed);
  RooExponential d0_part_reco_passed("d0_part_reco_passed",
                                     "d0_part_reco_passed", d0_m_var_shifted,
                                     k_part_reco_passed);

  RooExponential d0_comb_all_failed("d0_comb_all_failed", "d0_comb_all_failed",
                                    d0_m_var_shifted, k_comb_all_failed);
  RooExponential d0_comb_all_passed("d0_comb_all_passed", "d0_comb_all_passed",
                                    d0_m_var_shifted, k_comb_all_passed);

  ///////////////////////////
  // PDFs for MC-only fits //
  ///////////////////////////

  // D0 mass
  RooSimultaneous d0_model_mc(
      "d0_model_mc", "d0_model_mc",
      {{"mc_passed", &d0_model_passed_mc}, {"mc_failed", &d0_model_failed_mc}},
      sample);

  // dm
  RooSimultaneous dm_model_mc(
      "dm_model_mc", "dm_model_mc",
      {{"mc_passed", &dm_model_passed_mc}, {"mc_failed", &dm_model_failed_mc}},
      sample);

  // Full fit

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

  RooRealVar f1_bkg_passed("f1_bkg_passed", "f1_bkg_passed", 0, 1, "");
  RooRealVar f2_bkg_passed("f2_bkg_passed", "f2_bkg_passed", 0, 1, "");
  RooRealVar f1_bkg_failed("f1_bkg_failed", "f1_bkg_failed", 0, 1, "");
  RooRealVar f2_bkg_failed("f2_bkg_failed", "f2_bkg_failed", 0, 1, "");

  RooAddPdf bkg_passed(
      "bkg_passed", "bkg_passed",
      RooArgList(dm_comb_passed, d0_comb_passed, dmd0_comb_passed),
      RooArgList(f1_bkg_passed, f2_bkg_passed), true);
  RooAddPdf bkg_failed(
      "bkg_failed", "bkg_failed",
      RooArgList(dm_comb_failed, d0_comb_failed, dmd0_comb_failed),
      RooArgList(f1_bkg_failed, f2_bkg_failed), true);

  RooRealVar n_mc_passed("n_mc_passed", "n_mc_passed", 0, "");
  RooRealVar n_mc_failed("n_mc_failed", "n_mc_failed", 0, "");

  RooExtendPdf model_mc_passed("model_mc_passed", "model_mc_passed",
                               sig_passed_mc, n_mc_passed);
  RooExtendPdf model_mc_failed("model_mc_failed", "model_mc_failed",
                               sig_failed_mc, n_mc_failed);

  // PID efficiency
  RooRealVar eff("eff", "eff", 0, 1, "");
  if (fake_mu)
    eff.setMin(0.70);
  else
    eff.setMax(0.04);

  // PIDCalib mass window efficiencies
  RooRealVar eff_mw_passed("eff_mw_passed", "eff_mw_passed", 1, "");
  RooRealVar eff_mw_failed("eff_mw_failed", "eff_mw_failed", 1, "");

  // Normalizations
  RooRealVar n_sig("n_sig", "n_sig", 0, 1, "");
  RooRealVar n_bkg_passed("n_bkg_passed", "n_bkg_passed", 0, 1, "");
  RooRealVar n_bkg_failed("n_bkg_failed", "n_bkg_failed", 0, 1, "");

  RooFormulaVar n_sig_passed("n_sig_passed", "n_sig_passed", "x[0]*x[1]*x[2]",
                             RooArgList(n_sig, eff, eff_mw_passed));
  RooFormulaVar n_sig_failed("n_sig_failed", "n_sig_failed",
                             "x[0] * (1. - x[1]) * x[2]",
                             RooArgList(n_sig, eff, eff_mw_failed));

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
                        {{"mc_passed", &model_mc_passed},
                         {"mc_failed", &model_mc_failed},
                         {"calib_passed", &model_passed},
                         {"calib_failed", &model_failed}},
                        sample);

  // Define datasets to fit
  RooDataSet *datasets_mc_passed[N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P];
  RooDataSet *datasets_mc_failed[N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P];
  RooDataSet *datasets_calib_passed[N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P];
  RooDataSet *datasets_calib_failed[N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P];

  auto set_parameters_d0 = [&](const string &probe) {
    if (probe == "k") {
      mean_d0_failed.setVal(1865);
      mean_d0_passed.setVal(1865);
      width_d0.setVal(8.9);
      width_cb_d0.setVal(7.9);
      alpha_L_d0_failed.setVal(1.4);
      alpha_L_d0_passed.setVal(0.8);
      n_L_d0_failed.setVal(3.5);
      n_L_d0_passed.setVal(5);
      alpha_R_d0_failed.setVal(1.5);
      alpha_R_d0_passed.setVal(0.8);
      n_R_d0_failed.setVal(7);
      n_R_d0_passed.setVal(6.5);
      f_d0_passed.setVal(0.5);
      f_d0_failed.setVal(0.5);
    } else if (probe == "pi") {
      mean_d0_failed.setVal(1865.5);
      mean_d0_passed.setVal(1865);
      width_d0.setVal(7.1);
      width_cb_d0.setVal(10.5);
      alpha_L_d0_failed.setVal(1.4);
      alpha_L_d0_passed.setVal(0.6);
      n_L_d0_failed.setVal(5.5);
      n_L_d0_passed.setVal(5);
      alpha_R_d0_failed.setVal(1.7);
      alpha_R_d0_passed.setVal(1.8);
      n_R_d0_failed.setVal(7);
      n_R_d0_passed.setVal(6);
      f_d0_passed.setVal(0.2);
      f_d0_failed.setVal(0.3);
    }
  };

  auto set_parameters_dm = [&](const string &probe) {
    if (probe == "k") {
      mean_dm_failed.setVal(145.4);
      mean_dm_passed.setVal(145.4);
      width_dm.setVal(0.6);
      width_cb_dm.setVal(0.7);
      alpha_L_dm_failed.setVal(1.4);
      alpha_L_dm_passed.setVal(1.2);
      n_L_dm_failed.setVal(5.9);
      n_L_dm_passed.setVal(6.5);
      alpha_R_dm_failed.setVal(1.2);
      alpha_R_dm_passed.setVal(1.0);
      n_R_dm_failed.setVal(5.1);
      n_R_dm_passed.setVal(5.5);
      f_dm_passed.setVal(0.2);
      f_dm_failed.setVal(0.2);
    } else if (probe == "pi") {
      mean_dm_failed.setVal(145.45);
      mean_dm_passed.setVal(145.45);
      width_dm.setVal(0.5);
      width_cb_dm.setVal(0.8);
      alpha_L_dm_failed.setVal(1.6);
      alpha_L_dm_passed.setVal(1.5);
      n_L_dm_failed.setVal(6);
      n_L_dm_passed.setVal(6);
      alpha_R_dm_failed.setVal(1.3);
      alpha_R_dm_passed.setVal(1.2);
      n_R_dm_failed.setVal(5);
      n_R_dm_passed.setVal(5.5);
      f_dm_passed.setVal(0.2);
      f_dm_failed.setVal(0.2);
    }
  };

  auto set_parameters_calib = [&](const string &probe) {
    m_shift_d0_failed.setVal(0);
    m_shift_d0_passed.setVal(0);
    w_scale_d0_failed.setVal(1);
    w_scale_d0_passed.setVal(1);
    m_shift_dm_failed.setVal(0);
    m_shift_dm_passed.setVal(0);
    w_scale_dm_failed.setVal(1);
    w_scale_dm_passed.setVal(1);
    if (probe == "k") {
      c_comb_spi_failed.setVal(0.49);
      c_comb_spi_passed.setVal(0.50);
      c_comb_all_failed.setVal(0.32);
      c_comb_all_passed.setVal(0.40);
    } else if (probe == "pi") {
      c_comb_spi_failed.setVal(0.50);
      c_comb_spi_passed.setVal(0.50);
      c_comb_all_failed.setVal(0.33);
      c_comb_all_passed.setVal(0.35);
    }
  };

  RooArgSet fit_vars(d0_m_var, dm_var);

  RooArgSet argset_minos(eff);

  RooArgSet params_d0(mean_d0_failed, mean_d0_passed, width_d0, width_cb_d0,
                      alpha_L_d0_failed, alpha_L_d0_passed, n_L_d0_failed,
                      n_L_d0_passed, alpha_R_d0_failed, alpha_R_d0_passed,
                      n_R_d0_failed, n_R_d0_passed, f_d0_passed, f_d0_failed);

  RooArgSet params_dm(mean_dm_failed, mean_dm_passed, width_dm, width_cb_dm,
                      alpha_L_dm_failed, alpha_L_dm_passed, n_L_dm_failed,
                      n_L_dm_passed, alpha_R_dm_failed, alpha_R_dm_passed,
                      n_R_dm_failed, n_R_dm_passed, f_dm_passed, f_dm_failed);

  RooRealVar m_shift_d0_delta("m_shift_d0_delta", "m_shift_d0_delta", -5, 5,
                              "");
  RooRealVar m_shift_dm_delta("m_shift_dm_delta", "m_shift_dm_delta", -5, 5,
                              "");
  RooRealVar w_scale_d0_delta("w_scale_d0_delta", "w_scale_d0_delta", -5, 5,
                              "");
  RooRealVar w_scale_dm_delta("w_scale_dm_delta", "w_scale_dm_delta", -5, 5,
                              "");
  RooArgSet  params_calib(
      m_shift_d0_failed, m_shift_d0_passed, m_shift_d0_delta, w_scale_d0_failed,
      w_scale_d0_passed, w_scale_d0_delta, m_shift_dm_failed, m_shift_dm_passed,
      m_shift_dm_delta, w_scale_dm_failed, w_scale_dm_passed, w_scale_dm_delta,
      c_comb_spi_failed, c_comb_spi_passed, c_comb_all_failed,
      c_comb_all_passed, k_part_reco_failed, k_part_reco_passed,
      k_comb_all_failed, k_comb_all_passed, eff, f1_bkg_passed, f2_bkg_passed,
      f1_bkg_failed, f2_bkg_failed);

  for (auto probe : particles) {
    cout << "INFO Selecting " << probe << endl;

    TString tag;
    if (probe == "pi")
      tag = "k";
    else
      tag = "pi";

    // Initialize MC datasets
    cout << "INFO Initializing MC datasets " << endl;
    for (int ntrks_idx = 0; ntrks_idx < N_BINS_NTRACKS; ntrks_idx++) {
      for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
        for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
          TString suffix = TString::Format("%s_%d_%d_%d", probe.c_str(),
                                           ntrks_idx, eta_idx, p_idx);
          datasets_mc_passed[ntrks_idx][eta_idx][p_idx] =
              new RooDataSet("dataset_mc_passed_" + suffix,
                             "dataset_mc_passed_" + suffix, fit_vars);
          datasets_mc_failed[ntrks_idx][eta_idx][p_idx] =
              new RooDataSet("dataset_mc_failed_" + suffix,
                             "dataset_mc_failed_" + suffix, fit_vars);
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
      ch_mc.SetBranchStatus("*", 0);
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
      ch_mc.SetBranchAddress("TrackChi2PerDof", &probe_track_chi2ndof);
      ch_mc.SetBranchAddress("nTracks", &ntracks);

      int count_tm = 0, count_cond = 0, count_calib_sel = 0, count_range = 0,
          count_mw = 0;

      int    kin_bin        = -1;
      bool   in_mass_window = false, pid_ok = false;
      double dm = 0;

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
        // See also https://lhcb-pid-wgp-plots.web.cern.ch/lhcb-pid-wgp-plots/Run2/ and
        // https://gitlab.cern.ch/lhcb-datapkg/WG/SemilepConfig/-/blob/master/options/Filter_Dstar2D02KpiNoPID_2016MC.py?ref_type=heads
        if (tag_p < 2000 || tag_pt < 250) continue;
        if (probe_p < 2000 || probe_pt < 250) continue;
        if (spi_p < 1000 || spi_pt < 100) continue;
        if (probe_track_chi2ndof > 3) continue;
        if (std::max(probe_pt, tag_pt) < 1000 || d0_pt < 1500) continue;
        // TODO Missing track chi2ndof cut on tag and soft pion
        // Will need to reprocess step1 ntuples to include TupleToolTrackInfo
        count_calib_sel++;

        // MuonUnbiased equivalent to L0+HLT1 TIS. See MuonUnbiased defition at
        // https://gitlab.cern.ch/lhcb-datapkg/WG/PIDCalib/-/blob/master/scriptsR2/makeTuples_pp_2016_reprocessing.py#L71
        // and https://mattermost.web.cern.ch/lhcb/pl/893yre484jggigooti5u3gqb8w

        // Not applying MuonUnbiased and GHOSTPROB conditional cuts in MC since
        // some kinematical bins have very low statistics

        if (probe_p < 3000 || probe_p >= 100000 || ntracks >= 600) continue;

        probe_eta = 0.5 * log((probe_p + probe_pz) / (probe_p - probe_pz));
        if (probe_eta < 1.7 || probe_eta >= 5.0) continue;
        count_range++;

        // Determine kinematical bin
        kin_bin = histo_binnning.FindBin(probe_p, probe_eta, ntracks);
        histo_binnning.GetBinXYZ(kin_bin, p_bin, eta_bin, ntrks_bin);

        // Fill histograms
        dm             = dst_m - d0_m;
        in_mass_window = in_range(d0_m_var, d0_m) && in_range(dm_var, dm);
        if (in_mass_window) count_mw++;

        // PID
        if (fake_mu) {
          pid_ok = !probe_ismuon;
        } else if (vmu) {
          pid_ok = probe_ismuon && probe_dllmu > 2.0 && probe_dlle < 1 &&
                   probe_mu_ubdt < 0.25;
        } else {
          pid_ok = probe_ismuon && probe_dllmu > 2.0 && probe_dlle < 1 &&
                   probe_mu_ubdt > 0.25;
        }

        if (pid_ok) {
          count_total_passed[year_idx.at(year)][ntrks_bin - 1][eta_bin - 1]
                            [p_bin - 1]++;
          if (in_mass_window) {
            count_in_mw_passed[year_idx.at(year)][ntrks_bin - 1][eta_bin - 1]
                              [p_bin - 1]++;
            d0_m_var.setVal(d0_m);
            dm_var.setVal(dm);
            datasets_mc_passed[ntrks_bin - 1][eta_bin - 1][p_bin - 1]->addFast(
                fit_vars);
          }
        } else {
          count_total_failed[year_idx.at(year)][ntrks_bin - 1][eta_bin - 1]
                            [p_bin - 1]++;
          if (in_mass_window) {
            count_in_mw_failed[year_idx.at(year)][ntrks_bin - 1][eta_bin - 1]
                              [p_bin - 1]++;
            d0_m_var.setVal(d0_m);
            dm_var.setVal(dm);
            datasets_mc_failed[ntrks_bin - 1][eta_bin - 1][p_bin - 1]->addFast(
                fit_vars);
          }
        }
      }

      auto stop_mc_loop = high_resolution_clock::now();

      const int duration_mc_loop =
          duration_cast<seconds>(stop_mc_loop - start_mc_loop).count();
      cout << "\nMC loop took " << format_time(duration_mc_loop) << endl;

      cout << "\nINFO MC cutflow:" << endl;
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
      cout << "Overall: " << "(" << count_mw * 100. / entries_mc << "%)"
           << endl;
    }

    // Now loop over calib samples and make fits
    for (auto year : years) {
      // TODO Run all years
      if (year != "2016") continue;

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

          bool pid_ok  = false;
          int  kin_bin = -1;

          const int entries_calib = ch_calib.GetEntries();
          cout << "INFO Starting PIDCalib event loop over " << entries_calib
               << " entries in tree " << tree_name << endl;
          auto start_calib_loop = high_resolution_clock::now();
          for (int evt = 0; evt < entries_calib; evt++) {
            ch_calib.GetEntry(evt);

            // Conditional cuts

            if (!probe_hasmuon || probe_ghostprob > 0.5 || !probe_mu_unbiased)
              continue;

            if (probe_p < 3000. || probe_p >= 100000. || ntracks_calib >= 600)
              continue;

            if (probe_eta < 1.7 || probe_eta >= 5.0) continue;

            // Fill histograms

            d0_m_var.setVal(d0_m);
            dm_var.setVal(dst_m - d0_m);

            // Determine kinematical bin
            kin_bin = histo_binnning.FindBin(probe_p, probe_eta, ntracks_calib);

            histo_binnning.GetBinXYZ(kin_bin, p_bin, eta_bin, ntrks_bin);

            if (fake_mu) {
              pid_ok = !probe_ismuon;
            } else if (vmu) {
              pid_ok = probe_ismuon && probe_dllmu > 2.0 && probe_dlle < 1 &&
                       probe_mu_ubdt < 0.25;
            } else {
              pid_ok = probe_ismuon && probe_dllmu > 2.0 && probe_dlle < 1 &&
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
          cout << "\nCalib loop took " << format_time(duration_calib_loop)
               << endl;
        }
      }

      // Create histograms to store efficiency of PIDCalib mass windows
      cout << "INFO Initializing efficiency histograms " << endl;
      TString suffix = TString::Format("%s_%s", year.c_str(), probe.c_str());
      TH2D   *effs_passed =
          static_cast<TH2D *>(histo_binnning.Clone("effs_passed_" + suffix));
      TH2D *effs_passed_unc_hi = static_cast<TH2D *>(
          histo_binnning.Clone("effs_passed_" + suffix + "_unc_hi"));
      TH2D *effs_passed_unc_lo = static_cast<TH2D *>(
          histo_binnning.Clone("effs_passed_" + suffix + "_unc_lo"));
      TH2D *effs_failed =
          static_cast<TH2D *>(histo_binnning.Clone("effs_failed_" + suffix));
      TH2D *effs_failed_unc_hi = static_cast<TH2D *>(
          histo_binnning.Clone("effs_failed_" + suffix + "_unc_hi"));
      TH2D *effs_failed_unc_lo = static_cast<TH2D *>(
          histo_binnning.Clone("effs_failed_" + suffix + "_unc_lo"));

      // Check distribution of fitted parameters
      RooDataSet ds_params_d0("ds_params_d0", "ds_params_d0", params_d0);
      RooDataSet ds_params_dm("ds_params_dm", "ds_params_dm", params_dm);
      RooDataSet ds_params_calib("ds_params_calib", "ds_params_calib",
                                 params_calib);

      TH3D *fit_status_calib = static_cast<TH3D *>(
          histo_binnning.Clone("fit_status_calib_" + suffix));
      fit_status_calib->SetMinimum(-0.1);
      fit_status_calib->SetMaximum(12);
      TH3D *fit_cov_qual_calib = static_cast<TH3D *>(
          histo_binnning.Clone("fit_cov_qual_calib_" + suffix));
      fit_cov_qual_calib->SetMinimum(-0.1);
      fit_cov_qual_calib->SetMaximum(6);
      TH3D *fit_status_d0_mc = static_cast<TH3D *>(
          histo_binnning.Clone("fit_status_d0_mc_" + suffix));
      fit_status_d0_mc->SetMinimum(-0.1);
      fit_status_d0_mc->SetMaximum(12);
      TH3D *fit_cov_qual_d0_mc = static_cast<TH3D *>(
          histo_binnning.Clone("fit_cov_qual_d0_mc_" + suffix));
      fit_cov_qual_d0_mc->SetMinimum(-0.1);
      fit_cov_qual_d0_mc->SetMaximum(6);
      TH3D *fit_status_dm_mc = static_cast<TH3D *>(
          histo_binnning.Clone("fit_status_dm_mc_" + suffix));
      fit_status_dm_mc->SetMinimum(-0.1);
      fit_status_dm_mc->SetMaximum(12);
      TH3D *fit_cov_qual_dm_mc = static_cast<TH3D *>(
          histo_binnning.Clone("fit_cov_qual_dm_mc_" + suffix));
      fit_cov_qual_dm_mc->SetMinimum(-0.1);
      fit_cov_qual_dm_mc->SetMaximum(6);
      TH3D *fit_status =
          static_cast<TH3D *>(histo_binnning.Clone("fit_status"));
      fit_status->SetMinimum(-0.1);
      fit_status->SetMaximum(12);
      TH3D *fit_cov_qual =
          static_cast<TH3D *>(histo_binnning.Clone("fit_cov_qual"));
      fit_cov_qual->SetMinimum(-0.1);
      fit_cov_qual->SetMaximum(6);

      TH1D d0_retries("d0_retries", "d0_m MC fit #retries;#;;",
                      max_fix_reattempts + 1, 0, max_fix_reattempts + 1);
      TH1D dm_retries("dm_retries", "dm MC fit #retries;#;;",
                      max_fix_reattempts + 1, 0, max_fix_reattempts + 1);
      TH1D calib_retries("calib_retries", "calib fit #retries;#;;",
                         max_fix_reattempts + 1, 0, max_fix_reattempts + 1);

      ofile.cd();

      cout << "INFO Calculating efficiencies " << endl;

      // Create histogram to hold PID efficiency
      // Use "eff" as name to match PIDCalib output
      unique_ptr<TH3D> histo_pid(
          static_cast<TH3D *>(histo_binnning.Clone("eff")));

      for (int ntrks_idx = 0; ntrks_idx < N_BINS_NTRACKS; ntrks_idx++) {
        for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
          for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
            suffix.Form("%s_%s_%d_%d_%d", year.c_str(), probe.c_str(),
                        ntrks_idx, eta_idx, p_idx);
            const TString tag = TString::Format(
                "(%.0f < nTracks < %.0f, %.1f < #eta < %.1f, %.0f < p < %.0f)",
                BINS_NTRACKS[ntrks_idx], BINS_NTRACKS[ntrks_idx + 1],
                BINS_ETA[eta_idx], BINS_ETA[eta_idx + 1], BINS_P[p_idx],
                BINS_P[p_idx + 1]);

            // Get p, eta bin
            const int kin_bin =
                effs_passed->GetBin(p_idx + 1, eta_idx + 1, ntrks_idx + 1);

            // Get correction for "passed" events
            const double in_mw_passed =
                count_in_mw_passed[year_idx.at(year)][ntrks_idx][eta_idx]
                                  [p_idx];
            const double total_passed =
                count_total_passed[year_idx.at(year)][ntrks_idx][eta_idx]
                                  [p_idx];
            const double eff_passed        = in_mw_passed / total_passed;
            eff_mw_passed                  = eff_passed;
            const double eff_passed_unc_hi = TEfficiency::Wilson(
                total_passed, in_mw_passed, ONE_SIGMA, true);
            const double eff_passed_unc_lo = TEfficiency::Wilson(
                total_passed, in_mw_passed, ONE_SIGMA, false);
            effs_passed->SetBinContent(kin_bin, eff_passed);
            effs_passed_unc_hi->SetBinContent(kin_bin, eff_passed_unc_hi);
            effs_passed_unc_lo->SetBinContent(kin_bin, eff_passed_unc_lo);

            // Now get corrections for "failed" events
            const double in_mw_failed =
                count_in_mw_failed[year_idx.at(year)][ntrks_idx][eta_idx]
                                  [p_idx];
            const double total_failed =
                count_total_failed[year_idx.at(year)][ntrks_idx][eta_idx]
                                  [p_idx];
            const double eff_failed        = in_mw_failed / total_failed;
            eff_mw_failed                  = eff_failed;
            const double eff_failed_unc_hi = TEfficiency::Wilson(
                total_failed, in_mw_failed, ONE_SIGMA, true);
            const double eff_failed_unc_lo = TEfficiency::Wilson(
                total_failed, in_mw_failed, ONE_SIGMA, false);
            effs_failed->SetBinContent(kin_bin, eff_failed);
            effs_failed_unc_hi->SetBinContent(kin_bin, eff_failed_unc_hi);
            effs_failed_unc_lo->SetBinContent(kin_bin, eff_failed_unc_lo);

            // Print mass-window efficiencies
            cout << "\nINFO Mass-window efficiency for passed sample: "
                 << in_mw_passed << " / " << total_passed << " = " << eff_passed
                 << endl;
            cout << "\nINFO Mass-window efficiency for failed sample: "
                 << in_mw_failed << " / " << total_failed << " = " << eff_failed
                 << endl;

            //////////////////////////////////////////////////////////
            // Simultaneous fits to passed and failed distributions //
            //////////////////////////////////////////////////////////

            mean_d0_failed.setConstant(false);
            mean_d0_passed.setConstant(false);
            width_d0.setConstant(false);
            width_cb_d0.setConstant(false);
            alpha_L_d0_passed.setConstant(false);
            alpha_R_d0_passed.setConstant(false);
            n_L_d0_passed.setConstant(false);
            n_R_d0_passed.setConstant(false);
            alpha_L_d0_failed.setConstant(false);
            alpha_R_d0_failed.setConstant(false);
            n_L_d0_failed.setConstant(false);
            n_R_d0_failed.setConstant(false);
            f_d0_passed.setConstant(false);
            f_d0_failed.setConstant(false);
            mean_dm_failed.setConstant(false);
            mean_dm_passed.setConstant(false);
            width_dm.setConstant(false);
            width_cb_dm.setConstant(false);
            alpha_L_dm_passed.setConstant(false);
            alpha_R_dm_passed.setConstant(false);
            n_L_dm_passed.setConstant(false);
            n_R_dm_passed.setConstant(false);
            alpha_L_dm_failed.setConstant(false);
            alpha_R_dm_failed.setConstant(false);
            n_L_dm_failed.setConstant(false);
            n_R_dm_failed.setConstant(false);
            f_dm_passed.setConstant(false);
            f_dm_failed.setConstant(false);

            // Build combined dataset
            auto &dataset_mc_passed =
                datasets_mc_passed[ntrks_idx][eta_idx][p_idx];
            auto &dataset_mc_failed =
                datasets_mc_failed[ntrks_idx][eta_idx][p_idx];
            auto &dataset_calib_passed =
                datasets_calib_passed[ntrks_idx][eta_idx][p_idx];
            auto &dataset_calib_failed =
                datasets_calib_failed[ntrks_idx][eta_idx][p_idx];

            double       n_calib_passed    = dataset_calib_passed->numEntries();
            double       n_calib_failed    = dataset_calib_failed->numEntries();
            const double n_mc_guess_passed = dataset_mc_passed->numEntries();
            const double n_mc_guess_failed = dataset_mc_failed->numEntries();

            // Plot 2D distributions in calib sample
            TH2F *th2_calib_passed = dataset_calib_passed->createHistogram(
                d0_m_var, dm_var, 40, 40, "", "th2_calib_passed");
            TH2F *th2_calib_failed = dataset_calib_failed->createHistogram(
                d0_m_var, dm_var, 100, 100, "", "th2_calib_failed");
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
            cout << "  - " << n_mc_guess_passed << " passed MC events" << endl;
            cout << "  - " << n_mc_guess_failed << " failed MC events" << endl;

            RooDataSet dataset("dataset", "combined calib + MC data", fit_vars,
                               Index(sample),
                               Import("mc_passed", *dataset_mc_passed),
                               Import("mc_failed", *dataset_mc_failed),
                               Import("calib_passed", *dataset_calib_passed),
                               Import("calib_failed", *dataset_calib_failed));

            cout << "INFO Built combined dataset with " << dataset.numEntries()
                 << " calib + MC events" << endl;

            // d0_m fit
            // Use "minimize" instead of default "migrad".
            // "mimimize" calls MIGRAD, and in case of failure calls
            // SIMPLEX, before calling MIGRAD again. Has been quite helpful
            // (just search for SIMPLEX in the log).
            cout << "\nINFO Starting D0 MC fit " << suffix << endl;
            set_parameters_d0(probe);
            RooFitResult *fit_res_d0_mc = d0_model_mc.fitTo(
                dataset, Strategy(2), NumCPU(8), PrintLevel(0), Save(),
                Offset(), Minimizer("Minuit", "minimize"), Timer());

            int d0_fit_reattempts = 0;
            // If the fit does not converge, keep alternating simplex and
            // migrad until it converges or hits the limit of retries. Here we
            // hope that simplex can at least find a better starting point for
            // migrad.
            // TODO Make fits work with less retries (ideally none!), as those
            // take time
            while (!fit_ok(fit_res_d0_mc) &&
                   (d0_fit_reattempts < max_fix_reattempts)) {
              d0_fit_reattempts++;
              cout << "\nINFO Retrying D0 MC fit " << suffix << " (retry #"
                   << d0_fit_reattempts << ")" << endl;
              cout << "INFO Previous fit status: " << fit_res_d0_mc->status()
                   << endl;
              cout << "INFO Previous cov matrix status: "
                   << fit_res_d0_mc->covQual() << endl;
              d0_model_mc.fitTo(dataset, NumCPU(8), PrintLevel(0), Offset(),
                                Minimizer("Minuit", "scan"), Timer(),
                                Hesse(false));
              fit_res_d0_mc = d0_model_mc.fitTo(
                  dataset, Strategy(2), NumCPU(8), PrintLevel(0), Save(),
                  Offset(), Minimizer("Minuit", "minimize"), Timer());
            }

            d0_retries.Fill(d0_fit_reattempts);

            if (fit_ok(fit_res_d0_mc)) {
              ds_params_d0.addFast(params_d0);
            }

            fit_status_d0_mc->SetBinContent(kin_bin, fit_res_d0_mc->status());
            fit_cov_qual_d0_mc->SetBinContent(kin_bin,
                                              fit_res_d0_mc->covQual());

            // Plot fit results
            RooPlot *frame_d0_passed =
                d0_m_var.frame(Title("D0 M Passed " + tag));
            RooPlot *frame_d0_failed =
                d0_m_var.frame(Title("D0 M Failed " + tag));
            dataset.plotOn(frame_d0_passed, Binning(bins_d0_m_passed),
                           Cut("sample==sample::mc_passed"));
            dataset.plotOn(frame_d0_failed, Binning(bins_d0_m),
                           Cut("sample==sample::mc_failed"));
            d0_model_mc.plotOn(frame_d0_passed, Slice(sample, "mc_passed"),
                               ProjWData(sample, dataset));
            d0_model_mc.plotOn(frame_d0_failed, Slice(sample, "mc_failed"),
                               ProjWData(sample, dataset));

            c_double.cd(1);
            frame_d0_passed->Draw();
            c_double.cd(2);
            frame_d0_failed->Draw();
            c_double.SaveAs(opath + "/figs/fits/d0_m_" + suffix + ".pdf");

            // dm fit
            cout << "\nINFO Starting dm MC fit " << suffix << endl;
            set_parameters_dm(probe);
            RooFitResult *fit_res_dm_mc = dm_model_mc.fitTo(
                dataset, Strategy(2), NumCPU(8), Save(), PrintLevel(0),
                Offset(), Minimizer("Minuit", "minimize"), Timer());

            int dm_fit_reattempts = 0;

            while (!fit_ok(fit_res_dm_mc) &&
                   (dm_fit_reattempts < max_fix_reattempts)) {
              dm_fit_reattempts++;
              cout << "\nINFO Retrying dm MC fit " << suffix << " (retry #"
                   << dm_fit_reattempts << ")" << endl;
              cout << "INFO Previous fit status: " << fit_res_dm_mc->status()
                   << endl;
              cout << "INFO Previous cov matrix status: "
                   << fit_res_dm_mc->covQual() << endl;
              dm_model_mc.fitTo(dataset, NumCPU(8), PrintLevel(0), Offset(),
                                Minimizer("Minuit", "scan"), Timer(),
                                Hesse(false));
              fit_res_dm_mc = dm_model_mc.fitTo(
                  dataset, Strategy(2), NumCPU(8), PrintLevel(0), Save(),
                  Offset(), Minimizer("Minuit", "minimize"), Timer());
            }

            dm_retries.Fill(dm_fit_reattempts);

            if (fit_ok(fit_res_dm_mc)) {
              ds_params_dm.addFast(params_dm);
            }

            fit_status_dm_mc->SetBinContent(kin_bin, fit_res_dm_mc->status());
            fit_cov_qual_dm_mc->SetBinContent(kin_bin,
                                              fit_res_dm_mc->covQual());

            // Plot fit results
            RooPlot *frame_dm_passed = dm_var.frame(Title("dm Passed " + tag));
            RooPlot *frame_dm_failed = dm_var.frame(Title("dm Failed " + tag));
            dataset.plotOn(frame_dm_passed, Binning(bins_dm),
                           Cut("sample==sample::mc_passed"));
            dataset.plotOn(frame_dm_failed, Binning(bins_dm),
                           Cut("sample==sample::mc_failed"));
            dm_model_mc.plotOn(frame_dm_passed, Slice(sample, "mc_passed"),
                               ProjWData(sample, dataset));
            dm_model_mc.plotOn(frame_dm_failed, Slice(sample, "mc_failed"),
                               ProjWData(sample, dataset));

            c_double.cd(1);
            frame_dm_passed->Draw();
            c_double.cd(2);
            frame_dm_failed->Draw();
            c_double.SaveAs(opath + "/figs/fits/dm_" + suffix + ".pdf");

            mean_d0_failed.setConstant();
            mean_d0_passed.setConstant();
            width_d0.setConstant();
            width_cb_d0.setConstant();
            alpha_L_d0_passed.setConstant();
            alpha_R_d0_passed.setConstant();
            n_L_d0_passed.setConstant();
            n_R_d0_passed.setConstant();
            alpha_L_d0_failed.setConstant();
            alpha_R_d0_failed.setConstant();
            n_L_d0_failed.setConstant();
            n_R_d0_failed.setConstant();
            mean_dm_failed.setConstant();
            mean_dm_passed.setConstant();
            width_dm.setConstant();
            width_cb_dm.setConstant();
            alpha_L_dm_passed.setConstant();
            alpha_R_dm_passed.setConstant();
            n_L_dm_passed.setConstant();
            n_R_dm_passed.setConstant();
            alpha_L_dm_failed.setConstant();
            alpha_R_dm_failed.setConstant();
            n_L_dm_failed.setConstant();
            n_R_dm_failed.setConstant();
            f_d0_passed.setConstant();
            f_d0_failed.setConstant();
            f_dm_passed.setConstant();
            f_dm_failed.setConstant();

            // Fit only calib sample

            const TString sb_cut =
                "(dm_var > 149) || (dm_var < 143) || (d0_m_var > 1890) || "
                "(d0_m_var < 1840)";
            constexpr double A_sig = (1890 - 1840) * (149 - 143);
            constexpr double A_fit = (1910 - 1825) * (153 - 141);
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
                "((dm_var > 149) || (dm_var < 143)) && ((d0_m_var < 1890) && "
                "(d0_m_var > 1840))";
            const TString sb2_cut =
                "((dm_var < 149) && (dm_var > 143)) && ((d0_m_var > 1890) || "
                "(d0_m_var < 1840))";
            const TString sb3_cut =
                "((dm_var > 149) || (dm_var < 143)) && ((d0_m_var > 1890) || "
                "(d0_m_var < 1840))";
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

            constexpr double A_sb3 =
                (153 - 149) * (1840 - 1825) + (153 - 149) * (1910 - 1890) +
                (143 - 141) * (1840 - 1825) + (143 - 141) * (1910 - 1890);
            constexpr double A_sb2 =
                (149 - 143) * (1840 - 1825) + (149 - 143) * (1910 - 1890);
            constexpr double A_sb1 =
                (153 - 149) * (1890 - 1840) + (143 - 141) * (1890 - 1840);

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
            n_mc_passed.setVal(n_mc_guess_passed);
            n_mc_failed.setVal(n_mc_guess_failed);

            // Estimate exponential coefficients
            const double dx = 1897.5 - 1827.5;

            const double y1_passed = dataset_calib_passed->sumEntries(
                "(d0_m_var > 1825) && (d0_m_var < 1830)");
            const double y2_passed = dataset_calib_passed->sumEntries(
                "(d0_m_var > 1895) && (d0_m_var < 1900)");
            const double k_passed_guess = log(y2_passed / y1_passed) / dx;
            if (k_passed_guess < 0) {
              k_comb_all_passed.setVal(k_passed_guess);
              k_part_reco_passed.setVal(k_passed_guess);
            } else {
              k_comb_all_passed.setVal(-1e-4);
              k_part_reco_passed.setVal(-1e-4);
            }

            const double y1_failed = dataset_calib_failed->sumEntries(
                "(d0_m_var > 1825) && (d0_m_var < 1830)");
            const double y2_failed = dataset_calib_failed->sumEntries(
                "(d0_m_var > 1895) && (d0_m_var < 1900)");
            const double k_failed_guess = log(y2_failed / y1_failed) / dx;
            if (k_passed_guess < 0) {
              k_comb_all_failed.setVal(k_failed_guess);
              k_part_reco_failed.setVal(k_failed_guess);
            } else {
              k_comb_all_failed.setVal(-1e-4);
              k_part_reco_failed.setVal(-1e-4);
            }

            // Fine tuning
            if ((probe == "pi") && (p_idx <= 1)) {
              // Use reduced fit range excluding higher D0 mass region to avoid
              // mass threshold in certain kinematic bins
              d0_m_var.setRange("fit", 1825, 1900);
              dm_var.setRange("fit", 141, 153);
            } else {
              // Use full fit range
              d0_m_var.setRange("fit", 1825, 1910);
              dm_var.setRange("fit", 141, 153);
            }

            //

            RooCmdArg fit_strat  = Strategy(2);
            RooCmdArg numcpu     = NumCPU(8);
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
            fit_args.Add(&minos);

            cout << "\n\nINFO Starting calib only fit " << suffix << endl;

            set_parameters_calib(probe);

            if (!vmu && !fake_mu) {
              if ((probe == "pi") && (p_idx < 1)) {
                w_scale_d0_passed.setConstant();
                m_shift_d0_passed.setConstant();
                w_scale_d0_passed.setVal(1);
                m_shift_d0_passed.setVal(0);
              } else {
                w_scale_d0_passed.setConstant(false);
                m_shift_d0_passed.setConstant(false);
              }
            }

            // if (probe == "k") {
            //   f2_bkg_failed.setConstant();
            //   f2_bkg_passed.setConstant();
            //   k_part_reco_failed.setConstant();
            //   k_part_reco_passed.setConstant();
            //   f2_bkg_failed.setVal(0);
            //   f2_bkg_passed.setVal(0);
            // } else {
            //   f2_bkg_failed.setConstant(false);
            //   f2_bkg_passed.setConstant(false);
            //   k_part_reco_failed.setConstant(false);
            //   k_part_reco_passed.setConstant(false);
            // }

            RooFitResult *fit_res_calib = model_calib.fitTo(dataset, fit_args);

            int calib_fit_reattempts = 0;
            while (!fit_ok(fit_res_calib) &&
                   (calib_fit_reattempts < max_fix_reattempts)) {
              calib_fit_reattempts++;
              cout << "\n\nINFO Retrying calib only fit " << suffix
                   << " (retry #" << calib_fit_reattempts << ")" << endl;
              cout << "INFO Previous fit status: " << fit_res_calib->status()
                   << endl;
              cout << "INFO Previous cov matrix status: "
                   << fit_res_calib->covQual() << endl;

              model_calib.fitTo(dataset, Extended(), NumCPU(8), PrintLevel(0),
                                Offset(), Range("fit"),
                                Minimizer("Minuit", "scan"), Timer(),
                                Hesse(false));
              fit_res_calib = model_calib.fitTo(dataset, fit_args);
            }

            calib_retries.Fill(calib_fit_reattempts);

            if (fit_ok(fit_res_calib)) {
              m_shift_d0_delta =
                  -(m_shift_d0_passed.getVal() - m_shift_d0_failed.getVal()) /
                  (sqrt(m_shift_d0_passed.getError() *
                            m_shift_d0_passed.getError() +
                        m_shift_d0_failed.getError() *
                            m_shift_d0_failed.getError()));

              m_shift_dm_delta =
                  -(m_shift_dm_passed.getVal() - m_shift_dm_failed.getVal()) /
                  (sqrt(m_shift_dm_passed.getError() *
                            m_shift_dm_passed.getError() +
                        m_shift_dm_failed.getError() *
                            m_shift_dm_failed.getError()));

              w_scale_d0_delta =
                  -(w_scale_d0_passed.getVal() - w_scale_d0_failed.getVal()) /
                  (sqrt(w_scale_d0_passed.getError() *
                            w_scale_d0_passed.getError() +
                        w_scale_d0_failed.getError() *
                            w_scale_d0_failed.getError()));

              w_scale_dm_delta =
                  -(w_scale_dm_passed.getVal() - w_scale_dm_failed.getVal()) /
                  (sqrt(w_scale_dm_passed.getError() *
                            w_scale_dm_passed.getError() +
                        w_scale_dm_failed.getError() *
                            w_scale_dm_failed.getError()));

              ds_params_calib.addFast(params_calib);
            }

            fit_status_calib->SetBinContent(kin_bin, fit_res_calib->status());
            fit_cov_qual_calib->SetBinContent(kin_bin,
                                              fit_res_calib->covQual());

            // Plot fit results
            RooPlot *frame_d0_calib_only_passed =
                d0_m_var.frame(Title("D0 M Calib Passed " + tag));
            RooPlot *frame_d0_calib_only_failed =
                d0_m_var.frame(Title("D0 M Calib Failed " + tag));
            RooPlot *frame_dm_calib_only_passed =
                dm_var.frame(Title("dm Calib Passed " + tag));
            RooPlot *frame_dm_calib_only_failed =
                dm_var.frame(Title("dm Calib Failed " + tag));

            // TODO
            // https://root-forum.cern.ch/t/simultaneous-fit-normalization-issue/33965
            dataset.plotOn(frame_d0_calib_only_passed,
                           Binning(bins_d0_m_passed),
                           Cut("sample==sample::calib_passed"));
            dataset.plotOn(frame_d0_calib_only_failed, Binning(bins_d0_m),
                           Cut("sample==sample::calib_failed"));
            dataset.plotOn(frame_dm_calib_only_passed, Binning(bins_dm_passed),
                           Cut("sample==sample::calib_passed"),
                           CutRange("fit"));
            dataset.plotOn(frame_dm_calib_only_failed, Binning(bins_dm),
                           Cut("sample==sample::calib_failed"),
                           CutRange("fit"));

            model_calib.plotOn(frame_d0_calib_only_passed,
                               Slice(sample, "calib_passed"),
                               ProjWData(sample, dataset), LineWidth(1));
            model_calib.plotOn(frame_d0_calib_only_passed,
                               Slice(sample, "calib_passed"),
                               ProjWData(sample, dataset), LineWidth(1),
                               Components(d0_comb_passed), LineColor(kMagenta));
            model_calib.plotOn(frame_d0_calib_only_passed,
                               Slice(sample, "calib_passed"),
                               ProjWData(sample, dataset), LineWidth(1),
                               Components(dm_comb_passed), LineColor(kGreen));
            model_calib.plotOn(frame_d0_calib_only_passed,
                               Slice(sample, "calib_passed"),
                               ProjWData(sample, dataset), LineWidth(1),
                               Components(dmd0_comb_passed), LineColor(kRed));

            model_calib.plotOn(frame_d0_calib_only_failed,
                               Slice(sample, "calib_failed"),
                               ProjWData(sample, dataset), LineWidth(1));
            model_calib.plotOn(frame_d0_calib_only_failed,
                               Slice(sample, "calib_failed"),
                               ProjWData(sample, dataset), LineWidth(1),
                               Components(d0_comb_failed), LineColor(kMagenta));
            model_calib.plotOn(frame_d0_calib_only_failed,
                               Slice(sample, "calib_failed"),
                               ProjWData(sample, dataset), LineWidth(1),
                               Components(dm_comb_failed), LineColor(kGreen));
            model_calib.plotOn(frame_d0_calib_only_failed,
                               Slice(sample, "calib_failed"),
                               ProjWData(sample, dataset), LineWidth(1),
                               Components(dmd0_comb_failed), LineColor(kRed));

            model_calib.plotOn(frame_dm_calib_only_passed,
                               Slice(sample, "calib_passed"),
                               ProjWData(sample, dataset),
                               ProjectionRange("fit"), LineWidth(1));
            model_calib.plotOn(
                frame_dm_calib_only_passed, Slice(sample, "calib_passed"),
                ProjWData(sample, dataset), ProjectionRange("fit"),
                LineWidth(1), Components(d0_comb_passed), LineColor(kMagenta));
            model_calib.plotOn(
                frame_dm_calib_only_passed, Slice(sample, "calib_passed"),
                ProjWData(sample, dataset), ProjectionRange("fit"),
                LineWidth(1), Components(dm_comb_passed), LineColor(kGreen));
            model_calib.plotOn(
                frame_dm_calib_only_passed, Slice(sample, "calib_passed"),
                ProjWData(sample, dataset), ProjectionRange("fit"),
                LineWidth(1), Components(dmd0_comb_passed), LineColor(kRed));

            model_calib.plotOn(frame_dm_calib_only_failed,
                               Slice(sample, "calib_failed"),
                               ProjWData(sample, dataset),
                               ProjectionRange("fit"), LineWidth(1));
            model_calib.plotOn(
                frame_dm_calib_only_failed, Slice(sample, "calib_failed"),
                ProjWData(sample, dataset), ProjectionRange("fit"),
                LineWidth(1), Components(d0_comb_failed), LineColor(kMagenta));
            model_calib.plotOn(
                frame_dm_calib_only_failed, Slice(sample, "calib_failed"),
                ProjWData(sample, dataset), ProjectionRange("fit"),
                LineWidth(1), Components(dm_comb_failed), LineColor(kGreen));
            model_calib.plotOn(
                frame_dm_calib_only_failed, Slice(sample, "calib_failed"),
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

              mean_d0_failed.setConstant(false);
              mean_d0_passed.setConstant(false);
              width_d0.setConstant(false);
              width_cb_d0.setConstant(false);
              alpha_L_d0_passed.setConstant(false);
              alpha_R_d0_passed.setConstant(false);
              n_L_d0_passed.setConstant(false);
              n_R_d0_passed.setConstant(false);
              alpha_L_d0_failed.setConstant(false);
              alpha_R_d0_failed.setConstant(false);
              n_L_d0_failed.setConstant(false);
              n_R_d0_failed.setConstant(false);
              mean_dm_failed.setConstant(false);
              mean_dm_passed.setConstant(false);
              width_dm.setConstant(false);
              width_cb_dm.setConstant(false);
              alpha_L_dm_passed.setConstant(false);
              alpha_R_dm_passed.setConstant(false);
              n_L_dm_passed.setConstant(false);
              n_R_dm_passed.setConstant(false);
              alpha_L_dm_failed.setConstant(false);
              alpha_R_dm_failed.setConstant(false);
              n_L_dm_failed.setConstant(false);
              n_R_dm_failed.setConstant(false);
              f_d0_passed.setConstant(false);
              f_d0_failed.setConstant(false);
              f_dm_passed.setConstant(false);
              f_dm_failed.setConstant(false);

              RooFitResult *fit_res = model.fitTo(dataset, fit_args);

              int full_fit_reattempts = 0;
              while (!fit_ok(fit_res) &&
                     (full_fit_reattempts < max_fix_reattempts)) {
                full_fit_reattempts++;
                cout << "\n\nINFO Retrying full fit " << suffix << " (retry #"
                     << full_fit_reattempts << ")" << endl;
                cout << "INFO Previous fit status: " << fit_res->status()
                     << endl;
                cout << "INFO Previous cov matrix status: "
                     << fit_res->covQual() << endl;
                model.fitTo(dataset, Extended(), NumCPU(8), PrintLevel(0),
                            Offset(), Range("fit"), Minimizer("Minuit", "scan"),
                            Timer(), Hesse(false));
                fit_res = model_calib.fitTo(dataset, fit_args);
              }

              fit_status->SetBinContent(kin_bin, fit_res->status());
              fit_cov_qual->SetBinContent(kin_bin, fit_res->covQual());
            }

            histo_pid->SetBinContent(kin_bin, eff.getVal());
            // Asymmetric eff errors are calculated with MINOS. As in our
            // pidcalib wrapper, use largest as symmetric error.
            const double eff_error = std::max(std::abs(eff.getErrorHi()),
                                              std::abs(eff.getErrorLo()));
            histo_pid->SetBinError(kin_bin, eff_error);

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
            cout << " - Fitted k_comb_passed = " << k_comb_all_passed.getVal()
                 << " vs estimated " << k_passed_guess << " ("
                 << k_comb_all_passed.getVal() / k_passed_guess << ")" << endl;
            cout << " - Fitted k_part_reco_failed = "
                 << k_part_reco_failed.getVal() << " vs estimated "
                 << k_failed_guess << " ("
                 << k_part_reco_failed.getVal() / k_failed_guess << ")" << endl;
            cout << " - Fitted k_comb_failed = " << k_comb_all_failed.getVal()
                 << " vs estimated " << k_failed_guess << " ("
                 << k_comb_all_failed.getVal() / k_failed_guess << ")\n"
                 << endl;

            // Plot fit results
            RooPlot *frame_d0_mc_passed =
                d0_m_var.frame(Title("D0 M MC Passed " + tag));
            RooPlot *frame_d0_calib_passed =
                d0_m_var.frame(Title("D0 M Calib Passed " + tag));
            RooPlot *frame_d0_mc_failed =
                d0_m_var.frame(Title("D0 M MC Failed " + tag));
            RooPlot *frame_d0_calib_failed =
                d0_m_var.frame(Title("D0 M Calib Failed " + tag));
            RooPlot *frame_dm_mc_passed =
                dm_var.frame(Title("dm MC Passed " + tag));
            RooPlot *frame_dm_calib_passed =
                dm_var.frame(Title("dm Calib Passed " + tag));
            RooPlot *frame_dm_mc_failed =
                dm_var.frame(Title("dm MC Failed " + tag));
            RooPlot *frame_dm_calib_failed =
                dm_var.frame(Title("dm Calib Failed " + tag));

            dataset.plotOn(frame_d0_mc_passed, Binning(bins_d0_m_passed),
                           Cut("sample==sample::mc_passed"));
            dataset.plotOn(frame_d0_mc_failed, Binning(bins_d0_m),
                           Cut("sample==sample::mc_failed"));
            dataset.plotOn(frame_d0_calib_passed, Binning(bins_d0_m_passed),
                           Cut("sample==sample::calib_passed"));
            dataset.plotOn(frame_d0_calib_failed, Binning(bins_d0_m),
                           Cut("sample==sample::calib_failed"));
            dataset.plotOn(frame_dm_mc_passed, Binning(bins_dm_passed),
                           Cut("sample==sample::mc_passed"));
            dataset.plotOn(frame_dm_mc_failed, Binning(bins_dm),
                           Cut("sample==sample::mc_failed"));
            dataset.plotOn(frame_dm_calib_passed, Binning(bins_dm_passed),
                           Cut("sample==sample::calib_passed"),
                           CutRange("fit"));
            dataset.plotOn(frame_dm_calib_failed, Binning(bins_dm),
                           Cut("sample==sample::calib_failed"),
                           CutRange("fit"));

            d0_model_mc.plotOn(frame_d0_mc_passed, Slice(sample, "mc_passed"),
                               ProjWData(sample, dataset), LineWidth(1));
            d0_model_mc.plotOn(frame_d0_mc_failed, Slice(sample, "mc_failed"),
                               ProjWData(sample, dataset), LineWidth(1));
            model_calib.plotOn(frame_d0_calib_passed,
                               Slice(sample, "calib_passed"),
                               ProjWData(sample, dataset), LineWidth(1));
            model_calib.plotOn(frame_d0_calib_failed,
                               Slice(sample, "calib_failed"),
                               ProjWData(sample, dataset), LineWidth(1));
            dm_model_mc.plotOn(frame_dm_mc_passed, Slice(sample, "mc_passed"),
                               ProjWData(sample, dataset), LineWidth(1));
            dm_model_mc.plotOn(frame_dm_mc_failed, Slice(sample, "mc_failed"),
                               ProjWData(sample, dataset), LineWidth(1));
            model_calib.plotOn(frame_dm_calib_passed,
                               Slice(sample, "calib_passed"),
                               ProjWData(sample, dataset),
                               ProjectionRange("fit"), LineWidth(1));
            model_calib.plotOn(frame_dm_calib_failed,
                               Slice(sample, "calib_failed"),
                               ProjWData(sample, dataset),
                               ProjectionRange("fit"), LineWidth(1));

            c_mult.cd(1);
            frame_d0_mc_passed->Draw();
            c_mult.cd(2);
            frame_dm_mc_passed->Draw();
            c_mult.cd(3);
            frame_d0_calib_passed->Draw();
            c_mult.cd(4);
            frame_dm_calib_passed->Draw();
            c_mult.cd(5);
            frame_d0_mc_failed->Draw();
            c_mult.cd(6);
            frame_dm_mc_failed->Draw();
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
      histo_pid->Write();
      ofile_probe.Close();

      // Now save other distributions in common output file
      ofile.cd();
      fit_status_d0_mc->Write();
      fit_cov_qual_d0_mc->Write();
      fit_status_dm_mc->Write();
      fit_cov_qual_dm_mc->Write();
      fit_status_calib->Write();
      fit_cov_qual_calib->Write();

      // Save fit status
      suffix.Form("%s_%s", year.c_str(), probe.c_str());
      c_single.cd();
      fit_status_d0_mc->Draw("BOX2Z");
      c_single.SaveAs(opath + "/figs/fit_status_d0_mc_" + suffix + ".pdf");
      fit_cov_qual_d0_mc->Draw("BOX2Z");
      c_single.SaveAs(opath + "/figs/fit_cov_qual_d0_mc_" + suffix + ".pdf");
      fit_status_dm_mc->Draw("BOX2Z");
      c_single.SaveAs(opath + "/figs/fit_status_dm_mc_" + suffix + ".pdf");
      fit_cov_qual_dm_mc->Draw("BOX2Z");
      c_single.SaveAs(opath + "/figs/fit_cov_qual_dm_mc_" + suffix + ".pdf");
      fit_status_calib->Draw("BOX2Z");
      c_single.SaveAs(opath + "/figs/fit_status_calib_" + suffix + ".pdf");
      fit_cov_qual_calib->Draw("BOX2Z");
      c_single.SaveAs(opath + "/figs/fit_cov_qual_calib_" + suffix + ".pdf");
      if (run_simultaneous_mc_calib) {
        fit_status->Draw("BOX2Z");
        c_single.SaveAs(opath + "/figs/fit_status_" + suffix + ".pdf");
        fit_cov_qual->Draw("BOX2Z");
        c_single.SaveAs(opath + "/figs/fit_cov_qual_" + suffix + ".pdf");
      }

      // Save number of fit retries
      d0_retries.Draw();
      c_single.SaveAs(opath + "/figs/fit_retries_d0_mc_" + suffix + ".pdf");
      dm_retries.Draw();
      c_single.SaveAs(opath + "/figs/fit_retries_dm_mc_" + suffix + ".pdf");
      calib_retries.Draw();
      c_single.SaveAs(opath + "/figs/fit_retries_calib_" + suffix + ".pdf");

      // Plot distribution of fitted variables
      plot_dataset(ds_params_d0, opath + "/figs/params/" + suffix + "_");
      plot_dataset(ds_params_dm, opath + "/figs/params/" + suffix + "_");
      plot_dataset(ds_params_calib, opath + "/figs/params/" + suffix + "_");

      // Save efficiencies
      cout << "INFO Saving efficiencies " << endl;
      effs_passed->Write();
      effs_passed_unc_hi->Write();
      effs_passed_unc_lo->Write();
      effs_failed->Write();
      effs_failed_unc_hi->Write();
      effs_failed_unc_lo->Write();
    }
  }
  ofile.Close();

  auto stop = high_resolution_clock::now();

  const int duration = duration_cast<seconds>(stop - start).count();
  cout << "INFO Finished processing in " << format_time(duration) << endl;
}
