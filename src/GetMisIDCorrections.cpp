// Author: Lucas Meyer Garcia
// License: BSD 2-clause
//
// Description: Calculate misid efficiencies for kaons and pions

#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <Math/Vector4D.h>
#include "TCanvas.h"
#include "TChain.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THStack.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"

#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooBifurGauss.h"
#include "RooBinning.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooHelpers.h"
#include "RooKeysPdf.h"  // Customized version in root-curated
#include "RooMinimizer.h"
#include "RooPlot.h"
#include "RooPowerLaw.h"  // Custom PDF defined in root-curated
#include "RooProdPdf.h"
#include "RooRealVar.h"

#include <yaml-cpp/yaml.h>
#include <cxxopts.hpp>

#include "utils.h"

using namespace RooFit;
using namespace std::chrono;

using std::cout, std::endl;
using std::string, std::unique_ptr, std::unordered_map, std::vector, std::array,
    std::pair;

///////////////
// Constants //
///////////////

constexpr int Dst_ID = 413;
constexpr int D0_ID  = 421;
constexpr int PI_ID  = 211;
constexpr int K_ID   = 321;
constexpr int MU_ID  = 13;
constexpr int KS_ID  = 310;

constexpr double DM_min   = 0.141;       // GeV
constexpr double DM_max   = 0.153;       // GeV
constexpr double DM       = 0.1454258;   // GeV
constexpr double D0_M_min = 1.825;       // GeV
constexpr double D0_M_max = 1.910;       // GeV
constexpr double D0_M     = 1.86484;     // GeV
constexpr double PI_M     = 0.13957039;  // GeV
constexpr double K_M      = 0.493677;    // GeV

constexpr double ONE_SIGMA = 0.682689492137;

//////////////////////////////////
// Generator level efficiencies //
//////////////////////////////////

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

// D0 -> K- pi+
constexpr array<double, 6> gen_effs_kpi = {0.005125, 0.005134, 0.005121,
                                           0.005119, 0.005110, 0.005132};
constexpr array<double, 6> gen_uncs_kpi = {0.000011, 0.000011, 0.000011,
                                           0.000011, 0.000011, 0.000011};

constexpr double gen_eff_kpi = weighted_average(gen_effs_kpi, gen_uncs_kpi);

// D0 -> pi+ pi- pi0
constexpr array<double, 4> gen_effs_pipipi0_dalitz = {0.18537, 0.18541, 0.18533,
                                                      0.18564};
constexpr array<double, 4> gen_uncs_pipipi0_dalitz = {0.00040, 0.00040, 0.00038,
                                                      0.00038};

constexpr double gen_eff_pipipi0_dalitz =
    weighted_average(gen_effs_pipipi0_dalitz, gen_uncs_pipipi0_dalitz);
constexpr double gen_eff_pipipi0_phsp =
    gen_eff_pipipi0_dalitz;  // Efficiencies not available in STATS09 table

// D0 -> pi- mu+ nu
constexpr array<double, 5> gen_effs_pimunu = {0.21296, 0.21242, 0.21262,
                                              0.21323, 0.21205};
constexpr array<double, 5> gen_uncs_pimunu = {0.00068, 0.00072, 0.00047,
                                              0.00046, 0.00049};

constexpr double gen_eff_pimunu =
    weighted_average(gen_effs_pimunu, gen_uncs_pimunu);

// D0 -> K- mu+ nu
// Ran gauss locally to produce gen level efficiencies.
// XML files are saved in the step1 directory in lhcb-ntuples-gen.
constexpr double kmunu_passed = 4989 + 4993 + 5060 + 4965 + 5040 + 5029 + 4998;
constexpr double kmunu_total =
    22432 + 22471 + 22413 + 22117 + 22285 + 22098 + 22275;

constexpr double gen_eff_kmunu = kmunu_passed / kmunu_total;
// Uncertainty ~0.00106

// D0 -> pi+ pi-
constexpr array<double, 7> gen_effs_pipi = {0.20392, 0.20364, 0.20443, 0.20518,
                                            0.20489, 0.20479, 0.20508};
constexpr array<double, 7> gen_uncs_pipi = {0.00063, 0.00061, 0.00043, 0.00044,
                                            0.00044, 0.00044, 0.00045};

constexpr double gen_eff_pipi = weighted_average(gen_effs_pipi, gen_uncs_pipi);

// D0 -> pi- e+ nu
constexpr array<double, 2> gen_effs_pienu = {0.20752, 0.20833};
constexpr array<double, 2> gen_uncs_pienu = {0.00071, 0.00070};

constexpr double gen_eff_pienu =
    weighted_average(gen_effs_pienu, gen_uncs_pienu);

// D0 -> K- e+ nu
constexpr array<double, 2> gen_effs_kenu = {0.21909, 0.21992};
constexpr array<double, 2> gen_uncs_kenu = {0.00070, 0.00075};

constexpr double gen_eff_kenu = weighted_average(gen_effs_kenu, gen_uncs_kenu);

// D0 -> Ks pi+ pi-
constexpr array<double, 2> gen_effs_kspipi = {0.2204, 0.2197};
constexpr array<double, 2> gen_uncs_kspipi = {0.0015, 0.0023};

constexpr double gen_eff_kspipi =
    weighted_average(gen_effs_kspipi, gen_uncs_kspipi);

//////////////////////////////
// Number of events in disk //
//////////////////////////////

// Signal
constexpr int N_disk_kpi =
    15462222 + 15903204 + 15001155 + 15044631 + 15045803 + 15017328;

// D0 bkg
constexpr int N_disk_pipipi0_phsp =
    1756054 + 1751894 + 1752298 + 1737653 + 502883 + 504796;
constexpr int N_disk_pipipi0_dalitz =
    1995881 + 2000074 + 2080042 + 2000189 + 1008307 + 1003547;
constexpr int N_disk_pimunu = 1776455 + 1758546 + 1750705 + 1888217;
constexpr int N_disk_kmunu  = 1751390 + 1757357 + 1750989 + 1740356 + 2023764 +
                             2005513 + 1047769 + 1035712;
constexpr int N_disk_pipi =
    1008347 + 1011491 + 2028005 + 2000032 + 2002621 + 2006539;
constexpr int N_disk_pienu  = 1756658 + 1751546 + 1750623 + 1979092;
constexpr int N_disk_kenu   = 1875792 + 1754493 + 1748598 + 1875587;
constexpr int N_disk_kspipi = 2231960 + 2283875;
constexpr int N_disk_klpipi = N_disk_kspipi;

////////////////////////////////
// Number of simulated events //
////////////////////////////////

constexpr double N_sim_kpi = N_disk_kpi / gen_eff_kpi;
constexpr double N_sim_pipipi0_phsp =
    N_disk_pipipi0_phsp / gen_eff_pipipi0_phsp;
constexpr double N_sim_pipipi0_dalitz =
    N_disk_pipipi0_dalitz / gen_eff_pipipi0_dalitz;
constexpr double N_sim_pimunu = N_disk_pimunu / gen_eff_pimunu;
constexpr double N_sim_kmunu  = N_disk_kmunu / gen_eff_kmunu;
constexpr double N_sim_pipi   = N_disk_pipi / gen_eff_pipi;
constexpr double N_sim_pienu  = N_disk_pienu / gen_eff_pienu;
constexpr double N_sim_kenu   = N_disk_pipi / gen_eff_kenu;
constexpr double N_sim_kspipi = N_disk_kspipi / gen_eff_kspipi;
constexpr double N_sim_klpipi = N_sim_kspipi;
// Because they are almost identical kinematically, we may use the kspipi sample
// as a klpipi sample as long as we ignore events in which one of the tracks was
// produced by a Ks decay product.

////////////////////////////////
// Number of simulated events //
////////////////////////////////

constexpr double BF_D0_kpi     = 0.03945;
constexpr double BF_D0_pipipi0 = 0.0149;
constexpr double BF_D0_pimunu  = 0.00267;
constexpr double BF_D0_kmunu   = 0.03418;
constexpr double BF_D0_pipi    = 0.001453;
constexpr double BF_D0_pienu   = 0.00291;
constexpr double BF_D0_kenu    = 0.03538;
constexpr double BF_D0_kspipi  = 0.0286;
constexpr double BF_D0_klpipi  = BF_D0_kspipi;
constexpr double BF_ks_pi0pi0  = 0.3069;
constexpr double BF_ks_pipi    = 0.6920;

/////////////
// Weights //
/////////////

// Normalize w.r.t. number of simulated pipipi0 DALITZ events
constexpr double w_pipipi0_phsp = (N_sim_pipipi0_dalitz / BF_D0_pipipi0) /
                                  (N_sim_pipipi0_phsp / BF_D0_pipipi0);
constexpr double w_pipipi0_dalitz = (N_sim_pipipi0_dalitz / BF_D0_pipipi0) /
                                    (N_sim_pipipi0_dalitz / BF_D0_pipipi0);
constexpr double w_pimunu =
    (N_sim_pipipi0_dalitz / BF_D0_pipipi0) / (N_sim_pimunu / BF_D0_pimunu);
constexpr double w_kmunu =
    (N_sim_pipipi0_dalitz / BF_D0_pipipi0) / (N_sim_kmunu / BF_D0_kmunu);
constexpr double w_pipi =
    (N_sim_pipipi0_dalitz / BF_D0_pipipi0) / (N_sim_pipi / BF_D0_pipi);
constexpr double w_pienu =
    (N_sim_pipipi0_dalitz / BF_D0_pipipi0) / (N_sim_pienu / BF_D0_pienu);
constexpr double w_kenu =
    (N_sim_pipipi0_dalitz / BF_D0_pipipi0) / (N_sim_kenu / BF_D0_kenu);
constexpr double w_kspipi_pipi = (N_sim_pipipi0_dalitz / BF_D0_pipipi0) /
                                 (N_sim_kspipi / (BF_D0_kspipi * BF_ks_pipi));
constexpr double w_kspipi_pi0pi0 = w_kspipi_pipi * (BF_ks_pi0pi0 / BF_ks_pipi);
constexpr double w_klpipi =
    (N_sim_pipipi0_dalitz / BF_D0_pipipi0) / (N_sim_kspipi / BF_D0_klpipi);

const unordered_map<string, double> w_d0_decays{
    {"pipipi0", w_pipipi0_dalitz},
    {"pimunu", w_pimunu},
    {"kmunu", w_kmunu},
    {"pipi", w_pipi},
    {"pienu", w_pienu},
    {"kenu", w_kenu},
    {"kspipi_pipi", w_kspipi_pipi},
    {"kspipi_pi0pi0", w_kspipi_pi0pi0},
    {"klpipi", w_kspipi_pipi}};

////////////
// Colors //
////////////

// Color pallete (https://jfly.uni-koeln.de/color/)
const TColor col1(0.00f, 0.00f, 0.00f);  // Black
const TColor col2(0.90f, 0.60f, 0.00f);  // Orange
const TColor col3(0.35f, 0.70f, 0.90f);  // Sky blue
const TColor col4(0.00f, 0.60f, 0.50f);  // Bluish green
const TColor col5(0.95f, 0.90f, 0.25f);  // Yellow
const TColor col6(0.00f, 0.45f, 0.70f);  // Blue
const TColor col7(0.80f, 0.40f, 0.00f);  // Vermillion
const TColor col8(0.80f, 0.60f, 0.70f);  // Reddish purple

const unordered_map<string, int> color_ids_d0_decays{
    {"pipipi0", col2.GetNumber()},     {"pimunu", col4.GetNumber()},
    {"kmunu", col6.GetNumber()},       {"pipi", col7.GetNumber()},
    {"pienu", col3.GetNumber()},       {"kenu", col5.GetNumber()},
    {"kspipi_pipi", col8.GetNumber()}, {"kspipi_pi0pi0", col8.GetNumber()},
    {"klpipi", col1.GetNumber()}};

/////////////
// Binning //
/////////////

constexpr double BINS_P[] = {3e+3, 6e+3, 10e+3, 15.6e+3, 27e+3, 60e+3, 100e+3};
constexpr double BINS_ETA[]     = {1.7, 3.6, 5.0};
constexpr double BINS_NTRACKS[] = {0., 200., 600.};

constexpr int N_BINS_P       = sizeof(BINS_P) / sizeof(double) - 1;
constexpr int N_BINS_ETA     = sizeof(BINS_ETA) / sizeof(double) - 1;
constexpr int N_BINS_NTRACKS = sizeof(BINS_NTRACKS) / sizeof(double) - 1;

const unordered_map<string, int> year_idx{
    {"2016", 0}, {"2017", 1}, {"2018", 2}};

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

// Add contents of <ds_source> to <ds_target>
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

//////////
// Main //
//////////

int main(int argc, char **argv) {
  auto start = high_resolution_clock::now();

  cxxopts::Options argOpts("GetMisIDCorrections",
                           "Calculate K/pi misid efficiencies.");

  // clang-format off
  argOpts.add_options()
    ("h,help", "Print help")
    ("d,debug", "Enable debug mode",
     cxxopts::value<bool>()->default_value("false"))
    ("p,particles", "Specify probed particle",
     cxxopts::value<vector<string>>()->default_value("k,pi"))
    ("c,config", "Specify input YAML config file",
     cxxopts::value<string>())
    ("b,bkgfile", "Specify input YAML with D0 bkg constraints",
     cxxopts::value<string>())
    ("y,year", "Specify data-taking year",
     cxxopts::value<string>()->default_value("2016"))
    ("o,output", "Specify output folder",
     cxxopts::value<string>()->default_value("gen/"))
    ("n,dry-run", "Do not perform fits",
     cxxopts::value<bool>()->default_value("false"))
    ("vmu", "Flag misid validation PID cuts",
     cxxopts::value<bool>()->default_value("false"))
    ("fake_mu", "Flag fake muon sample PID cuts",
     cxxopts::value<bool>()->default_value("false"))
    ("minos", "Use MINOS in fits to calibration samples",
     cxxopts::value<bool>()->default_value("false"))
    ("float_dif", "Float ratio of signal events with k/pi decay-in-flight",
     cxxopts::value<bool>()->default_value("true"))
    ;
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  // Read arguments
  const auto ymlFile   = parsedArgs["config"].as<string>();
  const auto d0BkgFile = parsedArgs["bkgfile"].as<string>();
  const auto particles = parsedArgs["particles"].as<vector<string>>();
  const auto dry_run   = parsedArgs["dry-run"].as<bool>();
  const auto vmu       = parsedArgs["vmu"].as<bool>();
  const auto fake_mu   = parsedArgs["fake_mu"].as<bool>();
  const auto use_minos = parsedArgs["minos"].as<bool>();
  const auto debug     = parsedArgs["debug"].as<bool>();
  const auto float_dif = parsedArgs["float_dif"].as<bool>();
  const auto year      = parsedArgs["year"].as<string>();

  const vector<string> years_mc = {"2016", "2017", "2018"};

  const TString opath = parsedArgs["output"].as<string>();

  if (vmu && fake_mu) {
    cout << "WARNING Both VMU and FAKE_MU flags set. The VMU flag is redundant "
            "as the FAKE MU flag implies no muBDT cut."
         << endl;
  } else if (vmu) {
    cout << "INFO Using VMU pid cuts" << endl;
  } else if (fake_mu) {
    cout << "INFO Using FAKE_MU pid cuts" << endl;
  }

  const string sample = fake_mu ? "fake_mu" : (vmu ? "vmu" : "iso_ctrl");

  if (dry_run) {
    cout << "WARNING Dry run. Fits will not be performed." << endl;
  }

  if (use_minos) cout << "INFO Using MINOS in calib fits" << endl;

  if (float_dif) {
    cout << "INFO Floating amount of ratio of DiF events in signal PDFs"
         << endl;
  }

  RooFit::MsgLevel msg_level = debug ? RooFit::INFO : RooFit::PROGRESS;
  RooHelpers::LocalChangeMsgLevel msg_level_plot(msg_level, 0u, Plotting,
                                                 false);

  constexpr int max_fix_reattempts = 10;

  // Parse YAML config
  const auto ymlConfig = YAML::LoadFile(ymlFile)["misid_corrections"];

  // Open YAML with D0 bkg constraints
  const auto ymlBkg = YAML::LoadFile(d0BkgFile);

  // Define histogram to easily determine kinematical bins
  TH3D histo_binning("histo_binning", ";p;#eta;nTracks", N_BINS_P, BINS_P,
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

  // Other variables
  TCanvas c_single("c_single", "c_single", 1280, 960);
  TCanvas c_double("c_double", "c_double", 640, 960);
  TCanvas c3("c3", "c3", 1920, 480);
  TCanvas c_four("c_four", "c_four", 1280, 960);
  c_four.Divide(2, 2);
  c3.Divide(3, 1);
  c_double.Divide(1, 2);

  // Define output file
  const TString opath_full = opath + "/mc_misid_histos.root";
  cout << "INFO Creating output file: " << opath_full << endl;
  TFile ofile(opath_full, "RECREATE");

  // Define histogram ranges
  // To allow shifting/scaling of histograms, need to feel bins immediately
  // below/above the fit range

  // Use different nbins but with matching bin boundaries
  // For now, same binning for d0_m and dm
  constexpr int nbins_d0_bkg = 15;
  constexpr int nbins_passed = 2 * nbins_d0_bkg;
  constexpr int nbins_failed = 2 * nbins_passed;

  constexpr int extra_bins_d0_bkg = 7;
  constexpr int extra_bins_passed = 2 * extra_bins_d0_bkg;
  constexpr int extra_bins_failed = 2 * extra_bins_passed;

  constexpr double r = static_cast<double>(nbins_d0_bkg) / extra_bins_d0_bkg;

  // Calculate histograsm ranges with additional bins
  // Limits are the same for all samples
  constexpr double extended_d0_m_max = D0_M_max + (D0_M_max - D0_M_min) / r;
  constexpr double extended_d0_m_min = D0_M_min - (D0_M_max - D0_M_min) / r;
  constexpr double extended_dm_max   = DM_max + (DM_max - DM_min) / r;
  constexpr double extended_dm_min   = DM_min - (DM_max - DM_min) / r;

  // Calculate number of bins for RooKeysPdf. The ideai is to have ~2k bins in
  // the fit range, scale the count porperly to account for the additional
  // range.
  constexpr int n_bins_keys =
      (2000 * (nbins_d0_bkg + 2 * extra_bins_d0_bkg)) / nbins_d0_bkg;

  // Observables
  RooRealVar d0_m_var("d0_m_var", "d0_m_var", extended_d0_m_min,
                      extended_d0_m_max, "GeV");
  RooRealVar dm_var("dm_var", "dm_var", extended_dm_min, extended_dm_max,
                    "GeV");

  d0_m_var.setRange("data", D0_M_min, D0_M_max);
  dm_var.setRange("data", DM_min, DM_max);

  RooRealVar d0_m_shift("d0_m_shift", "d0_m_shift", 0., -0.003, 0.003, "GeV");
  RooRealVar dm_shift("dm_shift", "dm_shift", 0., -0.0002, 0.0002, "GeV");
  RooRealVar d0_m_scale("d0_m_scale", "d0_m_scale", 1., 0.9, 1.1, "");
  RooRealVar dm_scale("dm_scale", "dm_scale", 1., 0.9, 1.1, "");

  RooConstVar d0_m_pdg("d0_m_pdg", "d0_m_pdg", D0_M);  // GeV
  RooConstVar dm_pdg("dm_pdg", "dm_pdg", DM);          // GeV

  // Binning schemes with additional bins to build templates that can be shifted
  RooBinning bins_histos_d0_m_failed(nbins_failed + 2 * extra_bins_failed,
                                     extended_d0_m_min, extended_d0_m_max,
                                     "bins_histos_d0_m_failed");
  RooBinning bins_histos_dm_failed(nbins_failed + 2 * extra_bins_failed,
                                   extended_dm_min, extended_dm_max,
                                   "bins_histos_dm_failed");
  RooBinning bins_histos_d0_m_passed(nbins_passed + 2 * extra_bins_passed,
                                     extended_d0_m_min, extended_d0_m_max,
                                     "bins_histos_d0_m_passed");
  RooBinning bins_histos_dm_passed(nbins_passed + 2 * extra_bins_passed,
                                   extended_dm_min, extended_dm_max,
                                   "bins_histos_dm_passed");
  RooBinning bins_histos_d0_m_d0_bkg(nbins_d0_bkg + 2 * extra_bins_d0_bkg,
                                     extended_d0_m_min, extended_d0_m_max,
                                     "bins_histos_d0_m_d0_bkg");
  RooBinning bins_histos_dm_d0_bkg(nbins_d0_bkg + 2 * extra_bins_d0_bkg,
                                   extended_dm_min, extended_dm_max,
                                   "bins_histos_dm_d0_bkg");

  d0_m_var.setBinning(bins_histos_d0_m_failed, "bins_histos_d0_m_failed");
  d0_m_var.setBinning(bins_histos_d0_m_passed, "bins_histos_d0_m_passed");
  dm_var.setBinning(bins_histos_dm_failed, "bins_histos_dm_failed");
  dm_var.setBinning(bins_histos_dm_passed, "bins_histos_dm_passed");

  /////////////////////////////
  // Static fit components   //
  /////////////////////////////

  // dm comb model: threshold function
  RooConstVar pi_m_pdg("pi_m_pdg", "pi_m_pdg", PI_M);  // GeV

  RooRealVar c_comb_spi_failed("c_comb_spi_failed", "c_comb_spi_failed", 0., 1.,
                               "");
  RooRealVar c_comb_spi_passed("c_comb_spi_passed", "c_comb_spi_passed", 0., 1.,
                               "");

  RooRealVar c_comb_all_failed("c_comb_all_failed", "c_comb_all_failed", 0., 1.,
                               "");
  RooRealVar c_comb_all_passed("c_comb_all_passed", "c_comb_all_passed", 0., 1.,
                               "");

  RooPowerLaw dm_comb_spi_failed("dm_comb_spi_failed", "dm_comb_spi_failed",
                                 dm_var, pi_m_pdg, c_comb_spi_failed);
  RooPowerLaw dm_comb_spi_passed("dm_comb_spi_passed", "dm_comb_spi_passed",
                                 dm_var, pi_m_pdg, c_comb_spi_passed);

  RooPowerLaw dm_comb_all_failed("dm_comb_all_failed", "dm_comb_all_failed",
                                 dm_var, pi_m_pdg, c_comb_all_failed);
  RooPowerLaw dm_comb_all_passed("dm_comb_all_passed", "dm_comb_all_passed",
                                 dm_var, pi_m_pdg, c_comb_all_passed);

  // D0 comb model: exponential distribution
  RooRealVar k_comb_all_passed("k_comb_all_passed", "k_comb_all_passed", -20.,
                               0., "");
  RooRealVar k_comb_all_failed("k_comb_all_failed", "k_comb_all_failed", -20.,
                               0., "");

  RooExponential d0_comb_all_failed("d0_comb_all_failed", "d0_comb_all_failed",
                                    d0_m_var, k_comb_all_failed);
  RooExponential d0_comb_all_passed("d0_comb_all_passed", "d0_comb_all_passed",
                                    d0_m_var, k_comb_all_passed);

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

  // Scale controlling yield of mis-identified non-dif tracks in "passed"
  // sample. The nominal MC prediction is recovered when the scale is 1.
  RooRealVar scale_nondif_calib("scale_nondif_calib", "scale_nondif_calib", 1.,
                                0.1, 20., "");

  // Number of non-dif events shifted from failed to passed
  RooFormulaVar nondif_yield_add_passed(
      "nondif_yield_add_passed", "nondif_yield_add_passed",
      "(x[0] - 1.0) * x[1]",
      RooArgList(scale_nondif_calib, n_inmw_passed_nondif));

  RooFormulaVar f_dif_calib_passed(
      "f_dif_calib_passed", "f_dif_calib_passed", "x[0]/(x[0]+x[1]+x[2])",
      RooArgList(n_inmw_passed_dif, n_inmw_passed_nondif,
                 nondif_yield_add_passed));
  RooFormulaVar f_dif_calib_failed(
      "f_dif_calib_failed", "f_dif_calib_failed", "x[0]/(x[0]+x[1]-x[2])",
      RooArgList(n_inmw_failed_dif, n_inmw_failed_nondif,
                 nondif_yield_add_passed));

  // Built heuristic constraint for non-dif scale
  RooBifurGauss scale_nondif_constrain(
      "scale_nondif_constrain", "scale_nondif_constrain", scale_nondif_calib,
      RooConst(1.0), RooConst(0.1), RooConst(5. / 3.));

  ///////////////////
  // Fit fractions //
  ///////////////////

  // Fraction of "physics" events i.e. signal + D0 bkg
  // Consequently, fraction of combinatorial bkg = 1 - f_phys
  RooRealVar f_phys_passed("f_phys_passed", "f_phys_passed", 0., 1., "");
  RooRealVar f_phys_failed("f_phys_failed", "f_phys_failed", 0., 1., "");

  // Fraction of signal in "physics" events
  RooRealVar f_sig_passed("f_sig_passed", "f_sig_passed", 0., 1., "");
  RooRealVar f_sig_failed("f_sig_failed", "f_sig_failed", 0., 1., "");

  // Fraction of soft-pion bkg in total combinatorial bkg
  RooRealVar f_spi_passed("f_spi_passed", "f_spi_passed", 0., 1., "");
  RooRealVar f_spi_failed("f_spi_failed", "f_spi_failed", 0., 1., "");

  // Signal
  RooFormulaVar f_sig_phys_passed("f_sig_phys_passed", "f_sig_phys_passed",
                                  "x[0] * x[1]",
                                  RooArgList(f_phys_passed, f_sig_passed));
  RooFormulaVar f_sig_phys_failed("f_sig_phys_failed", "f_sig_phys_failed",
                                  "x[0] * x[1]",
                                  RooArgList(f_phys_failed, f_sig_failed));

  // Mis-reconstructed D0 bkg
  RooFormulaVar f_d0_bkg_passed("f_d0_bkg_passed", "f_d0_bkg_passed",
                                "x[0] * (1.0 - x[1])",
                                RooArgList(f_phys_passed, f_sig_passed));
  RooFormulaVar f_d0_bkg_failed("f_d0_bkg_failed", "f_d0_bkg_failed",
                                "x[0] * (1.0 - x[1])",
                                RooArgList(f_phys_failed, f_sig_failed));

  // Soft-pion combinatorial bkg
  RooFormulaVar f_comb_spi_passed("f_comb_spi_passed", "f_comb_spi_passed",
                                  "(1.0 - x[0]) * x[1]",
                                  RooArgList(f_phys_passed, f_spi_passed));
  RooFormulaVar f_comb_spi_failed("f_comb_spi_failed", "f_comb_spi_failed",
                                  "(1.0 - x[0]) * x[1]",
                                  RooArgList(f_phys_failed, f_spi_failed));

  // Complete combinatorial bkg
  RooFormulaVar f_comb_all_passed("f_comb_all_passed", "f_comb_all_passed",
                                  "(1.0 - x[0]) * (1.0 - x[1])",
                                  RooArgList(f_phys_passed, f_spi_passed));
  RooFormulaVar f_comb_all_failed("f_comb_all_failed", "f_comb_all_failed",
                                  "(1.0 - x[0]) * (1.0 - x[1])",
                                  RooArgList(f_phys_failed, f_spi_failed));

  ////////////////////
  // Normalizations //
  ////////////////////

  // Total
  RooRealVar n_passed("n_passed", "n_passed", 0., "");
  RooRealVar n_failed("n_failed", "n_failed", 0., "");

  // Signal
  RooFormulaVar n_sig_passed("n_sig_passed", "n_sig_passed", "x[0] * x[1]",
                             RooArgList(f_sig_phys_passed, n_passed));
  RooFormulaVar n_sig_failed("n_sig_failed", "n_sig_failed", "x[0] * x[1]",
                             RooArgList(f_sig_phys_failed, n_failed));

  // Constraints for amount of D0 bkg
  RooRealVar f_sig_passed_mean("f_sig_passed_mean", "f_sig_passed_mean", 1.,
                               "");
  RooRealVar f_sig_passed_unc_hi("f_sig_passed_unc_hi", "f_sig_passed_unc_hi",
                                 0., "");
  RooRealVar f_sig_passed_unc_lo("f_sig_passed_unc_lo", "f_sig_passed_unc_lo",
                                 0., "");

  RooBifurGauss f_sig_passed_constrain(
      "f_sig_passed_constrain", "f_sig_passed_constrain", f_sig_passed,
      f_sig_passed_mean, f_sig_passed_unc_lo, f_sig_passed_unc_hi);

  RooRealVar f_sig_failed_mean("f_sig_failed_mean", "f_sig_failed_mean", 1.,
                               "");
  RooRealVar f_sig_failed_unc_hi("f_sig_failed_unc_hi", "f_sig_failed_unc_hi",
                                 0., "");
  RooRealVar f_sig_failed_unc_lo("f_sig_failed_unc_lo", "f_sig_failed_unc_lo",
                                 0., "");

  RooBifurGauss f_sig_failed_constrain(
      "f_sig_failed_constrain", "f_sig_failed_constrain", f_sig_failed,
      f_sig_failed_mean, f_sig_failed_unc_lo, f_sig_failed_unc_hi);

  // PID efficiency
  RooFormulaVar eff("eff", "eff", "x[0] / (x[0] + x[1])",
                    RooArgList(n_sig_passed, n_sig_failed));

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

  auto set_parameters_calib = [&](const string &probe) {
    scale_nondif_calib.setVal(1.0);
    d0_m_shift.setVal(0.);
    dm_shift.setVal(0.);
    d0_m_scale.setVal(1.);
    dm_scale.setVal(1.);
    if (probe == "k") {
      c_comb_spi_failed.setVal(0.46);
      c_comb_spi_passed.setVal(0.46);
      c_comb_all_failed.setVal(0.41);
      c_comb_all_passed.setVal(0.41);
    } else if (probe == "pi") {
      c_comb_spi_failed.setVal(0.45);
      c_comb_spi_passed.setVal(0.45);
      c_comb_all_failed.setVal(0.44);
      c_comb_all_passed.setVal(0.37);
    }
  };

  // RooArgSet with fit observables (d0_m_var, dm_var)
  RooArgSet fit_vars(d0_m_var, dm_var);

  RooRealVar w_d0_bkg("w_d0_bkg", "w_d0_bkg", 1., "");

  // RooArgSet with fit observables and weight (d0_m_var, dm_var, w_d0_bkg)
  RooArgSet fit_vars_w(d0_m_var, dm_var, w_d0_bkg);

  RooArgSet argset_minos(eff, scale_nondif_calib);

  // Monitor how much f_sig deviates from constraint
  RooRealVar f_sig_passed_pull("f_sig_passed_pull", "f_sig_passed_pull", 0.,
                               -6., 6., "");
  RooRealVar f_sig_failed_pull("f_sig_failed_pull", "f_sig_failed_pull", 0.,
                               -6., 6., "");

  RooArgSet params_calib(
      c_comb_spi_failed, c_comb_spi_passed, c_comb_all_failed,
      c_comb_all_passed, k_comb_all_failed, k_comb_all_passed, f_phys_passed,
      f_phys_failed, f_sig_passed, f_sig_failed, f_spi_passed, f_spi_failed,
      scale_nondif_calib, f_sig_passed_pull, f_sig_failed_pull);

  // Check proportions in mis-reconstructed D0 bkg
  unordered_map<string, vector<unordered_map<string, int>>>
      d0_bkg_counts_passed{{"k", vector<unordered_map<string, int>>(6)},
                           {"pi", vector<unordered_map<string, int>>(6)}};
  unordered_map<string, vector<unordered_map<string, int>>>
      d0_bkg_counts_failed{{"k", vector<unordered_map<string, int>>(6)},
                           {"pi", vector<unordered_map<string, int>>(6)}};

  for (auto probe : particles) {
    cout << "INFO Selecting " << probe << endl;

    const TString tag = (probe == "pi") ? "k" : "pi";

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
    cout << "INFO Initializing signal MC datahists" << endl;
    for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
      for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
        TString suffix =
            TString::Format("%s_%d_%d", probe.c_str(), eta_idx, p_idx);
        datasets_mc_passed_dif[eta_idx][p_idx] =
            new RooDataSet("ds_mc_passed_dif_" + suffix,
                           "ds_mc_passed_dif_" + suffix, fit_vars);
        datasets_mc_failed_dif[eta_idx][p_idx] =
            new RooDataSet("ds_mc_failed_dif_" + suffix,
                           "ds_mc_failed_dif_" + suffix, fit_vars);
        for (int ntrks_idx = 0; ntrks_idx < N_BINS_NTRACKS; ntrks_idx++) {
          suffix.Form("%s_%d_%d_%d", probe.c_str(), ntrks_idx, eta_idx, p_idx);
          datasets_mc_passed_nondif[ntrks_idx][eta_idx][p_idx] =
              new RooDataSet("ds_mc_passed_nondif_" + suffix,
                             "ds_mc_passed_nondif_" + suffix, fit_vars);
          datasets_mc_failed_nondif[ntrks_idx][eta_idx][p_idx] =
              new RooDataSet("ds_mc_failed_nondif_" + suffix,
                             "ds_mc_failed_nondif_" + suffix, fit_vars);
        }
      }
    }

    // Loop over MC files first to build combined MC sample,
    // but calculating separate mass window efficiencies
    for (auto year : years_mc) {
      // Open and loop over MC files
      const auto mc_paths =
          ymlConfig["mc_ntps"]["signal"][year][probe].as<vector<string>>();
      TChain ch_mc("tree");
      for (const auto &mc_path : mc_paths) {
        cout << "INFO Opening " << year << " signal MC files: " << mc_path
             << endl;
        ch_mc.Add(mc_path.c_str());
      }
      // Alex's sample yeilds +10M events out of 60M
      // (+ 16%)
      cout << "INFO Opened MC files:" << endl;
      print_files(ch_mc);

      // Define variable to access input ntuples
      int dst_id, probe_trueid, probe_daughter0_trueid, ntracks;

      double dst_m, d0_m, probe_p, probe_pz, probe_pt, probe_dllmu, probe_dlle,
          spi_p, spi_pt, tag_p, tag_pz, tag_pt, k_track_chi2ndof,
          pi_track_chi2ndof, spi_track_chi2ndof, k_px, k_py, pi_px, pi_py;

      float probe_mu_ubdt;

      bool probe_ismuon, probe_hasmuon;

      cout << "INFO Setting input branches " << endl;
      ch_mc.SetBranchStatus("*", false);
      ch_mc.SetBranchAddress("dst_M", &dst_m);
      ch_mc.SetBranchAddress("dst_TRUEID", &dst_id);
      ch_mc.SetBranchAddress("d0_M", &d0_m);
      ch_mc.SetBranchAddress("spi_P", &spi_p);
      ch_mc.SetBranchAddress("spi_PT", &spi_pt);
      ch_mc.SetBranchAddress(tag + "_P", &tag_p);
      ch_mc.SetBranchAddress(tag + "_PT", &tag_pt);
      ch_mc.SetBranchAddress(tag + "_PZ", &tag_pz);
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
      ch_mc.SetBranchAddress("k_PX", &k_px);
      ch_mc.SetBranchAddress("k_PY", &k_py);
      ch_mc.SetBranchAddress("pi_PX", &pi_px);
      ch_mc.SetBranchAddress("pi_PY", &pi_py);
      ch_mc.SetBranchAddress("k_TRACK_CHI2NDOF", &k_track_chi2ndof);
      ch_mc.SetBranchAddress("pi_TRACK_CHI2NDOF", &pi_track_chi2ndof);
      ch_mc.SetBranchAddress("spi_TRACK_CHI2NDOF", &spi_track_chi2ndof);
      ch_mc.SetBranchAddress("nTracks", &ntracks);

      int count_tm = 0, count_cond = 0, count_calib_sel = 0, count_range = 0,
          count_mw = 0;

      const int entries_mc = ch_mc.GetEntries();
      cout << "INFO Starting MC event loop over " << entries_mc << " entries"
           << endl;

      int eta_bin = 0, p_bin = 0, ntrks_bin = 0;

      auto start_mc_loop = high_resolution_clock::now();
      for (int evt = 0; evt < entries_mc; evt++) {
        ch_mc.GetEntry(evt);

        // Offline truth-matching
        if (std::abs(dst_id) != Dst_ID) continue;
        count_tm++;

        // Conditional cuts
        if (!probe_hasmuon) continue;
        count_cond++;

        // Calib sample cuts
        // See tables 36 and 48 in LHCb-PUB-2016-005
        // See also
        // https://lhcb-pid-wgp-plots.web.cern.ch/lhcb-pid-wgp-plots/Run2/ and
        // https://gitlab.cern.ch/lhcb-datapkg/WG/SemilepConfig/-/blob/master/options/Filter_Dstar2D02KpiNoPID_2016MC.py?ref_type=heads
        if ((tag_p < 2000.) || (tag_pt < 250.)) continue;
        if ((probe_p < 2000.) || (probe_pt < 250.)) continue;
        if ((spi_p < 1000.) || (spi_pt < 100.)) continue;
        if ((k_track_chi2ndof > 3.) || (pi_track_chi2ndof > 3.) ||
            (spi_track_chi2ndof > 3.))
          continue;

        double k_pz, pi_pz;
        if (probe == "k") {
          k_pz  = probe_pz;
          pi_pz = tag_pz;
        } else {
          k_pz  = tag_pz;
          pi_pz = probe_pz;
        }

        ROOT::Math::PxPyPzMVector p_k(k_px * 0.001, k_py * 0.001, k_pz * 0.001,
                                      K_M);
        ROOT::Math::PxPyPzMVector p_pi(pi_px * 0.001, pi_py * 0.001,
                                       pi_pz * 0.001, PI_M);

        if ((std::max(probe_pt, tag_pt) < 1000.) || ((p_k + p_pi).Pt() < 1.5))
          continue;

        // Veto wrong-mass hypotheses
        // pi pi
        p_k.SetM(PI_M);
        p_pi.SetM(PI_M);
        if (std::abs((p_k + p_pi).M() - D0_M) < 0.025) continue;

        // pi K (wrong sign)
        p_k.SetM(PI_M);
        p_pi.SetM(K_M);
        if (std::abs((p_k + p_pi).M() - D0_M) < 0.025) continue;

        // K K
        p_k.SetM(K_M);
        p_pi.SetM(K_M);
        if (std::abs((p_k + p_pi).M() - D0_M) < 0.025) continue;

        count_calib_sel++;

        // MuonUnbiased equivalent to L0+HLT1 TIS. See MuonUnbiased defition at
        // https://gitlab.cern.ch/lhcb-datapkg/WG/PIDCalib/-/blob/master/scriptsR2/makeTuples_pp_2016_reprocessing.py#L71
        // and https://mattermost.web.cern.ch/lhcb/pl/893yre484jggigooti5u3gqb8w

        // Not applying MuonUnbiased and GHOSTPROB conditional cuts in MC since
        // some kinematical bins have very low statistics

        if (probe_p < 3000. || probe_p >= 100000. || ntracks >= 600) continue;

        const double probe_eta =
            0.5 * log((probe_p + probe_pz) / (probe_p - probe_pz));
        if (probe_eta < 1.7 || probe_eta >= 5.0) continue;
        count_range++;

        // Determine kinematical bin
        const int kin_bin = histo_binning.FindBin(probe_p, probe_eta, ntracks);
        histo_binning.GetBinXYZ(kin_bin, p_bin, eta_bin, ntrks_bin);

        // Fill histograms
        const double dm = (dst_m - d0_m) * 0.001;
        d0_m            = d0_m * 0.001;

        const bool reduced_fit_range = (probe == "pi") && (p_bin <= 2);

        const bool in_extended_window =
            in_var_range(d0_m_var, d0_m) && in_var_range(dm_var, dm);

        const bool in_fit_window = reduced_fit_range
                                       ? in_range(D0_M_min, d0_m, 1.900) &&
                                             in_range(DM_min, dm, DM_max)
                                       : in_range(D0_M_min, d0_m, D0_M_max) &&
                                             in_range(DM_min, dm, DM_max);
        if (in_fit_window) count_mw++;

        // Check for decay in flight of probe hadron
        const bool dif = (std::abs(probe_trueid) == MU_ID) ||
                         (std::abs(probe_daughter0_trueid) == MU_ID);

        // PID
        bool pid_ok;
        if (fake_mu) {
          // For FAKE_MU, we calculate the complementary efficiency so that
          // the "passed" sample always corresponds to the K/pi misid case
          pid_ok = probe_ismuon;
        } else if (vmu) {
          pid_ok = probe_ismuon && probe_dllmu > 2.0 && probe_dlle < 1.0 &&
                   probe_mu_ubdt < 0.25;
        } else {
          pid_ok = probe_ismuon && probe_dllmu > 2.0 && probe_dlle < 1.0 &&
                   probe_mu_ubdt > 0.25;
        }

        if (pid_ok) {
          // Fill "passed" samples
          if (dif) {
            count_total_passed_dif[year_idx.at(year)][ntrks_bin - 1]
                                  [eta_bin - 1][p_bin - 1]++;
            if (in_extended_window) {
              d0_m_var.setVal(d0_m);
              dm_var.setVal(dm);
              datasets_mc_passed_dif[eta_bin - 1][p_bin - 1]->addFast(fit_vars);
              if (in_fit_window) {
                count_in_mw_passed_dif[year_idx.at(year)][ntrks_bin - 1]
                                      [eta_bin - 1][p_bin - 1]++;
              }
            }
          } else {
            count_total_passed_nondif[year_idx.at(year)][ntrks_bin - 1]
                                     [eta_bin - 1][p_bin - 1]++;
            if (in_extended_window) {
              d0_m_var.setVal(d0_m);
              dm_var.setVal(dm);
              datasets_mc_passed_nondif[ntrks_bin - 1][eta_bin - 1][p_bin - 1]
                  ->addFast(fit_vars);
              if (in_fit_window) {
                count_in_mw_passed_nondif[year_idx.at(year)][ntrks_bin - 1]
                                         [eta_bin - 1][p_bin - 1]++;
              }
            }
          }
        } else {
          // Fill "failed" samples
          if (dif) {
            count_total_failed_dif[year_idx.at(year)][ntrks_bin - 1]
                                  [eta_bin - 1][p_bin - 1]++;
            if (in_extended_window) {
              d0_m_var.setVal(d0_m);
              dm_var.setVal(dm);
              datasets_mc_failed_dif[eta_bin - 1][p_bin - 1]->addFast(fit_vars);
              if (in_fit_window) {
                count_in_mw_failed_dif[year_idx.at(year)][ntrks_bin - 1]
                                      [eta_bin - 1][p_bin - 1]++;
              }
            }
          } else {
            count_total_failed_nondif[year_idx.at(year)][ntrks_bin - 1]
                                     [eta_bin - 1][p_bin - 1]++;
            if (in_extended_window) {
              d0_m_var.setVal(d0_m);
              dm_var.setVal(dm);
              datasets_mc_failed_nondif[ntrks_bin - 1][eta_bin - 1][p_bin - 1]
                  ->addFast(fit_vars);
              if (in_fit_window) {
                count_in_mw_failed_nondif[year_idx.at(year)][ntrks_bin - 1]
                                         [eta_bin - 1][p_bin - 1]++;
              }
            }
          }
        }
      }

      auto stop_mc_loop = high_resolution_clock::now();

      const int duration_mc_loop =
          duration_cast<seconds>(stop_mc_loop - start_mc_loop).count();
      cout << "INFO MC loop took " << format_time(duration_mc_loop) << endl;

      cout << "INFO MC cutflow:\n";
      cout << " - Total MC events:              " << entries_mc << "\n";
      cout << " - After offline truth-matching: " << count_tm << "("
           << count_tm * 100. / entries_mc << "%)\n";
      cout << " - After conditional cuts:       " << count_cond << "("
           << count_cond * 100. / count_tm << "%)\n";
      cout << " - After PIDCalib cuts:          " << count_calib_sel << "("
           << count_calib_sel * 100. / count_cond << "%)\n";
      cout << " - After range cuts:             " << count_range << "("
           << count_range * 100. / count_calib_sel << "%)\n";
      cout << " - After mass window cuts:       " << count_mw << "("
           << count_mw * 100. / count_range << "%)\n";
      cout << " Overall: " << "(" << count_mw * 100. / entries_mc << "%)\n"
           << endl;
    }

    // Deduce list of background sources from list of normalization scale
    // factors
    vector<string> d0_bkg_decays;
    d0_bkg_decays.reserve(w_d0_decays.size());
    for (const auto &kv : w_d0_decays) {
      d0_bkg_decays.push_back(kv.first);
    }

    // Datasets to build PDFs for D0 bkg
    unordered_map<string, vector<RooDataSet>> datasets_d0_bkg_mc_passed;
    unordered_map<string, vector<RooDataSet>> datasets_d0_bkg_mc_failed;

    // Initialize MC datasets
    cout << "INFO Initializing D0 background MC datahists" << endl;
    datasets_d0_bkg_mc_passed.reserve(d0_bkg_decays.size());
    datasets_d0_bkg_mc_failed.reserve(d0_bkg_decays.size());

    for (const auto &d0_decay : d0_bkg_decays) {
      datasets_d0_bkg_mc_passed[d0_decay].reserve(N_BINS_P);
      datasets_d0_bkg_mc_failed[d0_decay].reserve(N_BINS_P);

      for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
        TString suffix =
            TString::Format("%s_%s_%d", d0_decay.c_str(), probe.c_str(), p_idx);

        datasets_d0_bkg_mc_passed[d0_decay].emplace_back(
            RooDataSet("ds_d0_bkg_mc_passed_" + suffix,
                       "ds_d0_bkg_mc_passed_" + suffix, fit_vars));
        datasets_d0_bkg_mc_failed[d0_decay].emplace_back(
            RooDataSet("ds_d0_bkg_mc_failed_" + suffix,
                       "ds_d0_bkg_mc_failed_" + suffix, fit_vars));
      }
    }

    // Loop over MC files to build templates
    for (string d0_decay : d0_bkg_decays) {
      // To make my life easier, I create separate datasets for the D0 ->
      // (ks->pi+pi-)pi+pi- and for the emulated D0 -> (ks->pi0pi0)pi+pi- and D0
      // -> (kl->X)pi+pi- decays, to be combine with porper weights afterwards.
      // This requires two labels kspipi_pipi and kspipi_pi0pi0. This string
      // should evaluate to kspipi in those three cases, to pickup the correct
      // ntuple.
      TString d0_decay_ntuple = d0_decay;
      if (d0_decay_ntuple.Contains("kspipi") ||
          d0_decay_ntuple.Contains("klpipi")) {
        d0_decay_ntuple = "kspipi";
      }
      const bool veto_ks_daughters =
          (d0_decay == "klpipi") || (d0_decay == "kspipi_pi0pi0");

      // Open and loop over MC files
      string mc_path;
      if (ymlConfig["mc_ntps"][d0_decay_ntuple.Data()]) {
        mc_path =
            ymlConfig["mc_ntps"][d0_decay_ntuple.Data()][probe].as<string>();
      } else {
        // Some samples are not available for all years
        cout << "WARNING Could not find input path for " << d0_decay
             << " sample. Skipping." << endl;
        continue;
      }
      cout << "INFO Opening " << d0_decay << " MC files: " << mc_path << endl;
      TChain ch_mc("tree");
      ch_mc.Add(mc_path.c_str());
      cout << "INFO Opened MC files:" << endl;
      print_files(ch_mc);

      // Define variable to access input ntuples
      int dst_id, probe_trueid, ntracks, k_mother_id, k_gd_mother_id,
          k_gd_gd_mother_id, k_gd_gd_gd_mother_id, pi_mother_id,
          pi_gd_mother_id, pi_gd_gd_mother_id, pi_gd_gd_gd_mother_id;

      double dst_m, d0_m, probe_p, probe_pz, probe_pt, probe_dllmu, probe_dlle,
          spi_p, spi_pt, tag_p, tag_pz, tag_pt, k_track_chi2ndof,
          pi_track_chi2ndof, spi_track_chi2ndof, k_px, k_py, pi_px, pi_py;

      float probe_mu_ubdt;

      bool probe_ismuon, probe_hasmuon;

      cout << "INFO Setting input branches " << endl;
      ch_mc.SetBranchStatus("*", false);
      ch_mc.SetBranchAddress("dst_M", &dst_m);
      ch_mc.SetBranchAddress("dst_TRUEID", &dst_id);
      ch_mc.SetBranchAddress("d0_M", &d0_m);
      ch_mc.SetBranchAddress("spi_P", &spi_p);
      ch_mc.SetBranchAddress("spi_PT", &spi_pt);
      ch_mc.SetBranchAddress(tag + "_P", &tag_p);
      ch_mc.SetBranchAddress(tag + "_PT", &tag_pt);
      ch_mc.SetBranchAddress(tag + "_PZ", &tag_pz);
      ch_mc.SetBranchAddress((probe + "_P").c_str(), &probe_p);
      ch_mc.SetBranchAddress((probe + "_PT").c_str(), &probe_pt);
      ch_mc.SetBranchAddress((probe + "_PZ").c_str(), &probe_pz);
      ch_mc.SetBranchAddress((probe + "_isMuon").c_str(), &probe_ismuon);
      ch_mc.SetBranchAddress((probe + "_hasMuon").c_str(), &probe_hasmuon);
      ch_mc.SetBranchAddress((probe + "_PIDmu").c_str(), &probe_dllmu);
      ch_mc.SetBranchAddress((probe + "_PIDe").c_str(), &probe_dlle);
      ch_mc.SetBranchAddress((probe + "_bdt_mu").c_str(), &probe_mu_ubdt);
      ch_mc.SetBranchAddress((probe + "_TRUEID").c_str(), &probe_trueid);
      ch_mc.SetBranchAddress("k_PX", &k_px);
      ch_mc.SetBranchAddress("k_PY", &k_py);
      ch_mc.SetBranchAddress("k_MC_MOTHER_ID", &k_mother_id);
      ch_mc.SetBranchAddress("k_MC_GD_MOTHER_ID", &k_gd_mother_id);
      ch_mc.SetBranchAddress("k_MC_GD_GD_MOTHER_ID", &k_gd_gd_mother_id);
      ch_mc.SetBranchAddress("k_MC_GD_GD_GD_MOTHER_ID", &k_gd_gd_gd_mother_id);
      ch_mc.SetBranchAddress("pi_PX", &pi_px);
      ch_mc.SetBranchAddress("pi_PY", &pi_py);
      ch_mc.SetBranchAddress("pi_MC_MOTHER_ID", &pi_mother_id);
      ch_mc.SetBranchAddress("pi_MC_GD_MOTHER_ID", &pi_gd_mother_id);
      ch_mc.SetBranchAddress("pi_MC_GD_GD_MOTHER_ID", &pi_gd_gd_mother_id);
      ch_mc.SetBranchAddress("pi_MC_GD_GD_GD_MOTHER_ID",
                             &pi_gd_gd_gd_mother_id);
      ch_mc.SetBranchAddress("k_TRACK_CHI2NDOF", &k_track_chi2ndof);
      ch_mc.SetBranchAddress("pi_TRACK_CHI2NDOF", &pi_track_chi2ndof);
      ch_mc.SetBranchAddress("spi_TRACK_CHI2NDOF", &spi_track_chi2ndof);
      ch_mc.SetBranchAddress("nTracks", &ntracks);

      int count_tm = 0, count_cond = 0, count_calib_sel = 0, count_range = 0,
          count_mw = 0;

      int eta_bin = 0, p_bin = 0, ntrks_bin = 0;

      w_d0_bkg.setVal(w_d0_decays.at(d0_decay));

      const int entries_mc = ch_mc.GetEntries();
      cout << "INFO Starting MC event loop over " << entries_mc << " entries"
           << endl;

      auto start_mc_loop = high_resolution_clock::now();
      for (int evt = 0; evt < entries_mc; evt++) {
        ch_mc.GetEntry(evt);

        // Offline truth-matching
        if (std::abs(dst_id) != Dst_ID) continue;

        // For emulated D0 -> (ks->pi0pi0)pi+pi- and D0 -> (kl->X)pi+pi-
        // decays only, veto events where K or pi candidates are Ks decay
        // products.
        // We certainty don't need to check four generations, but why not
        // ¯\_(ツ)_/¯
        const bool from_ks = veto_ks_daughters
                                 ? (abs(k_mother_id) == KS_ID) ||
                                       (abs(pi_mother_id) == KS_ID) ||
                                       (abs(k_gd_mother_id) == KS_ID) ||
                                       (abs(pi_gd_mother_id) == KS_ID) ||
                                       (abs(k_gd_gd_mother_id) == KS_ID) ||
                                       (abs(pi_gd_gd_mother_id) == KS_ID) ||
                                       (abs(k_gd_gd_gd_mother_id) == KS_ID) ||
                                       (abs(pi_gd_gd_gd_mother_id) == KS_ID)
                                 : false;

        if (from_ks) continue;
        count_tm++;

        // Conditional cuts
        if (!probe_hasmuon) continue;
        count_cond++;

        // Calib sample cuts
        // See tables 36 and 48 in LHCb-PUB-2016-005
        // See also
        // https://lhcb-pid-wgp-plots.web.cern.ch/lhcb-pid-wgp-plots/Run2/ and
        // https://gitlab.cern.ch/lhcb-datapkg/WG/SemilepConfig/-/blob/master/options/Filter_Dstar2D02KpiNoPID_2016MC.py?ref_type=heads
        if ((tag_p < 2000.) || (tag_pt < 250.)) continue;
        if ((probe_p < 2000.) || (probe_pt < 250.)) continue;
        if ((spi_p < 1000.) || (spi_pt < 100.)) continue;
        if (k_track_chi2ndof > 3. || pi_track_chi2ndof > 3. ||
            spi_track_chi2ndof > 3.)
          continue;

        double k_pz, pi_pz;
        if (probe == "k") {
          k_pz  = probe_pz;
          pi_pz = tag_pz;
        } else {
          k_pz  = tag_pz;
          pi_pz = probe_pz;
        }

        ROOT::Math::PxPyPzMVector p_k(k_px * 0.001, k_py * 0.001, k_pz * 0.001,
                                      K_M);
        ROOT::Math::PxPyPzMVector p_pi(pi_px * 0.001, pi_py * 0.001,
                                       pi_pz * 0.001, PI_M);

        if ((std::max(probe_pt, tag_pt) < 1000.) || ((p_k + p_pi).Pt() < 1.5))
          continue;

        // Veto wrong-mass hypotheses
        // pi pi
        p_k.SetM(PI_M);
        p_pi.SetM(PI_M);
        if (std::abs((p_k + p_pi).M() - D0_M) < 0.025) continue;

        // pi K (wrong sign)
        p_k.SetM(PI_M);
        p_pi.SetM(K_M);
        if (std::abs((p_k + p_pi).M() - D0_M) < 0.025) continue;

        // K K
        p_k.SetM(K_M);
        p_pi.SetM(K_M);
        if (std::abs((p_k + p_pi).M() - D0_M) < 0.025) continue;

        count_calib_sel++;

        // MuonUnbiased equivalent to L0+HLT1 TIS. See MuonUnbiased defition
        // at
        // https://gitlab.cern.ch/lhcb-datapkg/WG/PIDCalib/-/blob/master/scriptsR2/makeTuples_pp_2016_reprocessing.py#L71
        // and
        // https://mattermost.web.cern.ch/lhcb/pl/893yre484jggigooti5u3gqb8w

        // Not applying MuonUnbiased and GHOSTPROB conditional cuts in MC
        // since some kinematical bins have very low statistics

        if (probe_p < 3000. || probe_p >= 100000. || ntracks >= 600) continue;

        const double probe_eta =
            0.5 * log((probe_p + probe_pz) / (probe_p - probe_pz));
        if (probe_eta < 1.7 || probe_eta >= 5.0) continue;
        count_range++;

        // Determine kinematical bin
        const int kin_bin = histo_binning.FindBin(probe_p, probe_eta, ntracks);
        histo_binning.GetBinXYZ(kin_bin, p_bin, eta_bin, ntrks_bin);

        // Fill histograms
        const double dm = (dst_m - d0_m) * 0.001;
        d0_m            = d0_m * 0.001;

        const bool in_extended_window =
            in_var_range(d0_m_var, d0_m) && in_var_range(dm_var, dm);
        if (!in_extended_window) continue;

        const bool reduced_fit_range = (probe == "pi") && (p_bin <= 2);

        const bool in_fit_window = reduced_fit_range
                                       ? in_range(D0_M_min, d0_m, 1.900) &&
                                             in_range(DM_min, dm, DM_max)
                                       : in_range(D0_M_min, d0_m, D0_M_max) &&
                                             in_range(DM_min, dm, DM_max);
        if (in_fit_window) count_mw++;

        // PID
        bool pid_ok;
        if (fake_mu) {
          // For FAKE_MU, we calculate the complementary efficiency so that
          // the "passed" sample always corresponds to the K/pi misid case
          pid_ok = probe_ismuon;
        } else if (vmu) {
          pid_ok = probe_ismuon && probe_dllmu > 2.0 && probe_dlle < 1.0 &&
                   probe_mu_ubdt < 0.25;
        } else {
          pid_ok = probe_ismuon && probe_dllmu > 2.0 && probe_dlle < 1.0 &&
                   probe_mu_ubdt > 0.25;
        }

        d0_m_var.setVal(d0_m);
        dm_var.setVal(dm);

        if (pid_ok) {
          // Fill "passed" samples
          datasets_d0_bkg_mc_passed.at(d0_decay)[(p_bin - 1)].addFast(fit_vars);
          if (in_fit_window) {
            d0_bkg_counts_passed[probe][p_bin - 1][d0_decay]++;
          }
        } else {
          // Fill "failed" samples
          datasets_d0_bkg_mc_failed.at(d0_decay)[(p_bin - 1)].addFast(fit_vars);
          if (in_fit_window) {
            d0_bkg_counts_failed[probe][p_bin - 1][d0_decay]++;
          }
        }
      }

      auto stop_mc_loop = high_resolution_clock::now();

      const int duration_mc_loop =
          duration_cast<seconds>(stop_mc_loop - start_mc_loop).count();
      cout << "INFO MC loop took " << format_time(duration_mc_loop) << endl;

      cout << "INFO MC cutflow:\n";
      cout << " - Total MC events:              " << entries_mc << "\n";
      cout << " - After offline truth-matching: " << count_tm << "("
           << count_tm * 100. / entries_mc << "%)\n";
      cout << " - After conditional cuts:       " << count_cond << "("
           << count_cond * 100. / count_tm << "%)\n";
      cout << " - After PIDCalib cuts:          " << count_calib_sel << "("
           << count_calib_sel * 100. / count_cond << "%)\n";
      cout << " - After range cuts:             " << count_range << "("
           << count_range * 100. / count_calib_sel << "%)\n";
      cout << " - After mass window cuts:       " << count_mw << "("
           << count_mw * 100. / count_range << "%)\n";
      cout << " Overall: " << "(" << count_mw * 100. / entries_mc << "%)\n"
           << endl;

      cout << "INFO Normalization scale factor: " << w_d0_decays.at(d0_decay)
           << "\n"
           << endl;
    }

    // Now loop over calib samples and make fits

    // Define calib datasets
    RooDataSet *datasets_calib_passed[N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] = {
        {{nullptr}}};
    RooDataSet *datasets_calib_failed[N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] = {
        {{nullptr}}};

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

    for (unsigned p = 0; p < calib_paths.size(); p++) {
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

        // Define variable to access input ntuples
        double dst_m, d0_m, probe_p, probe_eta, probe_dllmu, probe_dlle,
            probe_ghostprob, probe_mu_unbiased, ntracks_calib;

        float probe_mu_ubdt;

        bool probe_ismuon, probe_hasmuon;

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

        int count_tis = 0, count_cond = 0, count_range = 0;

        int eta_bin = 0, p_bin = 0, ntrks_bin = 0;

        const int entries_calib = ch_calib.GetEntries();
        cout << "INFO Starting PIDCalib event loop over " << entries_calib
             << " entries in tree " << tree_name << endl;
        auto start_calib_loop = high_resolution_clock::now();
        for (int evt = 0; evt < entries_calib; evt++) {
          ch_calib.GetEntry(evt);

          // Trigger decorrelation
          if (!probe_mu_unbiased) continue;
          count_tis++;

          // Conditional cuts
          if (!probe_hasmuon || probe_ghostprob > 0.5) continue;
          count_cond++;

          // Kinematical range cuts
          if (probe_p < 3000. || probe_p >= 100000. || ntracks_calib >= 600)
            continue;

          if (probe_eta < 1.7 || probe_eta >= 5.0) continue;
          count_range++;

          // Fill histograms

          d0_m_var.setVal(d0_m * 0.001);
          dm_var.setVal((dst_m - d0_m) * 0.001);

          // Determine kinematical bin
          const int kin_bin =
              histo_binning.FindBin(probe_p, probe_eta, ntracks_calib);

          histo_binning.GetBinXYZ(kin_bin, p_bin, eta_bin, ntrks_bin);

          bool pid_ok;
          if (fake_mu) {
            // For FAKE_MU, we calculate the complementary efficiency so
            // that the "passed" sample always corresponds to the K/pi misid
            // case
            pid_ok = probe_ismuon;
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
            duration_cast<seconds>(stop_calib_loop - start_calib_loop).count();
        cout << "INFO Calib loop took " << format_time(duration_calib_loop)
             << endl;

        cout << "INFO Calib cutflow:\n";
        cout << " - Total MC events:              " << entries_calib << "\n";
        cout << " - After probe unbiased  cut:    " << count_tis << "("
             << count_tis * 100. / entries_calib << "%)\n";
        cout << " - After conditional cuts:       " << count_cond << "("
             << count_cond * 100. / count_tis << "%)\n";
        cout << " - After range cuts:             " << count_range << "("
             << count_range * 100. / count_cond << "%)\n";
        cout << " Overall: " << "(" << count_range * 100. / entries_calib
             << "%)\n"
             << endl;
      }
    }

    // Create histograms to store efficiency of PIDCalib mass windows
    cout << "INFO Initializing efficiency histograms " << endl;
    TString suffix = TString::Format("%s_%s", year.c_str(), probe.c_str());

    TH3D effs_mw_passed(histo_binning);
    effs_mw_passed.SetName("effs_mw_passed_" + suffix);

    TH3D effs_mw_failed(histo_binning);
    effs_mw_failed.SetName("effs_mw_failed_" + suffix);

    // Check distribution of fitted parameters
    RooDataSet ds_params_calib("ds_params_calib", "ds_params_calib",
                               params_calib);

    TH3D fit_status_calib(histo_binning);
    fit_status_calib.SetName("fit_status_calib_" + suffix);
    fit_status_calib.SetMinimum(-0.1);
    fit_status_calib.SetMaximum(12.);

    TH3D fit_cov_qual_calib(histo_binning);
    fit_cov_qual_calib.SetName("fit_cov_qual_calib_" + suffix);
    fit_cov_qual_calib.SetMinimum(-0.1);
    fit_cov_qual_calib.SetMaximum(6.);

    TH1D calib_retries("calib_retries", "calib fit #retries;#retries;;",
                       max_fix_reattempts + 1, 0., max_fix_reattempts + 1.);

    ofile.cd();

    cout << "INFO Calculating efficiencies " << endl;

    // Create histogram to hold PID efficiency
    // Use "eff" as name to match PIDCalib output
    TH3D histo_pid(histo_binning);
    histo_pid.SetName("eff");

    TH3D histo_pid_raw(histo_binning);
    histo_pid_raw.SetName("eff_raw");

    // Create histogram to store fitted DiF fractions
    TH3D histo_f_dif(histo_binning);
    histo_f_dif.SetName("f_dif_" + TString(probe));

    for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
      for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
        for (int ntrks_idx = 0; ntrks_idx < N_BINS_NTRACKS; ntrks_idx++) {
          // Use ntrks_idx as inner loop because we combine ntrack bins in the
          // DiF samples

          suffix.Form("%s_%s_%d_%d_%d", year.c_str(), probe.c_str(), ntrks_idx,
                      eta_idx, p_idx);
          const TString tag = TString::Format(
              "(%.0f < nTracks < %.0f, %.1f < #eta < %.1f, %.0f < p < %.0f)",
              BINS_NTRACKS[ntrks_idx], BINS_NTRACKS[ntrks_idx + 1],
              BINS_ETA[eta_idx], BINS_ETA[eta_idx + 1], BINS_P[p_idx],
              BINS_P[p_idx + 1]);

          const TString fit_dir_path = opath + "/figs/fits/" + suffix;
          gSystem->mkdir(fit_dir_path);

          // Get kinematical bin
          const int kin_bin =
              histo_binning.GetBin(p_idx + 1, eta_idx + 1, ntrks_idx + 1);

          // Get datasets
          const auto &dataset_calib_passed =
              datasets_calib_passed[ntrks_idx][eta_idx][p_idx];
          const auto &dataset_calib_failed =
              datasets_calib_failed[ntrks_idx][eta_idx][p_idx];

          if ((probe == "pi") && (p_idx <= 1)) {
            // Use reduced fit range excluding higher D0 mass region to avoid
            // upper mass threshold in some kinematic bins
            d0_m_var.setRange("fitRange", 1.825, 1.900);
            dm_var.setRange("fitRange", 0.141, 0.153);
          } else {
            // Use full fit range
            d0_m_var.setRange("fitRange", 1.825, 1.910);
            dm_var.setRange("fitRange", 0.141, 0.153);
          }

          const int n_calib_passed =
              dataset_calib_passed->sumEntries("", "fitRange");
          const int n_calib_failed =
              dataset_calib_failed->sumEntries("", "fitRange");

          // Recover fit limits
          const double d0m_range_max = d0_m_var.getMax("fitRange");
          const double d0m_range_min = d0_m_var.getMin("fitRange");
          const double dm_range_max  = dm_var.getMax("fitRange");
          const double dm_range_min  = dm_var.getMin("fitRange");

          // Plot 2D distributions in calib sample
          d0_m_var.setRange(D0_M_min, D0_M_max);
          dm_var.setRange(DM_min, DM_max);

          unique_ptr<TH2F> th2_calib_passed(
              dataset_calib_passed->createHistogram(
                  d0_m_var, dm_var, 30, 30, "", "th2_calib_passed_" + suffix));
          unique_ptr<TH2F> th2_calib_failed(
              dataset_calib_failed->createHistogram(
                  d0_m_var, dm_var, 60, 60, "", "th2_calib_failed_" + suffix));

          c_single.cd();
          th2_calib_passed->Draw("LEGO2Z 0");
          c_single.SaveAs(fit_dir_path + "/th2_" + suffix + "_passed.pdf");
          th2_calib_failed->Draw("LEGO2Z 0");
          c_single.SaveAs(fit_dir_path + "/th2_" + suffix + "_failed.pdf");

          d0_m_var.setRange(extended_d0_m_min, extended_d0_m_max);
          dm_var.setRange(extended_dm_min, extended_dm_max);

          d0_m_shift.setVal(0.);
          dm_shift.setVal(0.);
          d0_m_scale.setVal(1.);
          dm_scale.setVal(1.);

          ////////////////////////////
          // Non-DiF d0_m template  //
          ////////////////////////////

          cout << "\nINFO Building D0 MC template without dif " << suffix
               << endl;

          const auto &ds_mc_passed_nondif =
              datasets_mc_passed_nondif[ntrks_idx][eta_idx][p_idx];
          const auto &ds_mc_failed_nondif =
              datasets_mc_failed_nondif[ntrks_idx][eta_idx][p_idx];

          // The mass distributions for non-dif events are not correlated with
          // muon PID, so let's merge the two datasets
          RooDataSet ds_mc_nondif_merged("ds_mc_nondif_merged_" + suffix,
                                         "ds_mc_nondif_merged_" + suffix,
                                         fit_vars);
          ds_mc_nondif_merged.append(*ds_mc_passed_nondif);
          ds_mc_nondif_merged.append(*ds_mc_failed_nondif);

          cout << "INFO Building RooKeysPdf d0_model_passed_nondif with "
               << ds_mc_nondif_merged.numEntries() << " entries" << endl;
          RooKeysPdf d0_model_passed_nondif(
              "d0_model_passed_nondif_" + suffix,
              "d0_model_passed_nondif_" + suffix, d0_m_var, d0_m_var, d0_m_pdg,
              d0_m_scale, d0_m_shift, ds_mc_nondif_merged, RooKeysPdf::NoMirror,
              1.6, true, n_bins_keys);

          cout << "INFO Building RooKeysPdf d0_model_failed_nondif with "
               << ds_mc_nondif_merged.numEntries() << " entries" << endl;
          RooKeysPdf d0_model_failed_nondif(
              "d0_model_failed_nondif_" + suffix,
              "d0_model_failed_nondif_" + suffix, d0_m_var, d0_m_var, d0_m_pdg,
              d0_m_scale, d0_m_shift, ds_mc_nondif_merged, RooKeysPdf::NoMirror,
              1.6, true, n_bins_keys);

          // Plot fit results
          // PDF from merged dataset is compared to separate distribution
          // in individual samples
          unique_ptr<RooPlot> frame_d0_passed_nondif(
              d0_m_var.frame(Title("D0 M Passed " + tag)));
          unique_ptr<RooPlot> frame_d0_failed_nondif(
              d0_m_var.frame(Title("D0 M Failed " + tag)));
          ds_mc_passed_nondif->plotOn(frame_d0_passed_nondif.get(),
                                      Binning("bins_histos_d0_m_passed"));
          ds_mc_failed_nondif->plotOn(frame_d0_failed_nondif.get(),
                                      Binning("bins_histos_d0_m_failed"));
          d0_model_passed_nondif.plotOn(frame_d0_passed_nondif.get(),
                                        LineWidth(1), Range("data"),
                                        LineColor(kRed), VLines());
          d0_model_passed_nondif.plotOn(frame_d0_passed_nondif.get(),
                                        LineWidth(2), NormRange("data"));
          d0_model_failed_nondif.plotOn(frame_d0_failed_nondif.get(),
                                        LineWidth(1), Range("data"),
                                        LineColor(kRed), VLines());
          d0_model_failed_nondif.plotOn(frame_d0_failed_nondif.get(),
                                        LineWidth(2), NormRange("data"));

          // Note: All MC PDFs built with RooKeysPdf are plotted with
          // NormRange("data") i.e. normalization is calculated only in fit
          // range. This is because the RooKeysPdf estimate gets weird near
          // the boundaries (and we don't use that region in the fit anyways).

          c_double.cd(1);
          frame_d0_passed_nondif->Draw();
          c_double.cd(2);
          frame_d0_failed_nondif->Draw();
          c_double.SaveAs(fit_dir_path + "/d0_m_" + suffix + "_nondif.pdf");

          ///////////////////////
          // DiF d0_m template //
          ///////////////////////

          cout << "\nINFO Building D0 MC template with dif " << suffix << endl;

          const auto &ds_mc_passed_dif = datasets_mc_passed_dif[eta_idx][p_idx];
          const auto &ds_mc_failed_dif = datasets_mc_failed_dif[eta_idx][p_idx];

          cout << "INFO Building RooKeysPdf d0_model_passed_dif with "
               << ds_mc_passed_dif->numEntries() << " entries" << endl;
          RooKeysPdf d0_model_passed_dif(
              "d0_model_passed_dif_" + suffix, "d0_model_passed_dif_" + suffix,
              d0_m_var, d0_m_var, d0_m_pdg, d0_m_scale, d0_m_shift,
              *ds_mc_passed_dif, RooKeysPdf::NoMirror, 1.8, true, n_bins_keys);

          cout << "INFO Building RooKeysPdf d0_model_failed_dif with "
               << ds_mc_failed_dif->numEntries() << " entries" << endl;
          RooKeysPdf d0_model_failed_dif(
              "d0_model_failed_dif_" + suffix, "d0_model_failed_dif_" + suffix,
              d0_m_var, d0_m_var, d0_m_pdg, d0_m_scale, d0_m_shift,
              *ds_mc_failed_dif, RooKeysPdf::NoMirror, 1.8, true, n_bins_keys);

          // Plot fit results
          unique_ptr<RooPlot> frame_d0_passed_dif(
              d0_m_var.frame(Title("D0 M Passed " + tag)));
          unique_ptr<RooPlot> frame_d0_failed_dif(
              d0_m_var.frame(Title("D0 M Failed " + tag)));
          ds_mc_passed_dif->plotOn(frame_d0_passed_dif.get(),
                                   Binning("bins_histos_d0_m_passed"));
          ds_mc_failed_dif->plotOn(frame_d0_failed_dif.get(),
                                   Binning("bins_histos_d0_m_failed"));
          d0_model_passed_dif.plotOn(frame_d0_passed_dif.get(), LineWidth(1),
                                     Range("data"), LineColor(kRed), VLines());
          d0_model_passed_dif.plotOn(frame_d0_passed_dif.get(), LineWidth(2),
                                     NormRange("data"));
          d0_model_failed_dif.plotOn(frame_d0_failed_dif.get(), LineWidth(1),
                                     Range("data"), LineColor(kRed), VLines());
          d0_model_failed_dif.plotOn(frame_d0_failed_dif.get(), LineWidth(2),
                                     NormRange("data"));

          c_double.cd(1);
          frame_d0_passed_dif->Draw();
          c_double.cd(2);
          frame_d0_failed_dif->Draw();
          c_double.SaveAs(fit_dir_path + "/d0_m_" + suffix + "_dif.pdf");

          /////////////////////////
          // Non-DiF dm template //
          /////////////////////////

          cout << "\nINFO Building dm MC template without dif " << suffix
               << endl;

          cout << "INFO Building RooKeysPdf dm_model_passed_nondif with "
               << ds_mc_nondif_merged.numEntries() << " entries" << endl;
          RooKeysPdf dm_model_passed_nondif(
              "dm_model_passed_nondif_" + suffix,
              "dm_model_passed_nondif_" + suffix, dm_var, dm_var, dm_pdg,
              dm_scale, dm_shift, ds_mc_nondif_merged, RooKeysPdf::NoMirror,
              1.6, true, n_bins_keys);

          cout << "INFO Building RooKeysPdf dm_model_failed_nondif with "
               << ds_mc_nondif_merged.numEntries() << " entries" << endl;
          RooKeysPdf dm_model_failed_nondif(
              "dm_model_failed_nondif_" + suffix,
              "dm_model_failed_nondif_" + suffix, dm_var, dm_var, dm_pdg,
              dm_scale, dm_shift, ds_mc_nondif_merged, RooKeysPdf::NoMirror,
              1.6, true, n_bins_keys);

          // Plot fit results
          // PDF from merged dataset is compared to separate distribution
          // in individual samples
          unique_ptr<RooPlot> frame_dm_passed_nondif(
              dm_var.frame(Title("dm Passed " + tag)));
          unique_ptr<RooPlot> frame_dm_failed_nondif(
              dm_var.frame(Title("dm Failed " + tag)));
          ds_mc_passed_nondif->plotOn(frame_dm_passed_nondif.get(),
                                      Binning("bins_histos_dm_passed"));
          ds_mc_failed_nondif->plotOn(frame_dm_failed_nondif.get(),
                                      Binning("bins_histos_dm_failed"));
          dm_model_passed_nondif.plotOn(frame_dm_passed_nondif.get(),
                                        LineWidth(1), Range("data"),
                                        LineColor(kRed), VLines());
          dm_model_passed_nondif.plotOn(frame_dm_passed_nondif.get(),
                                        LineWidth(2), NormRange("data"));
          dm_model_failed_nondif.plotOn(frame_dm_failed_nondif.get(),
                                        LineWidth(1), Range("data"),
                                        LineColor(kRed), VLines());
          dm_model_failed_nondif.plotOn(frame_dm_failed_nondif.get(),
                                        LineWidth(2), NormRange("data"));

          c_double.cd(1);
          frame_dm_passed_nondif->Draw();
          c_double.cd(2);
          frame_dm_failed_nondif->Draw();
          c_double.SaveAs(fit_dir_path + "/dm_" + suffix + "_nondif.pdf");

          /////////////////////
          // DiF dm template //
          /////////////////////

          cout << "\nINFO Building dm MC template with dif " << suffix << endl;

          cout << "INFO Building RooKeysPdf dm_model_passed_dif with "
               << ds_mc_passed_dif->numEntries() << " entries" << endl;
          RooKeysPdf dm_model_passed_dif(
              "dm_model_passed_dif_" + suffix, "dm_model_passed_dif_" + suffix,
              dm_var, dm_var, dm_pdg, dm_scale, dm_shift, *ds_mc_passed_dif,
              RooKeysPdf::NoMirror, 1.6, true, n_bins_keys);

          cout << "INFO Building RooKeysPdf dm_model_failed_dif with "
               << ds_mc_failed_dif->numEntries() << " entries" << endl;
          RooKeysPdf dm_model_failed_dif(
              "dm_model_failed_dif_" + suffix, "dm_model_failed_dif_" + suffix,
              dm_var, dm_var, dm_pdg, dm_scale, dm_shift, *ds_mc_failed_dif,
              RooKeysPdf::NoMirror, 1.6, true, n_bins_keys);

          // https://root-forum.cern.ch/t/fit-an-offset-for-roohistpdf/9641

          // Plot fit results
          unique_ptr<RooPlot> frame_dm_passed_dif(
              dm_var.frame(Title("dm Passed " + tag)));
          unique_ptr<RooPlot> frame_dm_failed_dif(
              dm_var.frame(Title("dm Failed " + tag)));
          ds_mc_passed_dif->plotOn(frame_dm_passed_dif.get(),
                                   Binning("bins_histos_dm_passed"));
          ds_mc_failed_dif->plotOn(frame_dm_failed_dif.get(),
                                   Binning("bins_histos_dm_failed"));
          dm_model_passed_dif.plotOn(frame_dm_passed_dif.get(), LineWidth(1),
                                     Range("data"), LineColor(kRed), VLines());
          dm_model_passed_dif.plotOn(frame_dm_passed_dif.get(), LineWidth(2),
                                     NormRange("data"));
          dm_model_failed_dif.plotOn(frame_dm_failed_dif.get(), LineWidth(1),
                                     Range("data"), LineColor(kRed), VLines());
          dm_model_failed_dif.plotOn(frame_dm_failed_dif.get(), LineWidth(2),
                                     NormRange("data"));

          c_double.cd(1);
          frame_dm_passed_dif->Draw();
          c_double.cd(2);
          frame_dm_failed_dif->Draw();
          c_double.SaveAs(fit_dir_path + "/dm_" + suffix + "_dif.pdf");

          ////////////////////////////////////////////////////
          // d0_m templates for mis-reconstructed D0 decays //
          ////////////////////////////////////////////////////

          cout << "\nINFO Building d0_m MC template for D0 bkg " << suffix
               << endl;

          unordered_map<string, RooKeysPdf> d0m_pdfs_d0_bkg_passed;
          unordered_map<string, RooKeysPdf> d0m_pdfs_d0_bkg_failed;
          unordered_map<string, RooRealVar> d0m_fs_d0_bkg_passed;
          unordered_map<string, RooRealVar> d0m_fs_d0_bkg_failed;
          unordered_map<string, RooRealVar> d0m_fs_extended_d0_bkg_passed;
          unordered_map<string, RooRealVar> d0m_fs_extended_d0_bkg_failed;
          RooArgList d0m_pdf_list_d0_bkg_passed, d0m_pdf_list_d0_bkg_failed,
              d0m_f_list_d0_bkg_passed, d0m_f_list_d0_bkg_failed,
              d0m_f_extended_list_d0_bkg_passed,
              d0m_f_extended_list_d0_bkg_failed;

          RooDataSet ds_d0m_d0_bkg_passed_sum(
              "ds_d0m_d0_bkg_passed_sum_" + suffix,
              "ds_d0m_d0_bkg_passed_sum_" + suffix, fit_vars_w,
              WeightVar(w_d0_bkg));
          RooDataSet ds_d0m_d0_bkg_failed_sum(
              "ds_d0m_d0_bkg_failed_sum_" + suffix,
              "ds_d0m_d0_bkg_failed_sum_" + suffix, fit_vars_w,
              WeightVar(w_d0_bkg));

          unordered_map<string, TH1D> d0m_histos_d0_bkg_passed;
          unordered_map<string, TH1D> d0m_histos_d0_bkg_failed;
          TH1D    h_d0m_d0_bkg_passed("h_d0m_d0_bkg_passed_" + suffix,
                                      "h_d0m_d0_bkg_passed_" + suffix,
                                      bins_histos_d0_m_d0_bkg.numBins(),
                                      bins_histos_d0_m_d0_bkg.lowBound(),
                                      bins_histos_d0_m_d0_bkg.highBound());
          TH1D    h_d0m_d0_bkg_failed("h_d0m_d0_bkg_failed_" + suffix,
                                      "h_d0m_d0_bkg_failed_" + suffix,
                                      bins_histos_d0_m_d0_bkg.numBins(),
                                      bins_histos_d0_m_d0_bkg.lowBound(),
                                      bins_histos_d0_m_d0_bkg.highBound());
          THStack ths_d0m_d0_bkg_passed("ths_d0m_d0_bkg_passed_" + suffix,
                                        "ths_d0m_d0_bkg_passed_" + suffix);
          THStack ths_d0m_d0_bkg_failed("ths_d0m_d0_bkg_failed_" + suffix,
                                        "ths_d0m_d0_bkg_failed_" + suffix);

          double sum_f_passed_d0m = 0., sum_f_failed_d0m = 0.,
                 sum_f_ext_passed_d0m = 0., sum_f_ext_failed_d0m = 0.;

          double d0_bkg_sum_passed = 0., d0_bkg_sum_failed = 0.,
                 d0_bkg_sum_extended_passed = 0.,
                 d0_bkg_sum_extended_failed = 0.;

          for (const auto &d0_decay : d0_bkg_decays) {
            // Total counts over extended range
            d0_bkg_sum_extended_passed +=
                w_d0_decays.at(d0_decay) *
                datasets_d0_bkg_mc_passed[d0_decay][p_idx].numEntries();
            d0_bkg_sum_extended_failed +=
                w_d0_decays.at(d0_decay) *
                datasets_d0_bkg_mc_failed[d0_decay][p_idx].numEntries();
            // Total counts in pidcalib mass window
            d0_bkg_sum_passed += w_d0_decays.at(d0_decay) *
                                 d0_bkg_counts_passed[probe][p_idx][d0_decay];
            d0_bkg_sum_failed += w_d0_decays.at(d0_decay) *
                                 d0_bkg_counts_failed[probe][p_idx][d0_decay];
          }

          for (const auto &d0_decay : d0_bkg_decays) {
            cout << "\nINFO Producing d0_M PDFs for " << d0_decay << endl;
            const TString suffix_bkg = d0_decay + "_" + suffix;

            RooDataSet &ds_passed = datasets_d0_bkg_mc_passed[d0_decay][p_idx];
            RooDataSet &ds_failed = datasets_d0_bkg_mc_failed[d0_decay][p_idx];

            const double n_passed_d0_decay =
                ds_passed.numEntries() * w_d0_decays.at(d0_decay);
            const double n_failed_d0_decay =
                ds_failed.numEntries() * w_d0_decays.at(d0_decay);

            w_d0_bkg.setVal(w_d0_decays.at(d0_decay) * n_passed_d0_decay /
                            (n_passed_d0_decay + n_failed_d0_decay));
            append_to_dataset(ds_passed, ds_d0m_d0_bkg_passed_sum, w_d0_bkg);
            append_to_dataset(ds_failed, ds_d0m_d0_bkg_passed_sum, w_d0_bkg);

            w_d0_bkg.setVal(w_d0_decays.at(d0_decay) * n_failed_d0_decay /
                            (n_passed_d0_decay + n_failed_d0_decay));
            append_to_dataset(ds_passed, ds_d0m_d0_bkg_failed_sum, w_d0_bkg);
            append_to_dataset(ds_failed, ds_d0m_d0_bkg_failed_sum, w_d0_bkg);

            RooDataSet ds_d0m_d0_bkg("ds_d0m_d0_bkg_" + suffix_bkg,
                                     "ds_d0m_d0_bkg_" + suffix_bkg, fit_vars);

            ds_d0m_d0_bkg.append(ds_passed);
            ds_d0m_d0_bkg.append(ds_failed);

            unique_ptr<TH1D> h_passed(static_cast<TH1D *>(
                ds_passed.createHistogram("h_passed_" + suffix_bkg, d0_m_var,
                                          Binning(bins_histos_d0_m_d0_bkg))));
            unique_ptr<TH1D> h_failed(static_cast<TH1D *>(
                ds_failed.createHistogram("h_failed_" + suffix_bkg, d0_m_var,
                                          Binning(bins_histos_d0_m_d0_bkg))));

            h_passed->SetFillColor(color_ids_d0_decays.at(d0_decay));
            h_failed->SetFillColor(color_ids_d0_decays.at(d0_decay));

            d0m_histos_d0_bkg_passed.emplace(
                d0_decay, TH1D("h_d0m_d0_bkg_passed_merged_" + suffix_bkg,
                               "h_d0m_d0_bkg_passed_merged_" + suffix_bkg,
                               bins_histos_d0_m_d0_bkg.numBins(),
                               bins_histos_d0_m_d0_bkg.lowBound(),
                               bins_histos_d0_m_d0_bkg.highBound()));
            d0m_histos_d0_bkg_failed.emplace(
                d0_decay, TH1D("h_d0m_d0_bkg_failed_merged_" + suffix_bkg,
                               "h_d0m_d0_bkg_failed_merged_" + suffix_bkg,
                               bins_histos_d0_m_d0_bkg.numBins(),
                               bins_histos_d0_m_d0_bkg.lowBound(),
                               bins_histos_d0_m_d0_bkg.highBound()));

            TH1D &h_passed_merged = d0m_histos_d0_bkg_passed[d0_decay];
            TH1D &h_failed_merged = d0m_histos_d0_bkg_failed[d0_decay];

            h_passed_merged.Add(h_passed.get(), h_failed.get());
            h_failed_merged.Add(h_passed.get(), h_failed.get());

            h_passed_merged.SetFillColor(color_ids_d0_decays.at(d0_decay));
            h_failed_merged.SetFillColor(color_ids_d0_decays.at(d0_decay));

            h_passed_merged.Scale(n_passed_d0_decay /
                                  h_passed_merged.GetEntries());
            h_failed_merged.Scale(n_failed_d0_decay /
                                  h_failed_merged.GetEntries());
            ths_d0m_d0_bkg_passed.Add(&h_passed_merged);
            ths_d0m_d0_bkg_failed.Add(&h_failed_merged);
            h_d0m_d0_bkg_passed.Add(&h_passed_merged);
            h_d0m_d0_bkg_failed.Add(&h_failed_merged);

            cout << "INFO Building RooKeysPdf d0_model_d0_bkg_passed with "
                 << ds_d0m_d0_bkg.numEntries() << " entries" << endl;
            d0m_pdfs_d0_bkg_passed.emplace(
                d0_decay,
                RooKeysPdf("d0_model_d0_bkg_passed_" + suffix_bkg,
                           "d0_model_d0_bkg_passed_" + suffix_bkg, d0_m_var,
                           d0_m_var, d0_m_pdg, d0_m_scale, d0_m_shift,
                           ds_d0m_d0_bkg, RooKeysPdf::NoMirror, 1.5, true,
                           n_bins_keys));

            cout << "INFO Building RooKeysPdf d0_model_d0_bkg_failed with "
                 << ds_d0m_d0_bkg.numEntries() << " entries" << endl;
            d0m_pdfs_d0_bkg_failed.emplace(
                d0_decay,
                RooKeysPdf("d0_model_d0_bkg_failed_" + suffix_bkg,
                           "d0_model_d0_bkg_failed_" + suffix_bkg, d0_m_var,
                           d0_m_var, d0_m_pdg, d0_m_scale, d0_m_shift,
                           ds_d0m_d0_bkg, RooKeysPdf::NoMirror, 1.5, true,
                           n_bins_keys));

            unique_ptr<RooPlot> frame_d0_mc_d0_bkg_passed(
                d0_m_var.frame(Title("d0_m Passed " + tag)));
            unique_ptr<RooPlot> frame_d0_mc_d0_bkg_failed(
                d0_m_var.frame(Title("d0_m Failed " + tag)));
            unique_ptr<RooPlot> frame_d0_mc_d0_bkg_total(
                d0_m_var.frame(Title("d0_m Total " + tag)));

            ds_passed.plotOn(frame_d0_mc_d0_bkg_passed.get(),
                             Binning(bins_histos_d0_m_d0_bkg));
            ds_failed.plotOn(frame_d0_mc_d0_bkg_failed.get(),
                             Binning(bins_histos_d0_m_d0_bkg));
            ds_d0m_d0_bkg.plotOn(frame_d0_mc_d0_bkg_total.get(),
                                 Binning(bins_histos_d0_m_d0_bkg));
            d0m_pdfs_d0_bkg_passed.at(d0_decay).plotOn(
                frame_d0_mc_d0_bkg_passed.get(), LineWidth(1), LineColor(kRed),
                Range("fitRange"), VLines());
            d0m_pdfs_d0_bkg_passed.at(d0_decay).plotOn(
                frame_d0_mc_d0_bkg_passed.get(), LineWidth(2),
                LineColor(color_ids_d0_decays.at(d0_decay)),
                NormRange("fitRange"));
            d0m_pdfs_d0_bkg_failed.at(d0_decay).plotOn(
                frame_d0_mc_d0_bkg_failed.get(), LineWidth(1), LineColor(kRed),
                Range("fitRange"), VLines());
            d0m_pdfs_d0_bkg_failed.at(d0_decay).plotOn(
                frame_d0_mc_d0_bkg_failed.get(), LineWidth(2),
                LineColor(color_ids_d0_decays.at(d0_decay)),
                NormRange("fitRange"));
            d0m_pdfs_d0_bkg_failed.at(d0_decay).plotOn(
                frame_d0_mc_d0_bkg_total.get(), LineWidth(1), LineColor(kRed),
                Range("fitRange"), VLines());
            d0m_pdfs_d0_bkg_failed.at(d0_decay).plotOn(
                frame_d0_mc_d0_bkg_total.get(), LineWidth(2),
                LineColor(color_ids_d0_decays.at(d0_decay)),
                NormRange("fitRange"));

            c3.cd(1);
            frame_d0_mc_d0_bkg_passed->Draw();
            c3.cd(2);
            frame_d0_mc_d0_bkg_failed->Draw();
            c3.cd(3);
            frame_d0_mc_d0_bkg_total->Draw();
            c3.SaveAs(fit_dir_path + "/d0_m_" + suffix + "_d0_bkg_keys_" +
                      d0_decay + ".pdf");

            d0m_pdf_list_d0_bkg_passed.add(d0m_pdfs_d0_bkg_passed.at(d0_decay));
            d0m_pdf_list_d0_bkg_failed.add(d0m_pdfs_d0_bkg_failed.at(d0_decay));

            if (d0_decay != d0_bkg_decays.back()) {
              // Here we calculate fit fractions, so the this must be omitted
              // for the last component

              // Fit fractions within fit range
              d0m_fs_d0_bkg_passed.emplace(
                  d0_decay,
                  RooRealVar("d0_f_d0_bkg_passed_" + suffix_bkg,
                             "d0_f_d0_bkg_passed_" + suffix_bkg,
                             w_d0_decays.at(d0_decay) *
                                 d0_bkg_counts_passed[probe][p_idx][d0_decay] /
                                 d0_bkg_sum_passed,
                             ""));

              d0m_fs_d0_bkg_failed.emplace(
                  d0_decay,
                  RooRealVar("d0_f_d0_bkg_failed_" + suffix_bkg,
                             "d0_f_d0_bkg_failed_" + suffix_bkg,
                             w_d0_decays.at(d0_decay) *
                                 d0_bkg_counts_failed[probe][p_idx][d0_decay] /
                                 d0_bkg_sum_failed,
                             ""));

              d0m_f_list_d0_bkg_passed.add(d0m_fs_d0_bkg_passed.at(d0_decay));
              d0m_f_list_d0_bkg_failed.add(d0m_fs_d0_bkg_failed.at(d0_decay));

              // Fit fractions over extended range
              d0m_fs_extended_d0_bkg_passed.emplace(
                  d0_decay,
                  RooRealVar("d0_f_extended_d0_bkg_passed_" + suffix_bkg,
                             "d0_f_extended_d0_bkg_passed_" + suffix_bkg,
                             n_passed_d0_decay / d0_bkg_sum_extended_passed,
                             ""));

              d0m_fs_extended_d0_bkg_failed.emplace(
                  d0_decay,
                  RooRealVar("d0_f_extended_d0_bkg_failed_" + suffix_bkg,
                             "d0_f_extended_d0_bkg_failed_" + suffix_bkg,
                             n_failed_d0_decay / d0_bkg_sum_extended_failed,
                             ""));

              d0m_f_extended_list_d0_bkg_passed.add(
                  d0m_fs_extended_d0_bkg_passed.at(d0_decay));
              d0m_f_extended_list_d0_bkg_failed.add(
                  d0m_fs_extended_d0_bkg_failed.at(d0_decay));

              cout << "INFO " << d0_decay << " f passed: "
                   << d0m_fs_d0_bkg_passed.at(d0_decay).getVal() << endl;
              cout << "INFO " << d0_decay << " f failed: "
                   << d0m_fs_d0_bkg_failed.at(d0_decay).getVal() << endl;
              cout << "INFO " << d0_decay << " f ext passed: "
                   << d0m_fs_extended_d0_bkg_passed.at(d0_decay).getVal()
                   << endl;
              cout << "INFO " << d0_decay << " f ext failed: "
                   << d0m_fs_extended_d0_bkg_failed.at(d0_decay).getVal()
                   << endl;

              sum_f_passed_d0m += d0m_fs_d0_bkg_passed.at(d0_decay).getVal();
              sum_f_failed_d0m += d0m_fs_d0_bkg_failed.at(d0_decay).getVal();
              sum_f_ext_passed_d0m +=
                  d0m_fs_extended_d0_bkg_passed.at(d0_decay).getVal();
              sum_f_ext_failed_d0m +=
                  d0m_fs_extended_d0_bkg_failed.at(d0_decay).getVal();
            }
          }

          cout << "INFO " << d0_bkg_decays.back()
               << " f passed: " << 1. - sum_f_passed_d0m << endl;
          cout << "INFO " << d0_bkg_decays.back()
               << " f failed: " << 1. - sum_f_failed_d0m << endl;
          cout << "INFO " << d0_bkg_decays.back()
               << " f ext passed: " << 1. - sum_f_ext_passed_d0m << endl;
          cout << "INFO " << d0_bkg_decays.back()
               << " f ext failed: " << 1. - sum_f_ext_failed_d0m << endl;

          // Plot stacked histograms
          h_d0m_d0_bkg_passed.SetLineColor(kBlack);
          h_d0m_d0_bkg_failed.SetLineColor(kBlack);

          c_double.cd(1);
          ths_d0m_d0_bkg_passed.Draw("B HIST");
          h_d0m_d0_bkg_passed.Draw("L SAME");

          c_double.cd(2);
          ths_d0m_d0_bkg_failed.Draw("B HIST");
          h_d0m_d0_bkg_failed.Draw("L SAME");

          c_double.SaveAs(fit_dir_path + "/d0_m_" + suffix +
                          "_d0_bkg_histos_stacked.pdf");

          // Build pdfs with fit fractions calculated within the fit range
          // This is the one we use for the fits!
          RooAddPdf d0_model_d0_bkg_passed("d0_model_d0_bkg_passed_" + suffix,
                                           "d0_model_d0_bkg_passed_" + suffix,
                                           d0m_pdf_list_d0_bkg_passed,
                                           d0m_f_list_d0_bkg_passed);

          RooAddPdf d0_model_d0_bkg_failed("d0_model_d0_bkg_failed_" + suffix,
                                           "d0_model_d0_bkg_failed_" + suffix,
                                           d0m_pdf_list_d0_bkg_failed,
                                           d0m_f_list_d0_bkg_failed);

          // Build pdfs fit fractions calculated within the extended range
          RooAddPdf d0_model_d0_bkg_passed_ext(
              "d0_model_d0_bkg_passed_ext_" + suffix,
              "d0_model_d0_bkg_passed_ext_" + suffix,
              d0m_pdf_list_d0_bkg_passed, d0m_f_extended_list_d0_bkg_passed);

          RooAddPdf d0_model_d0_bkg_failed_ext(
              "d0_model_d0_bkg_failed_ext_" + suffix,
              "d0_model_d0_bkg_failed_ext_" + suffix,
              d0m_pdf_list_d0_bkg_failed, d0m_f_extended_list_d0_bkg_failed);

          // Plot PDF over extended range
          unique_ptr<RooPlot> frame_d0_mc_d0_bkg_passed_ext(
              d0_m_var.frame(Title("d0_m Passed " + tag)));
          unique_ptr<RooPlot> frame_d0_mc_d0_bkg_failed_ext(
              d0_m_var.frame(Title("d0_m Failed " + tag)));
          ds_d0m_d0_bkg_passed_sum.plotOn(frame_d0_mc_d0_bkg_passed_ext.get(),
                                          Binning(bins_histos_d0_m_d0_bkg));
          ds_d0m_d0_bkg_failed_sum.plotOn(frame_d0_mc_d0_bkg_failed_ext.get(),
                                          Binning(bins_histos_d0_m_d0_bkg));

          d0_model_d0_bkg_passed_ext.plotOn(frame_d0_mc_d0_bkg_passed_ext.get(),
                                            LineWidth(1), LineColor(kRed),
                                            Range("fitRange"), VLines());
          d0_model_d0_bkg_passed_ext.plotOn(frame_d0_mc_d0_bkg_passed_ext.get(),
                                            LineWidth(2),
                                            NormRange("fitRange"));
          d0_model_d0_bkg_failed_ext.plotOn(frame_d0_mc_d0_bkg_failed_ext.get(),
                                            LineWidth(1), LineColor(kRed),
                                            Range("fitRange"), VLines());
          d0_model_d0_bkg_failed_ext.plotOn(frame_d0_mc_d0_bkg_failed_ext.get(),
                                            LineWidth(2),
                                            NormRange("fitRange"));

          c_double.cd(1);
          frame_d0_mc_d0_bkg_passed_ext->Draw();
          c_double.cd(2);
          frame_d0_mc_d0_bkg_failed_ext->Draw();
          c_double.SaveAs(fit_dir_path + "/d0_m_" + suffix + "_d0_bkg_ext.pdf");

          // Plot PDF in fit range
          d0_m_var.setRange(d0m_range_min, d0m_range_max);
          dm_var.setRange(dm_range_min, dm_range_max);

          unique_ptr<RooPlot> frame_d0_mc_d0_bkg_passed(
              d0_m_var.frame(Title("d0_m Passed " + tag), Range("fitRange")));
          unique_ptr<RooPlot> frame_d0_mc_d0_bkg_failed(
              d0_m_var.frame(Title("d0_m Failed " + tag), Range("fitRange")));
          ds_d0m_d0_bkg_passed_sum.plotOn(frame_d0_mc_d0_bkg_passed.get(),
                                          Binning(bins_histos_d0_m_d0_bkg),
                                          CutRange("fitRange"));
          ds_d0m_d0_bkg_failed_sum.plotOn(frame_d0_mc_d0_bkg_failed.get(),
                                          Binning(bins_histos_d0_m_d0_bkg),
                                          CutRange("fitRange"));

          const double d0_bkg_norm_passed_d0m =
              ds_d0m_d0_bkg_passed_sum.sumEntries("", "fitRange");
          const double d0_bkg_norm_failed_d0m =
              ds_d0m_d0_bkg_failed_sum.sumEntries("", "fitRange");
          d0_model_d0_bkg_passed_ext.plotOn(
              frame_d0_mc_d0_bkg_passed.get(), LineWidth(2), Range("fitRange"),
              Normalization(d0_bkg_norm_passed_d0m, RooAbsReal::NumEvent),
              ProjectionRange("fitRange"));
          d0_model_d0_bkg_failed_ext.plotOn(
              frame_d0_mc_d0_bkg_failed.get(), LineWidth(2), Range("fitRange"),
              Normalization(d0_bkg_norm_failed_d0m, RooAbsReal::NumEvent),
              ProjectionRange("fitRange"));

          c_double.cd(1);
          frame_d0_mc_d0_bkg_passed->Draw();
          c_double.cd(2);
          frame_d0_mc_d0_bkg_failed->Draw();
          c_double.SaveAs(fit_dir_path + "/d0_m_" + suffix + "_d0_bkg.pdf");

          d0_m_var.setRange(extended_d0_m_min, extended_d0_m_max);
          dm_var.setRange(extended_dm_min, extended_dm_max);

          /////////////////////////////////////////////////
          // dm template for mis-reconstructed D0 decays //
          /////////////////////////////////////////////////

          cout << "\nINFO Building dm MC template for D0 bkg " << suffix
               << endl;

          unordered_map<string, RooKeysPdf> dm_pdfs_d0_bkg_passed;
          unordered_map<string, RooKeysPdf> dm_pdfs_d0_bkg_failed;
          unordered_map<string, RooRealVar> dm_fs_d0_bkg_passed;
          unordered_map<string, RooRealVar> dm_fs_d0_bkg_failed;
          unordered_map<string, RooRealVar> dm_fs_extended_d0_bkg_passed;
          unordered_map<string, RooRealVar> dm_fs_extended_d0_bkg_failed;
          RooArgList dm_pdf_list_d0_bkg_passed, dm_pdf_list_d0_bkg_failed,
              dm_f_list_d0_bkg_passed, dm_f_list_d0_bkg_failed,
              dm_f_extended_list_d0_bkg_passed,
              dm_f_extended_list_d0_bkg_failed;

          RooDataSet ds_dm_d0_bkg_passed_sum(
              "ds_dm_d0_bkg_passed_sum_" + suffix,
              "ds_dm_d0_bkg_passed_sum_" + suffix, fit_vars_w,
              WeightVar(w_d0_bkg));
          RooDataSet ds_dm_d0_bkg_failed_sum(
              "ds_dm_d0_bkg_failed_sum_" + suffix,
              "ds_dm_d0_bkg_failed_sum_" + suffix, fit_vars_w,
              WeightVar(w_d0_bkg));

          unordered_map<string, TH1D> dm_histos_d0_bkg_passed;
          unordered_map<string, TH1D> dm_histos_d0_bkg_failed;
          TH1D                        h_dm_d0_bkg_passed(
              "h_dm_d0_bkg_passed_" + suffix, "h_dm_d0_bkg_passed_" + suffix,
              bins_histos_dm_d0_bkg.numBins(), bins_histos_dm_d0_bkg.lowBound(),
              bins_histos_dm_d0_bkg.highBound());
          TH1D h_dm_d0_bkg_failed(
              "h_dm_d0_bkg_failed_" + suffix, "h_dm_d0_bkg_failed_" + suffix,
              bins_histos_dm_d0_bkg.numBins(), bins_histos_dm_d0_bkg.lowBound(),
              bins_histos_dm_d0_bkg.highBound());
          THStack ths_dm_d0_bkg_passed("ths_dm_d0_bkg_passed_" + suffix,
                                       "ths_dm_d0_bkg_passed_" + suffix);
          THStack ths_dm_d0_bkg_failed("ths_dm_d0_bkg_failed_" + suffix,
                                       "ths_dm_d0_bkg_failed_" + suffix);

          double sum_f_passed_dm = 0., sum_f_failed_dm = 0.,
                 sum_f_ext_passed_dm = 0., sum_f_ext_failed_dm = 0.;

          for (const auto &d0_decay : d0_bkg_decays) {
            cout << "\nINFO Producing dm PDFs for " << d0_decay << endl;
            const TString suffix_bkg = d0_decay + "_" + suffix;

            RooDataSet &ds_passed = datasets_d0_bkg_mc_passed[d0_decay][p_idx];
            RooDataSet &ds_failed = datasets_d0_bkg_mc_failed[d0_decay][p_idx];

            const double n_passed_d0_decay =
                ds_passed.numEntries() * w_d0_decays.at(d0_decay);
            const double n_failed_d0_decay =
                ds_failed.numEntries() * w_d0_decays.at(d0_decay);

            w_d0_bkg.setVal(w_d0_decays.at(d0_decay) * n_passed_d0_decay /
                            (n_passed_d0_decay + n_failed_d0_decay));
            append_to_dataset(ds_passed, ds_dm_d0_bkg_passed_sum, w_d0_bkg);
            append_to_dataset(ds_failed, ds_dm_d0_bkg_passed_sum, w_d0_bkg);

            w_d0_bkg.setVal(w_d0_decays.at(d0_decay) * n_failed_d0_decay /
                            (n_passed_d0_decay + n_failed_d0_decay));
            append_to_dataset(ds_passed, ds_dm_d0_bkg_failed_sum, w_d0_bkg);
            append_to_dataset(ds_failed, ds_dm_d0_bkg_failed_sum, w_d0_bkg);

            RooDataSet ds_dm_d0_bkg("ds_dm_d0_bkg_" + suffix,
                                    "ds_dm_d0_bkg_" + suffix, fit_vars);

            ds_dm_d0_bkg.append(ds_passed);
            ds_dm_d0_bkg.append(ds_failed);

            unique_ptr<TH1D> h_passed(static_cast<TH1D *>(
                ds_passed.createHistogram("h_passed_" + suffix_bkg, dm_var,
                                          Binning(bins_histos_dm_d0_bkg))));
            unique_ptr<TH1D> h_failed(static_cast<TH1D *>(
                ds_failed.createHistogram("h_failed_" + suffix_bkg, dm_var,
                                          Binning(bins_histos_dm_d0_bkg))));

            h_passed->SetFillColor(color_ids_d0_decays.at(d0_decay));
            h_failed->SetFillColor(color_ids_d0_decays.at(d0_decay));

            dm_histos_d0_bkg_passed.emplace(
                d0_decay, TH1D("h_dm_d0_bkg_passed_merged_" + suffix_bkg,
                               "h_dm_d0_bkg_passed_merged_" + suffix_bkg,
                               bins_histos_dm_d0_bkg.numBins(),
                               bins_histos_dm_d0_bkg.lowBound(),
                               bins_histos_dm_d0_bkg.highBound()));
            dm_histos_d0_bkg_failed.emplace(
                d0_decay, TH1D("h_dm_d0_bkg_failed_merged_" + suffix_bkg,
                               "h_dm_d0_bkg_failed_merged_" + suffix_bkg,
                               bins_histos_dm_d0_bkg.numBins(),
                               bins_histos_dm_d0_bkg.lowBound(),
                               bins_histos_dm_d0_bkg.highBound()));

            TH1D &h_passed_merged = dm_histos_d0_bkg_passed[d0_decay];
            TH1D &h_failed_merged = dm_histos_d0_bkg_failed[d0_decay];

            h_passed_merged.Add(h_passed.get(), h_failed.get());
            h_failed_merged.Add(h_passed.get(), h_failed.get());

            h_passed_merged.SetFillColor(color_ids_d0_decays.at(d0_decay));
            h_failed_merged.SetFillColor(color_ids_d0_decays.at(d0_decay));

            h_passed_merged.Scale(n_passed_d0_decay /
                                  h_passed_merged.GetEntries());
            h_failed_merged.Scale(n_failed_d0_decay /
                                  h_failed_merged.GetEntries());
            ths_dm_d0_bkg_passed.Add(&h_passed_merged);
            ths_dm_d0_bkg_failed.Add(&h_failed_merged);
            h_dm_d0_bkg_passed.Add(&h_passed_merged);
            h_dm_d0_bkg_failed.Add(&h_failed_merged);

            cout << "INFO Building RooKeysPdf dm_model_d0_bkg_passed with "
                 << ds_dm_d0_bkg.numEntries() << " entries" << endl;
            dm_pdfs_d0_bkg_passed.emplace(
                d0_decay,
                RooKeysPdf("dm_model_d0_bkg_passed_" + suffix_bkg,
                           "dm_model_d0_bkg_passed_" + suffix_bkg, dm_var,
                           dm_var, dm_pdg, dm_scale, dm_shift, ds_dm_d0_bkg,
                           RooKeysPdf::NoMirror, 1.6, true, n_bins_keys));

            cout << "INFO Building RooKeysPdf dm_model_d0_bkg_failed with "
                 << ds_dm_d0_bkg.numEntries() << " entries" << endl;
            dm_pdfs_d0_bkg_failed.emplace(
                d0_decay,
                RooKeysPdf("dm_model_d0_bkg_failed_" + suffix_bkg,
                           "dm_model_d0_bkg_failed_" + suffix_bkg, dm_var,
                           dm_var, dm_pdg, dm_scale, dm_shift, ds_dm_d0_bkg,
                           RooKeysPdf::NoMirror, 1.6, true, n_bins_keys));

            unique_ptr<RooPlot> frame_dm_mc_d0_bkg_passed(
                dm_var.frame(Title("dm Passed " + tag)));
            unique_ptr<RooPlot> frame_dm_mc_d0_bkg_failed(
                dm_var.frame(Title("dm Failed " + tag)));
            unique_ptr<RooPlot> frame_dm_mc_d0_bkg_total(
                dm_var.frame(Title("dm Total " + tag)));

            ds_passed.plotOn(frame_dm_mc_d0_bkg_passed.get(),
                             Binning(bins_histos_dm_d0_bkg));
            ds_failed.plotOn(frame_dm_mc_d0_bkg_failed.get(),
                             Binning(bins_histos_dm_d0_bkg));
            ds_dm_d0_bkg.plotOn(frame_dm_mc_d0_bkg_total.get(),
                                Binning(bins_histos_dm_d0_bkg));

            dm_pdfs_d0_bkg_passed.at(d0_decay).plotOn(
                frame_dm_mc_d0_bkg_passed.get(), LineWidth(1), LineColor(kRed),
                Range("fitRange"), VLines());
            dm_pdfs_d0_bkg_passed.at(d0_decay).plotOn(
                frame_dm_mc_d0_bkg_passed.get(), LineWidth(2),
                LineColor(color_ids_d0_decays.at(d0_decay)),
                NormRange("fitRange"));
            dm_pdfs_d0_bkg_failed.at(d0_decay).plotOn(
                frame_dm_mc_d0_bkg_failed.get(), LineWidth(1), LineColor(kRed),
                Range("fitRange"), VLines());
            dm_pdfs_d0_bkg_failed.at(d0_decay).plotOn(
                frame_dm_mc_d0_bkg_failed.get(), LineWidth(2),
                LineColor(color_ids_d0_decays.at(d0_decay)),
                NormRange("fitRange"));
            dm_pdfs_d0_bkg_failed.at(d0_decay).plotOn(
                frame_dm_mc_d0_bkg_total.get(), LineWidth(1), LineColor(kRed),
                Range("fitRange"), VLines());
            dm_pdfs_d0_bkg_failed.at(d0_decay).plotOn(
                frame_dm_mc_d0_bkg_total.get(), LineWidth(2),
                LineColor(color_ids_d0_decays.at(d0_decay)),
                NormRange("fitRange"));

            c3.cd(1);
            frame_dm_mc_d0_bkg_passed->Draw();
            c3.cd(2);
            frame_dm_mc_d0_bkg_failed->Draw();
            c3.cd(3);
            frame_dm_mc_d0_bkg_total->Draw();
            c3.SaveAs(fit_dir_path + "/dm_" + suffix + "_d0_bkg_keys_" +
                      d0_decay + ".pdf");

            dm_pdf_list_d0_bkg_passed.add(dm_pdfs_d0_bkg_passed.at(d0_decay));
            dm_pdf_list_d0_bkg_failed.add(dm_pdfs_d0_bkg_failed.at(d0_decay));

            if (d0_decay != d0_bkg_decays.back()) {
              // Here we calculate fit fractions, so the this must be omitted
              // for the last component

              // Fit fractions within fit range
              dm_fs_d0_bkg_passed.emplace(
                  d0_decay,
                  RooRealVar("dm_f_d0_bkg_passed_" + suffix_bkg,
                             "dm_f_d0_bkg_passed_" + suffix_bkg,
                             w_d0_decays.at(d0_decay) *
                                 d0_bkg_counts_passed[probe][p_idx][d0_decay] /
                                 d0_bkg_sum_passed,
                             ""));

              dm_fs_d0_bkg_failed.emplace(
                  d0_decay,
                  RooRealVar("dm_f_d0_bkg_failed_" + suffix_bkg,
                             "dm_f_d0_bkg_failed_" + suffix_bkg,
                             w_d0_decays.at(d0_decay) *
                                 d0_bkg_counts_failed[probe][p_idx][d0_decay] /
                                 d0_bkg_sum_failed,
                             ""));
              dm_f_list_d0_bkg_passed.add(dm_fs_d0_bkg_passed.at(d0_decay));
              dm_f_list_d0_bkg_failed.add(dm_fs_d0_bkg_failed.at(d0_decay));

              // Fit fractions over extended range
              dm_fs_extended_d0_bkg_passed.emplace(
                  d0_decay,
                  RooRealVar("dm_f_extended_d0_bkg_passed_" + suffix_bkg,
                             "dm_f_extended_d0_bkg_passed_" + suffix_bkg,
                             n_passed_d0_decay / d0_bkg_sum_extended_passed,
                             ""));

              dm_fs_extended_d0_bkg_failed.emplace(
                  d0_decay,
                  RooRealVar("dm_f_extended_d0_bkg_failed_" + suffix_bkg,
                             "dm_f_extended_d0_bkg_failed_" + suffix_bkg,
                             n_failed_d0_decay / d0_bkg_sum_extended_failed,
                             ""));
              dm_f_extended_list_d0_bkg_passed.add(
                  dm_fs_extended_d0_bkg_passed.at(d0_decay));
              dm_f_extended_list_d0_bkg_failed.add(
                  dm_fs_extended_d0_bkg_failed.at(d0_decay));

              cout << "INFO " << d0_decay
                   << " f passed: " << dm_fs_d0_bkg_passed.at(d0_decay).getVal()
                   << endl;
              cout << "INFO " << d0_decay
                   << " f failed: " << dm_fs_d0_bkg_failed.at(d0_decay).getVal()
                   << endl;
              cout << "INFO " << d0_decay << " f ext passed: "
                   << dm_fs_extended_d0_bkg_passed.at(d0_decay).getVal()
                   << endl;
              cout << "INFO " << d0_decay << " f ext failed: "
                   << dm_fs_extended_d0_bkg_failed.at(d0_decay).getVal()
                   << endl;

              sum_f_passed_dm += dm_fs_d0_bkg_passed.at(d0_decay).getVal();
              sum_f_failed_dm += dm_fs_d0_bkg_failed.at(d0_decay).getVal();
              sum_f_ext_passed_dm +=
                  dm_fs_extended_d0_bkg_passed.at(d0_decay).getVal();
              sum_f_ext_failed_dm +=
                  dm_fs_extended_d0_bkg_failed.at(d0_decay).getVal();
            }
          }

          cout << "INFO " << d0_bkg_decays.back()
               << " f passed: " << 1. - sum_f_passed_dm << endl;
          cout << "INFO " << d0_bkg_decays.back()
               << " f failed: " << 1. - sum_f_failed_dm << endl;
          cout << "INFO " << d0_bkg_decays.back()
               << " f ext passed: " << 1. - sum_f_ext_passed_dm << endl;
          cout << "INFO " << d0_bkg_decays.back()
               << " f ext failed: " << 1. - sum_f_ext_failed_dm << endl;

          // Plot stacked histograms
          h_dm_d0_bkg_passed.SetLineColor(kBlack);
          h_dm_d0_bkg_failed.SetLineColor(kBlack);

          c_double.cd(1);
          ths_dm_d0_bkg_passed.Draw("B HIST");
          h_dm_d0_bkg_passed.Draw("L SAME");

          c_double.cd(2);
          ths_dm_d0_bkg_failed.Draw("B HIST");
          h_dm_d0_bkg_failed.Draw("L SAME");

          c_double.SaveAs(fit_dir_path + "/dm_" + suffix +
                          "_d0_bkg_histos_stacked.pdf");

          // Build pdfs with fit fractions calculated within the fit range
          // This is the one we use for the fits!
          RooAddPdf dm_model_d0_bkg_passed("dm_model_d0_bkg_passed_" + suffix,
                                           "dm_model_d0_bkg_passed_" + suffix,
                                           dm_pdf_list_d0_bkg_passed,
                                           dm_f_list_d0_bkg_passed);

          RooAddPdf dm_model_d0_bkg_failed("dm_model_d0_bkg_failed_" + suffix,
                                           "dm_model_d0_bkg_failed_" + suffix,
                                           dm_pdf_list_d0_bkg_failed,
                                           dm_f_list_d0_bkg_failed);

          // Build pdfs fit fractions calculated within the extended range
          RooAddPdf dm_model_d0_bkg_passed_ext(
              "dm_model_d0_bkg_passed_ext_" + suffix,
              "dm_model_d0_bkg_passed_ext_" + suffix, dm_pdf_list_d0_bkg_passed,
              dm_f_extended_list_d0_bkg_passed);

          RooAddPdf dm_model_d0_bkg_failed_ext(
              "dm_model_d0_bkg_failed_ext_" + suffix,
              "dm_model_d0_bkg_failed_ext_" + suffix, dm_pdf_list_d0_bkg_failed,
              dm_f_extended_list_d0_bkg_failed);

          // Plot PDF over extended range
          unique_ptr<RooPlot> frame_dm_mc_d0_bkg_passed_ext(
              dm_var.frame(Title("dm Passed " + tag)));
          unique_ptr<RooPlot> frame_dm_mc_d0_bkg_failed_ext(
              dm_var.frame(Title("dm Failed " + tag)));
          ds_dm_d0_bkg_passed_sum.plotOn(frame_dm_mc_d0_bkg_passed_ext.get(),
                                         Binning(bins_histos_dm_d0_bkg));
          ds_dm_d0_bkg_failed_sum.plotOn(frame_dm_mc_d0_bkg_failed_ext.get(),
                                         Binning(bins_histos_dm_d0_bkg));

          dm_model_d0_bkg_passed_ext.plotOn(frame_dm_mc_d0_bkg_passed_ext.get(),
                                            LineWidth(1), LineColor(kRed),
                                            Range("fitRange"), VLines());
          dm_model_d0_bkg_passed_ext.plotOn(frame_dm_mc_d0_bkg_passed_ext.get(),
                                            LineWidth(2),
                                            NormRange("fitRange"));
          dm_model_d0_bkg_failed_ext.plotOn(frame_dm_mc_d0_bkg_failed_ext.get(),
                                            LineWidth(1), LineColor(kRed),
                                            Range("fitRange"), VLines());
          dm_model_d0_bkg_failed_ext.plotOn(frame_dm_mc_d0_bkg_failed_ext.get(),
                                            LineWidth(2),
                                            NormRange("fitRange"));

          c_double.cd(1);
          frame_dm_mc_d0_bkg_passed_ext->Draw();
          c_double.cd(2);
          frame_dm_mc_d0_bkg_failed_ext->Draw();
          c_double.SaveAs(fit_dir_path + "/dm_" + suffix + "_d0_bkg_ext.pdf");

          // Plot PDF in fit range
          d0_m_var.setRange(d0m_range_min, d0m_range_max);
          dm_var.setRange(dm_range_min, dm_range_max);

          unique_ptr<RooPlot> frame_dm_mc_d0_bkg_passed(
              dm_var.frame(Title("dm Passed " + tag), Range("fitRange")));
          unique_ptr<RooPlot> frame_dm_mc_d0_bkg_failed(
              dm_var.frame(Title("dm Failed " + tag), Range("fitRange")));
          ds_dm_d0_bkg_passed_sum.plotOn(frame_dm_mc_d0_bkg_passed.get(),
                                         Binning(bins_histos_dm_d0_bkg),
                                         CutRange("fitRange"));
          ds_dm_d0_bkg_failed_sum.plotOn(frame_dm_mc_d0_bkg_failed.get(),
                                         Binning(bins_histos_dm_d0_bkg),
                                         CutRange("fitRange"));

          const double d0_bkg_norm_passed_dm =
              ds_dm_d0_bkg_passed_sum.sumEntries("", "fitRange");
          const double d0_bkg_norm_failed_dm =
              ds_dm_d0_bkg_failed_sum.sumEntries("", "fitRange");
          dm_model_d0_bkg_passed_ext.plotOn(
              frame_dm_mc_d0_bkg_passed.get(), LineWidth(2), Range("fitRange"),
              Normalization(d0_bkg_norm_passed_dm, RooAbsReal::NumEvent),
              ProjectionRange("fitRange"));
          dm_model_d0_bkg_failed_ext.plotOn(
              frame_dm_mc_d0_bkg_failed.get(), LineWidth(2), Range("fitRange"),
              Normalization(d0_bkg_norm_failed_dm, RooAbsReal::NumEvent),
              ProjectionRange("fitRange"));

          c_double.cd(1);
          frame_dm_mc_d0_bkg_passed->Draw();
          c_double.cd(2);
          frame_dm_mc_d0_bkg_failed->Draw();
          c_double.SaveAs(fit_dir_path + "/dm_" + suffix + "_d0_bkg.pdf");

          d0_m_var.setRange(extended_d0_m_min, extended_d0_m_max);
          dm_var.setRange(extended_dm_min, extended_dm_max);

          //////////////////////
          // Build full model //
          //////////////////////

          RooAddPdf d0_model_passed(
              "d0_model_passed_" + suffix, "d0_model_passed_" + suffix,
              d0_model_passed_dif, d0_model_passed_nondif, f_dif_calib_passed);

          RooAddPdf dm_model_passed(
              "dm_model_passed_" + suffix, "dm_model_passed_" + suffix,
              dm_model_passed_dif, dm_model_passed_nondif, f_dif_calib_passed);

          RooAddPdf d0_model_failed(
              "d0_model_failed_" + suffix, "d0_model_failed_" + suffix,
              d0_model_failed_dif, d0_model_failed_nondif, f_dif_calib_failed);

          RooAddPdf dm_model_failed(
              "dm_model_failed_" + suffix, "dm_model_failed_" + suffix,
              dm_model_failed_dif, dm_model_failed_nondif, f_dif_calib_failed);

          RooProdPdf sig_passed("sig_passed_" + suffix, "sig_passed_" + suffix,
                                d0_model_passed, dm_model_passed);
          RooProdPdf sig_failed("sig_failed_" + suffix, "sig_failed_" + suffix,
                                d0_model_failed, dm_model_failed);

          RooProdPdf comb_spi_passed("comb_spi_passed_" + suffix,
                                     "comb_spi_passed_" + suffix,
                                     d0_model_passed, dm_comb_spi_passed);
          RooProdPdf comb_spi_failed("comb_spi_failed_" + suffix,
                                     "comb_spi_failed_" + suffix,
                                     d0_model_failed, dm_comb_spi_failed);

          RooProdPdf d0_bkg_passed(
              "d0_bkg_passed_" + suffix, "d0_bkg_passed_" + suffix,
              d0_model_d0_bkg_passed, dm_model_d0_bkg_passed);
          RooProdPdf d0_bkg_failed(
              "d0_bkg_failed_" + suffix, "d0_bkg_failed_" + suffix,
              d0_model_d0_bkg_failed, dm_model_d0_bkg_failed);

          RooProdPdf comb_all_passed("comb_all_passed_" + suffix,
                                     "comb_all_passed_" + suffix,
                                     d0_comb_all_passed, dm_comb_all_passed);
          RooProdPdf comb_all_failed("comb_all_failed_" + suffix,
                                     "comb_all_failed_" + suffix,
                                     d0_comb_all_failed, dm_comb_all_failed);

          RooAddPdf model_passed("model_passed_" + suffix,
                                 "model_passed_" + suffix,
                                 RooArgList(sig_passed, d0_bkg_passed,
                                            comb_spi_passed, comb_all_passed),
                                 RooArgList(f_sig_phys_passed, f_d0_bkg_passed,
                                            f_comb_spi_passed));

          RooAddPdf model_failed("model_failed_" + suffix,
                                 "model_failed_" + suffix,
                                 RooArgList(sig_failed, d0_bkg_failed,
                                            comb_spi_failed, comb_all_failed),
                                 RooArgList(f_sig_phys_failed, f_d0_bkg_failed,
                                            f_comb_spi_failed));

          ////////////////////
          // Set parameters //
          ////////////////////

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

          // Define approximate signal region
          const double d0m_sig_max = d0m_range_max - 0.015;  // GeV
          const double d0m_sig_min = d0m_range_min + 0.015;  // GeV
          const double dm_sig_max  = dm_range_max - 0.002;   // GeV
          const double dm_sig_min  = dm_range_min + 0.002;   // GeV

          // Define ranges for each component
          // Signal
          d0_m_var.setRange("sig", d0m_sig_min, d0m_sig_max);
          dm_var.setRange("sig", dm_sig_min, dm_sig_max);
          // SB1: soft pi comb (+ all comb)
          d0_m_var.setRange("sb1_1", d0m_sig_min, d0m_sig_max);
          d0_m_var.setRange("sb1_2", d0m_sig_min, d0m_sig_max);
          dm_var.setRange("sb1_1", dm_range_min, dm_sig_min);
          dm_var.setRange("sb1_2", dm_sig_max, dm_range_max);
          // SB2: part. reco. D0 comb (+ all comb)
          d0_m_var.setRange("sb2_1", d0m_range_min, d0m_sig_min);
          d0_m_var.setRange("sb2_2", d0m_sig_max, d0m_range_max);
          dm_var.setRange("sb2_1", dm_sig_min, dm_sig_max);
          dm_var.setRange("sb2_2", dm_sig_min, dm_sig_max);
          // SB3: all comb
          d0_m_var.setRange("sb3_1", d0m_range_min, d0m_sig_min);
          d0_m_var.setRange("sb3_2", d0m_range_min, d0m_sig_min);
          d0_m_var.setRange("sb3_3", d0m_sig_max, d0m_range_max);
          d0_m_var.setRange("sb3_4", d0m_sig_max, d0m_range_max);
          dm_var.setRange("sb3_1", dm_range_min, dm_sig_min);
          dm_var.setRange("sb3_2", dm_sig_max, dm_range_max);
          dm_var.setRange("sb3_3", dm_range_min, dm_sig_min);
          dm_var.setRange("sb3_4", dm_sig_max, dm_range_max);

          const double A_fit =
              (d0m_range_max - d0m_range_min) * (dm_range_max - dm_range_min);
          const double A_sig =
              (d0m_sig_max - d0m_sig_min) * (dm_sig_max - dm_sig_min);
          const double A_sb3 =
              (dm_range_max - dm_sig_max) * (d0m_sig_min - d0m_range_min) +
              (dm_range_max - dm_sig_max) * (d0m_range_max - d0m_sig_max) +
              (dm_sig_min - dm_range_min) * (d0m_sig_min - d0m_range_min) +
              (dm_sig_min - dm_range_min) * (d0m_range_max - d0m_sig_max);
          // const double A_sb2 = (dm_sig_max - dm_sig_min) * (d0m_sig_min -
          // d0m_range_min) +
          //                          (dm_sig_max - dm_sig_min) *
          //                          (d0m_range_max - d0m_sig_max);
          const double A_sb1 =
              (dm_range_max - dm_sig_max) * (d0m_sig_max - d0m_sig_min) +
              (dm_sig_min - dm_range_min) * (d0m_sig_max - d0m_sig_min);

          // Yields in each SB
          const double n_sb1_calib_passed =
              dataset_calib_passed->sumEntries("", "sb1_1,sb1_2");
          const double n_sb1_calib_failed =
              dataset_calib_failed->sumEntries("", "sb1_1,sb1_2");
          // const double n_sb2_calib_passed =
          //     dataset_calib_passed->sumEntries(sb2_cut);
          // const double n_sb2_calib_failed =
          //     dataset_calib_failed->sumEntries(sb2_cut);
          const double n_sb3_calib_passed =
              dataset_calib_passed->sumEntries("", "sb3_1,sb3_2,sb3_3,sb3_4");
          const double n_sb3_calib_failed =
              dataset_calib_failed->sumEntries("", "sb3_1,sb3_2,sb3_3,sb3_4");

          // Estimated yields of each bkg component in full fit range
          const double n_bkg1_calib_passed =
              (n_sb1_calib_passed - n_sb3_calib_passed * A_sb1 / A_sb3) *
              (A_sb1 + A_sig) / A_sb1;
          const double n_bkg1_calib_failed =
              (n_sb1_calib_failed - n_sb3_calib_failed * A_sb1 / A_sb3) *
              (A_sb1 + A_sig) / A_sb1;
          // const double n_bkg2_calib_passed =
          //     (n_sb2_calib_passed - n_sb3_calib_passed * A_sb2 / A_sb3) *
          //     (A_sb2 + A_sig) / A_sb2;
          // const double n_bkg2_calib_failed =
          //     (n_sb2_calib_failed - n_sb3_calib_failed * A_sb2 / A_sb3) *
          //     (A_sb2 + A_sig) / A_sb2;
          const double n_bkg3_calib_passed = n_sb3_calib_passed * A_fit / A_sb3;
          const double n_bkg3_calib_failed = n_sb3_calib_failed * A_fit / A_sb3;

          const double n_comb_guess_passed =
              n_bkg1_calib_passed + n_bkg3_calib_passed;
          const double n_comb_guess_failed =
              n_bkg1_calib_failed + n_bkg3_calib_failed;

          const double f_phys_passed_guess =
              1.0 - n_comb_guess_passed / n_calib_passed;
          const double f_phys_failed_guess =
              1.0 - n_comb_guess_failed / n_calib_failed;

          f_phys_passed.setVal(f_phys_passed_guess);
          f_phys_failed.setVal(f_phys_failed_guess);

          const double f_spi_passed_guess =
              n_bkg1_calib_passed / n_comb_guess_passed;
          const double f_spi_failed_guess =
              n_bkg1_calib_failed / n_comb_guess_failed;

          f_spi_passed.setVal(f_spi_passed_guess);
          f_spi_failed.setVal(f_spi_failed_guess);

          // Constrain yield of D0 bkg based on MC
          const double d0_bkg_passed_signal =
              ymlBkg[sample][probe]["passed"][p_idx]["signal"].as<double>();
          const double d0_bkg_passed_total =
              ymlBkg[sample][probe]["passed"][p_idx]["total"].as<double>();
          const double d0_bkg_failed_signal =
              ymlBkg[sample][probe]["failed"][p_idx]["signal"].as<double>();
          const double d0_bkg_failed_total =
              ymlBkg[sample][probe]["failed"][p_idx]["total"].as<double>();

          // Additional width on constraints to account for MC-data
          // discrepancies and outdated BFs in decay.dec
          constexpr double d0_bkg_syst_unc = 0.02;

          const double f_sig_passed_guess =
              d0_bkg_passed_signal / d0_bkg_passed_total;
          const double f_sig_passed_guess_unc_hi =
              TEfficiency::Wilson(d0_bkg_passed_total, d0_bkg_passed_signal,
                                  ONE_SIGMA, true) -
              f_sig_passed_guess;
          const double f_sig_passed_guess_unc_lo =
              f_sig_passed_guess - TEfficiency::Wilson(d0_bkg_passed_total,
                                                       d0_bkg_passed_signal,
                                                       ONE_SIGMA, false);

          const double f_sig_passed_guess_unc_total_hi =
              sqrt_sum_sq(f_sig_passed_guess_unc_hi, d0_bkg_syst_unc);
          const double f_sig_passed_guess_unc_total_lo =
              sqrt_sum_sq(f_sig_passed_guess_unc_lo, d0_bkg_syst_unc);

          f_sig_passed_mean.setVal(f_sig_passed_guess);
          f_sig_passed_unc_hi.setVal(std::min(f_sig_passed_guess_unc_total_hi,
                                              1.0 - f_sig_passed_guess));
          f_sig_passed_unc_lo.setVal(
              std::min(f_sig_passed_guess_unc_total_lo, f_sig_passed_guess));

          const double f_sig_failed_guess =
              d0_bkg_failed_signal / d0_bkg_failed_total;
          const double f_sig_failed_guess_unc_hi =
              TEfficiency::Wilson(d0_bkg_failed_total, d0_bkg_failed_signal,
                                  ONE_SIGMA, true) -
              f_sig_failed_guess;
          const double f_sig_failed_guess_unc_lo =
              f_sig_failed_guess - TEfficiency::Wilson(d0_bkg_failed_total,
                                                       d0_bkg_failed_signal,
                                                       ONE_SIGMA, false);

          const double f_sig_failed_guess_unc_total_hi =
              sqrt_sum_sq(f_sig_failed_guess_unc_hi, d0_bkg_syst_unc);
          const double f_sig_failed_guess_unc_total_lo =
              sqrt_sum_sq(f_sig_failed_guess_unc_lo, d0_bkg_syst_unc);

          f_sig_failed_mean.setVal(f_sig_failed_guess);
          f_sig_failed_unc_hi.setVal(std::min(f_sig_failed_guess_unc_total_hi,
                                              1.0 - f_sig_failed_guess));
          f_sig_failed_unc_lo.setVal(
              std::min(f_sig_failed_guess_unc_total_lo, f_sig_failed_guess));

          f_sig_passed.setVal(f_sig_passed_guess);
          f_sig_failed.setVal(f_sig_failed_guess);

          cout << "\nINFO Constrained ratios of D0 bkg to signal:\n";
          cout << "  - f_sig_passed: " << f_sig_passed_guess << " + "
               << f_sig_passed_unc_hi.getVal() << " - "
               << f_sig_passed_unc_lo.getVal() << "\n";
          cout << "  - f_sig_failed: " << f_sig_failed_guess << " + "
               << f_sig_failed_unc_hi.getVal() << " - "
               << f_sig_failed_unc_lo.getVal() << endl;

          // Set known normalizations
          n_passed.setVal(n_calib_passed);
          n_failed.setVal(n_calib_failed);

          // Estimate exponential coefficients
          constexpr double dx = (1.900 + 1.890) / 2. - (1.835 + 1.825) / 2.;

          const double y1_comb_all_passed = dataset_calib_passed->sumEntries(
              "((dm_var > 0.149) || (dm_var < 0.143)) && (d0_m_var > "
              "1.825) && (d0_m_var < 1.835)");
          const double y1_comb_all_failed = dataset_calib_failed->sumEntries(
              "((dm_var > 0.149) || (dm_var < 0.143)) && (d0_m_var > "
              "1.825) && (d0_m_var < 1.835)");
          const double y2_comb_all_passed = dataset_calib_passed->sumEntries(
              "((dm_var > 0.149) || (dm_var < 0.143)) && (d0_m_var > "
              "1.890) && (d0_m_var < 1.900)");
          const double y2_comb_all_failed = dataset_calib_failed->sumEntries(
              "((dm_var > 0.149) || (dm_var < 0.143)) && (d0_m_var > "
              "1.890) && (d0_m_var < 1.900)");

          const double k_comb_all_failed_guess =
              (y2_comb_all_failed > 0.) && (y1_comb_all_failed > 0.)
                  ? log(y2_comb_all_failed / y1_comb_all_failed) / dx
                  : 1.;
          const double k_comb_all_passed_guess =
              (y2_comb_all_passed > 0.) && (y1_comb_all_passed > 0.)
                  ? log(y2_comb_all_passed / y1_comb_all_passed) / dx
                  : 1.;

          if (k_comb_all_failed_guess < 0.) {
            k_comb_all_failed.setVal(k_comb_all_failed_guess);
          } else {
            k_comb_all_failed.setVal(-1e-4);
          }

          if (k_comb_all_passed_guess < 0.) {
            k_comb_all_passed.setVal(k_comb_all_passed_guess);
          } else {
            k_comb_all_passed.setVal(-1e-4);
          }

          /////////////
          // Run fit //
          /////////////

          // As of 6.32.16, RooSimultaneous is broken for ranged fits
          // https://github.com/root-project/root/issues/18718

          cout << "\nINFO Starting calib fit " << suffix << endl;

          d0_m_var.setRange(d0m_range_min, d0m_range_max);
          dm_var.setRange(dm_range_min, dm_range_max);
          set_parameters_calib(probe);

          // Create RooFitResult pointer to track lastest data fit
          unique_ptr<RooFitResult> fit_result(nullptr);

          if (!dry_run) {
            auto start_calib = high_resolution_clock::now();

            // Prefit failed subsample floating all parameters except dif
            // fraction
            cout << "\nINFO Prefit failed calib sample " << suffix << endl;

            scale_nondif_calib.setConstant();
            d0_m_shift.setConstant(false);
            d0_m_scale.setConstant(false);
            dm_shift.setConstant(false);
            dm_scale.setConstant(false);

            unique_ptr<RooAbsReal> nll_failed(model_failed.createNLL(
                *dataset_calib_failed, Range("fitRange"),
                ExternalConstraints(
                    RooArgSet(f_sig_failed_constrain, scale_nondif_constrain)),
                GlobalObservables(
                    RooArgSet(f_sig_failed, scale_nondif_calib))));
            // Default nll name causes issues when combining the nlls below
            nll_failed->SetName("nll_failed_" + suffix);

            // In recent root releases (6.30+?), RooAbdPdf::fitTo() and
            // RooAbdPdf::createNLL() accidentaly ignore some flags, such as
            // Offset(). Setting it with RooMinizer::setOffsetting() works
            // though (flags are ignored at NLL construction, but RooMinimizer
            // back-propagates them).

            RooMinimizer minuit_failed(*nll_failed);
            minuit_failed.setPrintLevel(0);
            minuit_failed.setStrategy(2);
            minuit_failed.setOffsetting(true);
            minuit_failed.optimizeConst(true);
            minuit_failed.setMinimizerType("Minuit");

            minuit_failed.migrad();

            // Now prefit passed subsample also fixing d0_m and dm shifts and
            // scalings
            cout << "\nINFO Prefit passed calib sample " << suffix << endl;

            scale_nondif_calib.setConstant();
            d0_m_shift.setConstant();
            d0_m_scale.setConstant();
            dm_shift.setConstant();
            dm_scale.setConstant();

            unique_ptr<RooAbsReal> nll_passed(model_passed.createNLL(
                *dataset_calib_passed, Range("fitRange"),
                ExternalConstraints(
                    RooArgSet(f_sig_passed_constrain, scale_nondif_constrain)),
                GlobalObservables(
                    RooArgSet(f_sig_passed, scale_nondif_calib))));
            nll_passed->SetName("nll_passed_" + suffix);

            RooMinimizer minuit_passed(*nll_passed);
            minuit_passed.setPrintLevel(0);
            minuit_passed.setStrategy(2);
            minuit_passed.setOffsetting(true);
            minuit_passed.optimizeConst(true);
            minuit_passed.setMinimizerType("Minuit");

            minuit_passed.migrad();

            // Now build combined NLL and perform full fit
            cout << "\nINFO Fit full calib sample " << suffix << endl;

            RooAddition nll_combined("nll_combined", "nll_combined",
                                     RooArgSet(*nll_passed, *nll_failed));

            if (float_dif) scale_nondif_calib.setConstant(false);
            d0_m_shift.setConstant(false);
            d0_m_scale.setConstant(false);
            dm_shift.setConstant(false);
            dm_scale.setConstant(false);

            cout << "INFO Initial nondif_yield_add_passed value: "
                 << nondif_yield_add_passed.getVal() << endl;

            RooMinimizer minuit(nll_combined);

            minuit.setPrintLevel(0);
            minuit.setStrategy(2);
            minuit.setOffsetting(true);
            minuit.optimizeConst(true);
            minuit.setMinimizerType("Minuit");

            minuit.migrad();
            int fit_status = minuit.save()->status();

            int cov_matrix_status = -1;
            if (fit_status == 0 || fit_status == 4000) {
              minuit.hesse();
              cov_matrix_status = minuit.save()->covQual();
            }

            int minos_status = -1;
            if (use_minos && cov_matrix_status == 3) {
              minos_status = minuit.minos();
            }

            cout << "\nINFO Fit status: " << fit_status << "\n";
            cout << "INFO Covariance matrix status: " << cov_matrix_status
                 << "\n";
            cout << "INFO MINOS status: " << minos_status << "\n" << endl;

            int calib_fit_reattempts = 0;
            while ((cov_matrix_status != 3) &&
                   (calib_fit_reattempts < max_fix_reattempts)) {
              calib_fit_reattempts++;
              cov_matrix_status = -1;
              minos_status      = -1;
              cout << "INFO Retrying calib fit " << suffix << " (retry #"
                   << calib_fit_reattempts << ")" << endl;

              // Call SIMPLEX
              minuit.simplex();

              // Call MIGRAD
              minuit.migrad();
              fit_status = minuit.save()->status();

              if (fit_status == 0 || fit_status == 4000) {
                minuit.hesse();
                cov_matrix_status = minuit.save()->covQual();
              }

              if (use_minos && cov_matrix_status == 3) {
                minos_status = minuit.minos();
              }

              cout << "\nINFO Fit status: " << fit_status << "\n";
              cout << "INFO Covariance matrix status: " << cov_matrix_status
                   << "\n";
              cout << "INFO MINOS status: " << minos_status << "\n" << endl;
            }

            fit_result.reset(minuit.save());

            auto stop_calib = high_resolution_clock::now();

            const int duration_calib =
                duration_cast<seconds>(stop_calib - start_calib).count();

            cout << "INFO Calib fits took " << format_time(duration_calib)
                 << "\n";
            cout << "INFO Needed " << calib_fit_reattempts << " retries with "
                 << format_time(duration_calib / (calib_fit_reattempts + 1))
                 << " on average\n"
                 << endl;

            calib_retries.Fill(calib_fit_reattempts);

            if (cov_matrix_status == 3) {
              // Monitor how much f_sig deviates from constraint
              double f_sig_passed_sigma =
                  f_sig_passed.getVal() - f_sig_passed_guess;
              if (f_sig_passed_sigma >= 0.) {
                f_sig_passed_sigma /= f_sig_passed_unc_hi.getVal();
              } else {
                f_sig_passed_sigma /= f_sig_passed_unc_lo.getVal();
              }

              double f_sig_failed_sigma =
                  f_sig_failed.getVal() - f_sig_failed_guess;
              if (f_sig_failed_sigma >= 0.) {
                f_sig_failed_sigma /= f_sig_failed_unc_hi.getVal();
              } else {
                f_sig_failed_sigma /= f_sig_failed_unc_lo.getVal();
              }

              f_sig_passed_pull.setVal(f_sig_passed_sigma);
              f_sig_failed_pull.setVal(f_sig_failed_sigma);

              ds_params_calib.addFast(params_calib);
            } else {
              f_sig_passed_pull.setVal(-100.);
              f_sig_failed_pull.setVal(-100.);
            }

            fit_status_calib.SetBinContent(kin_bin, fit_status);
            fit_cov_qual_calib.SetBinContent(kin_bin, cov_matrix_status);
          }

          //////////////////////
          // Plot fit results //
          //////////////////////

          // Build modified binning for data to match shifted/scaled PDFs
          RooBinning bins_histos_d0_m_failed_mod(
              nbins_failed + 2 * extra_bins_failed,
              d0_m_scale.getVal() * (extended_d0_m_min - D0_M) + D0_M +
                  d0_m_shift.getVal(),
              d0_m_scale.getVal() * (extended_d0_m_max - D0_M) + D0_M +
                  d0_m_shift.getVal(),
              "bins_histos_d0_m_failed_mod");
          RooBinning bins_histos_dm_failed_mod(
              nbins_failed + 2 * extra_bins_failed,
              dm_scale.getVal() * (extended_dm_min - DM) + DM +
                  dm_shift.getVal(),
              dm_scale.getVal() * (extended_dm_max - DM) + DM +
                  dm_shift.getVal(),
              "bins_histos_dm_failed_mod");
          RooBinning bins_histos_d0_m_passed_mod(
              nbins_passed + 2 * extra_bins_passed,
              d0_m_scale.getVal() * (extended_d0_m_min - D0_M) + D0_M +
                  d0_m_shift.getVal(),
              d0_m_scale.getVal() * (extended_d0_m_max - D0_M) + D0_M +
                  d0_m_shift.getVal(),
              "bins_histos_d0_m_passed_mod");
          RooBinning bins_histos_dm_passed_mod(
              nbins_passed + 2 * extra_bins_passed,
              dm_scale.getVal() * (extended_dm_min - DM) + DM +
                  dm_shift.getVal(),
              dm_scale.getVal() * (extended_dm_max - DM) + DM +
                  dm_shift.getVal(),
              "bins_histos_dm_passed_mod");

          // Define lambda to plot fit results in different ranges
          auto plot_fit_results = [&](const TString range,
                                      const TString name_suffix) {
            unique_ptr<RooPlot> frame_d0_calib_passed(d0_m_var.frame(
                Title("D0 M Calib Passed " + tag), Range("data")));
            unique_ptr<RooPlot> frame_d0_calib_failed(d0_m_var.frame(
                Title("D0 M Calib Failed " + tag), Range("data")));
            unique_ptr<RooPlot> frame_dm_calib_passed(
                dm_var.frame(Title("dm Calib Passed " + tag), Range("data")));
            unique_ptr<RooPlot> frame_dm_calib_failed(
                dm_var.frame(Title("dm Calib Failed " + tag), Range("data")));

            // TODO:
            // https://root-forum.cern.ch/t/simultaneous-fit-normalization-issue/33965
            if (range == "fitRange") {
              // For plots over full range, also show high d0_M region in
              // low-p pion fits
              dataset_calib_passed->plotOn(frame_d0_calib_passed.get(),
                                           Binning(bins_histos_d0_m_passed));
              dataset_calib_failed->plotOn(frame_d0_calib_failed.get(),
                                           Binning(bins_histos_d0_m_failed));
            } else {
              dataset_calib_passed->plotOn(frame_d0_calib_passed.get(),
                                           Binning(bins_histos_d0_m_passed),
                                           CutRange(range));
              dataset_calib_failed->plotOn(frame_d0_calib_failed.get(),
                                           Binning(bins_histos_d0_m_failed),
                                           CutRange(range));
            }
            dataset_calib_passed->plotOn(frame_dm_calib_passed.get(),
                                         Binning(bins_histos_dm_passed),
                                         CutRange(range));
            dataset_calib_failed->plotOn(frame_dm_calib_failed.get(),
                                         Binning(bins_histos_dm_failed),
                                         CutRange(range));

            model_passed.plotOn(frame_d0_calib_passed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2));
            model_passed.plotOn(frame_d0_calib_passed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2), Components(sig_passed),
                                LineColor(kTeal + 2));
            model_passed.plotOn(frame_d0_calib_passed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2), Components(d0_bkg_passed),
                                LineColor(kViolet));
            model_passed.plotOn(frame_d0_calib_passed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2), Components(comb_spi_passed),
                                LineColor(kOrange));
            model_passed.plotOn(frame_d0_calib_passed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2), Components(comb_all_passed),
                                LineColor(kRed));

            model_failed.plotOn(frame_d0_calib_failed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2));
            model_failed.plotOn(frame_d0_calib_failed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2), Components(sig_failed),
                                LineColor(kTeal + 2));
            model_failed.plotOn(frame_d0_calib_failed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2), Components(d0_bkg_failed),
                                LineColor(kViolet));
            model_failed.plotOn(frame_d0_calib_failed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2), Components(comb_spi_failed),
                                LineColor(kOrange));
            model_failed.plotOn(frame_d0_calib_failed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2), Components(comb_all_failed),
                                LineColor(kRed));

            model_passed.plotOn(frame_dm_calib_passed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2));
            model_passed.plotOn(frame_dm_calib_passed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2), Components(sig_passed),
                                LineColor(kTeal + 2));
            model_passed.plotOn(frame_dm_calib_passed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2), Components(d0_bkg_passed),
                                LineColor(kViolet));
            model_passed.plotOn(frame_dm_calib_passed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2), Components(comb_spi_passed),
                                LineColor(kOrange));
            model_passed.plotOn(frame_dm_calib_passed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2), Components(comb_all_passed),
                                LineColor(kRed));

            model_failed.plotOn(frame_dm_calib_failed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2));
            model_failed.plotOn(frame_dm_calib_failed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2), Components(sig_failed),
                                LineColor(kTeal + 2));
            model_failed.plotOn(frame_dm_calib_failed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2), Components(d0_bkg_failed),
                                LineColor(kViolet));
            model_failed.plotOn(frame_dm_calib_failed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2), Components(comb_spi_failed),
                                LineColor(kOrange));
            model_failed.plotOn(frame_dm_calib_failed.get(), Range(range),
                                ProjectionRange(range), NormRange("fitRange"),
                                LineWidth(2), Components(comb_all_failed),
                                LineColor(kRed));

            c_four.cd(1);
            frame_d0_calib_passed->Draw();
            c_four.cd(2);
            frame_dm_calib_passed->Draw();
            c_four.cd(3);
            frame_d0_calib_failed->Draw();
            c_four.cd(4);
            frame_dm_calib_failed->Draw();

            c_four.SaveAs(fit_dir_path + "/fit_calib_" + suffix + name_suffix +
                          ".pdf");
          };

          cout << "\nINFO Plotting full fit region" << endl;
          plot_fit_results("fitRange", "");

          cout << "\nINFO Plotting signal region" << endl;
          plot_fit_results("sig", "_sig");

          cout << "\nINFO Plotting soft pion bkg regions" << endl;
          plot_fit_results("sb1_1", "_sb1_1");
          plot_fit_results("sb1_2", "_sb1_2");

          cout << "\nINFO Plotting d0 bkg regions" << endl;
          plot_fit_results("sb2_1", "_sb2_1");
          plot_fit_results("sb2_2", "_sb2_2");

          cout << "\nINFO Plotting all comb regions" << endl;
          plot_fit_results("sb3_1", "_sb3_1");
          plot_fit_results("sb3_2", "_sb3_2");
          plot_fit_results("sb3_3", "_sb3_3");
          plot_fit_results("sb3_4", "_sb3_4");

          // Save fitted PID efficiency
          if (fake_mu) {
            // For FAKE_MU, we fit calculate the complementary efficiency so
            // that the "passed" sample always corresponds to the K/pi misid
            // case. We recover the proper FAKE_MU PIF efficiency by storing
            // 1.0 - eff instead of eff itself.
            histo_pid_raw.SetBinContent(kin_bin, 1.0 - eff.getVal());
            histo_pid.SetBinContent(kin_bin, 1.0 - eff_corrected.getVal());
          } else {
            histo_pid_raw.SetBinContent(kin_bin, eff.getVal());
            histo_pid.SetBinContent(kin_bin, eff_corrected.getVal());
          }
          // Save fitted DiF fractions
          histo_f_dif.SetBinContent(kin_bin, f_dif_calib_passed.getVal());
          // Save errors if fit was performed
          if (fit_result.get()) {
            histo_pid_raw.SetBinError(kin_bin,
                                      eff.getPropagatedError(*fit_result));
            histo_pid.SetBinError(
                kin_bin, eff_corrected.getPropagatedError(*fit_result));
            histo_f_dif.SetBinError(
                kin_bin, f_dif_calib_passed.getPropagatedError(*fit_result));
          }

          // Store final mass-window efficiencies in histogram
          effs_mw_passed.SetBinContent(kin_bin, eff_mw_passed.getVal());
          effs_mw_failed.SetBinContent(kin_bin, eff_mw_failed.getVal());

          // Print mass-window efficiencies
          cout << "\nINFO Mass-window efficiencies\n";
          cout << " - Passed: " << eff_mw_passed.getVal() * 100. << "%\n";
          cout << " - Failed: " << eff_mw_failed.getVal() * 100. << "%\n";

          // Check normalization estimates
          cout << "\nINFO Fitted vs estimated:\n";
          cout << " - Fitted n_passed = " << n_passed.getVal()
               << " vs estimated " << n_calib_passed << " ("
               << n_passed.getVal() / n_calib_passed << ")\n";
          cout << " - Fitted n_failed = " << n_failed.getVal()
               << " vs estimated " << n_calib_failed << " ("
               << n_failed.getVal() / n_calib_failed << ")\n";

          cout << " - Fitted f_phys_passed = " << f_phys_passed.getVal()
               << " vs estimated " << f_phys_passed_guess << " ("
               << f_phys_passed.getVal() / f_phys_passed_guess << ")\n";
          cout << " - Fitted f_phys_failed = " << f_phys_failed.getVal()
               << " vs estimated " << f_phys_failed_guess << " ("
               << f_phys_failed.getVal() / f_phys_failed_guess << ")\n";

          cout << " - Fitted f_sig_passed = " << f_sig_passed.getVal()
               << " vs estimated " << f_sig_passed_guess << " ("
               << f_sig_passed.getVal() / f_sig_passed_guess << "), "
               << f_sig_passed_pull.getVal() << " sigma away from constrain\n";
          cout << " - Fitted f_sig_failed = " << f_sig_failed.getVal()
               << " vs estimated " << f_sig_failed_guess << " ("
               << f_sig_failed.getVal() / f_sig_failed_guess << "), "
               << f_sig_failed_pull.getVal() << " sigma away from constrain\n";

          cout << " - Fitted f_spi_passed = " << f_spi_passed.getVal()
               << " vs estimated " << f_spi_passed_guess << " ("
               << f_spi_passed.getVal() / f_spi_passed_guess << ")\n";
          cout << " - Fitted f_spi_failed = " << f_spi_failed.getVal()
               << " vs estimated " << f_spi_failed_guess << " ("
               << f_spi_failed.getVal() / f_spi_failed_guess << ")\n";

          cout << " - Fitted k_comb_all_failed = " << k_comb_all_failed.getVal()
               << " vs estimated " << k_comb_all_failed_guess << " ("
               << k_comb_all_failed.getVal() / k_comb_all_failed_guess << ")\n";
          cout << " - Fitted k_comb_all_passed = " << k_comb_all_passed.getVal()
               << " vs estimated " << k_comb_all_passed_guess << " ("
               << k_comb_all_passed.getVal() / k_comb_all_passed_guess << ")\n";

          cout << " - Fitted c_comb_all_failed = " << c_comb_all_failed.getVal()
               << "\n";
          cout << " - Fitted c_comb_all_passed = " << c_comb_all_passed.getVal()
               << "\n";

          cout << " - Fitted c_comb_spi_failed = " << c_comb_spi_failed.getVal()
               << "\n";
          cout << " - Fitted c_comb_spi_passed = " << c_comb_all_passed.getVal()
               << "\n";

          cout << " - Fitted non-dif passed yield = "
               << n_inmw_passed_nondif.getVal() +
                      nondif_yield_add_passed.getVal()
               << " vs estimated " << n_inmw_passed_nondif.getVal() << " ("
               << 1. + nondif_yield_add_passed.getVal() /
                           n_inmw_passed_nondif.getVal()
               << ")\n";
          cout << " - Fitted non-dif failed yield = "
               << n_inmw_failed_nondif.getVal() -
                      nondif_yield_add_passed.getVal()
               << " vs estimated " << n_inmw_failed_nondif.getVal() << " ("
               << 1. - nondif_yield_add_passed.getVal() /
                           n_inmw_failed_nondif.getVal()
               << ")\n";

          cout << " - Fitted f_dif_calib_passed = "
               << f_dif_calib_passed.getVal() << " vs estimated "
               << n_inmw_passed_dif.getVal() / (n_inmw_passed_dif.getVal() +
                                                n_inmw_passed_nondif.getVal())
               << "\n";
          cout << " - Fitted f_dif_calib_failed = "
               << f_dif_calib_failed.getVal() << " vs estimated "
               << n_inmw_failed_dif.getVal() / (n_inmw_failed_dif.getVal() +
                                                n_inmw_failed_nondif.getVal());
          cout << endl;

          d0_m_var.setRange(extended_d0_m_min, extended_d0_m_max);
          dm_var.setRange(extended_dm_min, extended_dm_max);
        }
      }
    }

    // Define output file por PID efficiency
    // File name matches the output of PIDCalib
    const TString fname_suffix =
        fake_mu ? "_denom" : (vmu ? "_nom_vmu" : "_nom");
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
    fit_status_calib.Write();
    fit_cov_qual_calib.Write();

    // Save fit status
    suffix.Form("%s_%s", year.c_str(), probe.c_str());
    c_single.cd();
    fit_status_calib.Draw("BOX2Z");
    c_single.SaveAs(opath + "/figs/fit_status_calib_" + suffix + ".pdf");
    fit_cov_qual_calib.Draw("BOX2Z");
    c_single.SaveAs(opath + "/figs/fit_cov_qual_calib_" + suffix + ".pdf");

    // Save number of fit retries
    calib_retries.Draw();
    c_single.SaveAs(opath + "/figs/fit_retries_calib_" + suffix + ".pdf");

    // Plot distribution of fitted variables
    if (!dry_run) {
      plot_dataset(ds_params_calib, opath + "/figs/params/" + suffix + "_");
    } else {
      cout << "INFO Dry run: will not plot distributions of fit parameters"
           << endl;
    }

    // Save efficiencies
    cout << "INFO Saving mass-window efficiencies " << endl;
    effs_mw_passed.Write();
    effs_mw_failed.Write();

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

  cout << "INFO Dumping counts of PASSED D0 bkg MC events" << endl;
  for (const auto &counts_probe : d0_bkg_counts_passed) {
    const auto &probe  = counts_probe.first;
    const auto &counts = counts_probe.second;
    cout << "INFO probe " << probe << "\n" << endl;
    for (unsigned p_bin = 0; p_bin < counts.size(); p_bin++) {
      cout << "INFO - p bin " << p_bin << endl;
      double bin_total = 0;
      for (const auto &counts_decay : counts[p_bin]) {
        const auto &decay = counts_decay.first;
        bin_total += counts_decay.second * w_d0_decays.at(decay);
      }
      for (const auto &counts_decay : counts[p_bin]) {
        const auto  &decay = counts_decay.first;
        const double count = counts_decay.second * w_d0_decays.at(decay);
        cout << "INFO -- " << decay << ": " << count / bin_total * 100. << "%"
             << endl;
      }
      cout << endl;
    }
  }

  cout << "INFO Dumping counts of FAILED D0 bkg MC events" << endl;
  for (const auto &counts_probe : d0_bkg_counts_failed) {
    const auto &probe  = counts_probe.first;
    const auto &counts = counts_probe.second;
    cout << "INFO probe " << probe << "\n" << endl;
    for (unsigned p_bin = 0; p_bin < counts.size(); p_bin++) {
      cout << "INFO - p bin " << p_bin << endl;
      double bin_total = 0;
      for (const auto &counts_decay : counts[p_bin]) {
        const auto &decay = counts_decay.first;
        bin_total += counts_decay.second * w_d0_decays.at(decay);
      }
      for (const auto &counts_decay : counts[p_bin]) {
        const auto  &decay = counts_decay.first;
        const double count = counts_decay.second * w_d0_decays.at(decay);
        cout << "INFO -- " << decay << ": " << count / bin_total * 100. << "%"
             << endl;
      }
      cout << endl;
    }
  }

  auto stop = high_resolution_clock::now();

  const int duration = duration_cast<seconds>(stop - start).count();
  cout << "INFO Finished processing in " << format_time(duration) << endl;
}
