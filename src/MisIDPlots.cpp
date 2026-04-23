// Author: Lucas Meyer Garcia
// License: BSD 2-clause
//
// Description: Make some useful plots to study DiF events in MC.

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include "TCanvas.h"
#include "TChain.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLine.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include "TText.h"

#include "RooArgSet.h"
#include "RooDataSet.h"

#include "misid.h"
#include "utils.h"

#include <yaml-cpp/yaml.h>
#include <cxxopts.hpp>

using ROOT::Math::PxPyPzMVector;
using ROOT::Math::XYZVector;
using std::tuple;

void plot_sd_lines(TCanvas                                       &c,
                   const map<string, tuple<double, double, int>> &sds,
                   vector<TLine *> &sd_lines, const double &y_min,
                   const double &y_max) {
  c.cd();
  for (const auto &sd : sds) {
    TLine *sd_line_before =
        new TLine(std::get<0>(sd.second), y_min, std::get<0>(sd.second), y_max);
    TLine *sd_line_after =
        new TLine(std::get<1>(sd.second), y_min, std::get<1>(sd.second), y_max);
    sd_line_before->SetLineColor(std::get<2>(sd.second));
    sd_line_after->SetLineColor(std::get<2>(sd.second));
    sd_line_before->SetLineWidth(2);
    sd_line_after->SetLineWidth(2);
    sd_line_before->SetLineStyle(kDashed);
    sd_line_after->SetLineStyle(kDashed);
    sd_lines.push_back(sd_line_before);
    sd_lines.push_back(sd_line_after);
    sd_line_before->Draw();
    sd_line_after->Draw();
  }
}

void plot_d0m_lines(TCanvas &c, vector<TLine *> &sd_lines, const TH2D *h2) {
  c.cd();

  const double x_min = h2->GetXaxis()->GetBinLowEdge(1);
  const double x_max = h2->GetXaxis()->GetBinUpEdge(h2->GetNbinsX());

  TLine *d0m_min_line = new TLine(x_min, D0_M_min, x_max, D0_M_min);
  TLine *d0m_pdg_line = new TLine(x_min, D0_M, x_max, D0_M);
  TLine *d0m_max_line = new TLine(x_min, D0_M_max, x_max, D0_M_max);

  d0m_min_line->SetLineColor(kViolet + 1);
  d0m_pdg_line->SetLineColor(kOrange + 10);
  d0m_max_line->SetLineColor(kViolet + 1);

  d0m_min_line->SetLineWidth(1);
  d0m_pdg_line->SetLineWidth(1);
  d0m_max_line->SetLineWidth(1);

  d0m_min_line->SetLineStyle(kSolid);
  d0m_pdg_line->SetLineStyle(kSolid);
  d0m_max_line->SetLineStyle(kSolid);

  sd_lines.push_back(d0m_min_line);
  sd_lines.push_back(d0m_pdg_line);
  sd_lines.push_back(d0m_max_line);

  d0m_min_line->Draw();
  d0m_pdg_line->Draw();
  d0m_max_line->Draw();
}

void plot_sd_lines(TCanvas                                       &c,
                   const map<string, tuple<double, double, int>> &sds,
                   vector<TLine *> &sd_lines, TH1D *h) {
  h->SetStats(false);
  const double y_max = h->GetMaximum();
  const double y_min = h->GetMinimum();
  plot_sd_lines(c, sds, sd_lines, y_min, y_max);
}

void plot_sd_lines(TCanvas                                       &c,
                   const map<string, tuple<double, double, int>> &sds,
                   vector<TLine *> &sd_lines, TH2D *h2) {
  h2->SetStats(false);
  const double y_max = h2->GetYaxis()->GetBinUpEdge(h2->GetNbinsY());
  const double y_min = h2->GetYaxis()->GetBinLowEdge(1);
  plot_sd_lines(c, sds, sd_lines, y_min, y_max);
}

void plot_sd_labels(TCanvas                                       &c,
                    const map<string, tuple<double, double, int>> &sds,
                    vector<TText *> &sd_labels, const double &y) {
  c.cd();
  for (const auto &sd : sds) {
    const double x = (std::get<0>(sd.second) + std::get<1>(sd.second)) / 2.;
    TText       *sd_label = new TText(x, y, sd.first.c_str());
    sd_label->SetTextColor(std::get<2>(sd.second));
    sd_label->SetTextFont(62);
    sd_label->SetTextSize(0.03);
    sd_label->SetTextAlign(kHAlignCenter + kVAlignBottom);
    sd_labels.push_back(sd_label);
    sd_label->Draw();
  }
}

void plot_sd_labels(TCanvas                                       &c,
                    const map<string, tuple<double, double, int>> &sds,
                    vector<TText *> &sd_labels, TH1D *h) {
  const double y = h->GetMaximum() * 1.002;
  plot_sd_labels(c, sds, sd_labels, y);
}

void plot_sd_labels(TCanvas                                       &c,
                    const map<string, tuple<double, double, int>> &sds,
                    vector<TText *> &sd_labels, TH2D *h2) {
  const double y = h2->GetYaxis()->GetBinUpEdge(h2->GetNbinsY()) * 1.002;
  plot_sd_labels(c, sds, sd_labels, y);
}
auto                                    a = 9500 + 2332;
map<string, tuple<double, double, int>> sds{
    {"", {RICH1_Z, RICH1_Z, kTeal}},
    {"R1", {990., 2165., kViolet}},
    {"R2", {9500., 11832., kViolet}},
    {"P", {12320., 12320., kGreen + 2}},
    {"EC", {12500., 13335. - 25, kGreen}},
    {"HC", {13335. + 25, 14990. - 25, kGray + 1}},
    {"M2", {15000. + 25, 15400., kRed}},
    {"M3", {16200., 16600., kRed}},
    {"M4", {17400., 17800., kRed}},
    {"M5", {18600., 19000., kRed}},
    {"T1", {7673., 8038., kOrange}},
    {"T2", {8360., 8725., kOrange}},
    {"T3", {9050., 9415., kOrange}},
    {"TT", {2320., 2620., kOrange}},
    {"VE", {-175., 750., kOrange}}};

int main(int argc, char **argv) {
  cxxopts::Options argOpts("GetMisIDCorrections",
                           "Calculate K/pi misid efficiencies.");

  // clang-format off
  argOpts.add_options()
    ("h,help", "Print help")
    ("p,particles", "Specify probed particle",
     cxxopts::value<vector<string>>()->default_value("pi,k"))
    ("c,config", "Specify input YAML config file",
     cxxopts::value<string>())
    ("y,year", "Specify data-taking year",
     cxxopts::value<string>()->default_value("2016"))
    ("o,output", "Specify output folder",
     cxxopts::value<string>()->default_value("gen/"))
    ("vmu", "Flag misid validation PID cuts",
     cxxopts::value<bool>()->default_value("false"))
    ("fake_mu", "Flag fake muon sample PID cuts",
     cxxopts::value<bool>()->default_value("false"))
    ("r,run2ang", "Use Run2Ang PID cuts",
     cxxopts::value<bool>()->default_value("false"))
    ;
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  // Read arguments
  const auto ymlFile   = parsedArgs["config"].as<string>();
  const auto particles = parsedArgs["particles"].as<vector<string>>();
  const auto vmu       = parsedArgs["vmu"].as<bool>();
  const auto fake_mu   = parsedArgs["fake_mu"].as<bool>();
  const auto year      = parsedArgs["year"].as<string>();
  const auto run2ang   = parsedArgs["run2ang"].as<bool>();

  if (run2ang) {
    cout << "INFO Using Run 2 Angular Analysis cuts" << endl;
  } else {
    cout << "INFO Using Run 1 RDx Analysis cuts" << endl;
  }

  // Parse YAML config
  const auto ymlConfig = YAML::LoadFile(ymlFile)["misid_corrections"];

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

  const Sample sample      = fake_mu ? FAKE_MU : (vmu ? VMU : ISO_CTRL);
  const string sample_name = SAMPLES.at(sample);

  // Observables
  RooRealVar d0_m_var("d0_m_var", "M_{D^{0}}", D0_M_min, D0_M_max, "GeV");
  RooRealVar dm_var("dm_var", "dm", DM_min, DM_max, "GeV");

  RooArgSet fit_vars(d0_m_var, dm_var);

  // Define histogram to easily determine kinematical bins
  TH3D histo_binning("histo_binning", ";p;#eta;nTracks", N_BINS_P, BINS_P,
                     N_BINS_ETA, BINS_ETA, N_BINS_NTRACKS, BINS_NTRACKS);

  TCanvas c("c", "c", 1280, 960);

  for (auto probe : particles) {
    cout << "INFO Selecting " << probe << endl;

    const TString tag = (probe == "pi") ? "k" : "pi";

    RooDataSet *datasets_mc_passed_dif[N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] = {
        {{nullptr}}};
    RooDataSet *datasets_mc_failed_dif[N_BINS_NTRACKS][N_BINS_ETA][N_BINS_P] = {
        {{nullptr}}};
    RooDataSet *datasets_mc_passed_nondif[N_BINS_NTRACKS][N_BINS_ETA]
                                         [N_BINS_P] = {{{nullptr}}};
    RooDataSet *datasets_mc_failed_nondif[N_BINS_NTRACKS][N_BINS_ETA]
                                         [N_BINS_P] = {{{nullptr}}};

    // Initialize MC datasets
    // Note: In the fits, we merge ntrack bins for dif events, but here we keep
    // them separate for plotting
    cout << "INFO Initializing signal MC datahists" << endl;
    for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
      for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
        for (int ntrks_idx = 0; ntrks_idx < N_BINS_NTRACKS; ntrks_idx++) {
          TString suffix = TString::Format("%s_%d_%d_%d", probe.c_str(),
                                           ntrks_idx, eta_idx, p_idx);
          datasets_mc_passed_dif[ntrks_idx][eta_idx][p_idx] =
              new RooDataSet("ds_mc_passed_dif_" + suffix,
                             "ds_mc_passed_dif_" + suffix, fit_vars);
          datasets_mc_failed_dif[ntrks_idx][eta_idx][p_idx] =
              new RooDataSet("ds_mc_failed_dif_" + suffix,
                             "ds_mc_failed_dif_" + suffix, fit_vars);
          datasets_mc_passed_nondif[ntrks_idx][eta_idx][p_idx] =
              new RooDataSet("ds_mc_passed_nondif_" + suffix,
                             "ds_mc_passed_nondif_" + suffix, fit_vars);
          datasets_mc_failed_nondif[ntrks_idx][eta_idx][p_idx] =
              new RooDataSet("ds_mc_failed_nondif_" + suffix,
                             "ds_mc_failed_nondif_" + suffix, fit_vars);
        }
      }
    }

    // Create additional histograms
    TH1D *h_dif_z_all_nocuts[N_BINS_ETA][N_BINS_P];
    TH1D *h_dif_z_all[N_BINS_ETA][N_BINS_P];
    TH1D *h_dif_z_ismuon[N_BINS_ETA][N_BINS_P];
    TH1D *h_dif_z_ismuon_dllmu[N_BINS_ETA][N_BINS_P];
    TH1D *h_dif_z_ismuon_mubdt[N_BINS_ETA][N_BINS_P];
    TH1D *h_dif_z_ismuon_dllmu_mubdt[N_BINS_ETA][N_BINS_P];

    TH1D *h_mu_p_all_nocuts[N_BINS_ETA][N_BINS_P];
    TH1D *h_mu_p_all[N_BINS_ETA][N_BINS_P];
    TH1D *h_mu_p_ismuon[N_BINS_ETA][N_BINS_P];
    TH1D *h_mu_p_ismuon_dllmu[N_BINS_ETA][N_BINS_P];
    TH1D *h_mu_p_ismuon_mubdt[N_BINS_ETA][N_BINS_P];
    TH1D *h_mu_p_ismuon_dllmu_mubdt[N_BINS_ETA][N_BINS_P];

    TH1D *h_dif_d0m_hadron_norich[N_BINS_ETA][N_BINS_P];
    TH1D *h_dif_d0m_hadron_rich1[N_BINS_ETA][N_BINS_P];
    TH1D *h_dif_d0m_hadron_rich12[N_BINS_ETA][N_BINS_P];
    TH1D *h_nondif_d0m[N_BINS_ETA][N_BINS_P];

    TH2D *h2_dif_d0m_z_all_nocuts[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_d0m_z_all[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_d0m_z_ismuon[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_d0m_z_ismuon_dllmu[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_d0m_z_ismuon_mubdt[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_d0m_z_ismuon_dllmu_mubdt[N_BINS_ETA][N_BINS_P];

    TH2D *h2_dif_alpha_z_all_nocuts[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_alpha_z_all_nomw[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_alpha_z_all[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_alpha_z_ismuon[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_alpha_z_ismuon_dllmu[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_alpha_z_ismuon_mubdt[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_alpha_z_ismuon_dllmu_mubdt[N_BINS_ETA][N_BINS_P];

    TH2D *h2_dif_dira_z_all_nocuts[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_dira_z_all_nomw[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_dira_z_all[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_dira_z_ismuon[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_dira_z_ismuon_dllmu[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_dira_z_ismuon_mubdt[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_dira_z_ismuon_dllmu_mubdt[N_BINS_ETA][N_BINS_P];

    TH2D *h2_dif_track_chi2ndof_z_all_nocuts[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_track_chi2ndof_z_all_nomw[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_track_chi2ndof_z_all[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_track_chi2ndof_z_ismuon[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_track_chi2ndof_z_ismuon_dllmu[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_track_chi2ndof_z_ismuon_mubdt[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_track_chi2ndof_z_ismuon_dllmu_mubdt[N_BINS_ETA][N_BINS_P];

    TH2D *h2_dif_had_p_z_all_nocuts[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_had_p_z_all_nomw[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_had_p_z_all[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_had_p_z_ismuon[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_had_p_z_ismuon_dllmu[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_had_p_z_ismuon_mubdt[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_had_p_z_ismuon_dllmu_mubdt[N_BINS_ETA][N_BINS_P];

    TH2D *h2_dif_mu_p_z_all_nocuts[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_mu_p_z_all_nomw[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_mu_p_z_all[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_mu_p_z_ismuon[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_mu_p_z_ismuon_dllmu[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_mu_p_z_ismuon_mubdt[N_BINS_ETA][N_BINS_P];
    TH2D *h2_dif_mu_p_z_ismuon_dllmu_mubdt[N_BINS_ETA][N_BINS_P];

    for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
      for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
        const TString suffix =
            TString::Format("%s_%d_%d", probe.c_str(), eta_idx, p_idx);
        const TString tag =
            TString::Format("(%.1f < #eta < %.1f, %.0f < p < %.0f)",

                            BINS_ETA[eta_idx], BINS_ETA[eta_idx + 1],
                            BINS_P[p_idx], BINS_P[p_idx + 1]);

        h_dif_z_all_nocuts[eta_idx][p_idx] =
            new TH1D("h_dif_z_all_nocuts_" + suffix, tag + ";z_{DiF} [mm]", 160,
                     -500., 19500.);
        h_dif_z_all[eta_idx][p_idx] = new TH1D(
            "h_dif_z_all_" + suffix, tag + ";z_{DiF} [mm]", 160, -500., 19500.);
        h_dif_z_ismuon[eta_idx][p_idx] =
            new TH1D("h_dif_z_ismuon_" + suffix, tag + ";z_{DiF} [mm]", 160,
                     -500., 19500.);
        h_dif_z_ismuon_dllmu[eta_idx][p_idx] =
            new TH1D("h_dif_z_ismuon_dllmu_" + suffix, tag + ";z_{DiF} [mm]",
                     160, -500., 19500.);
        h_dif_z_ismuon_mubdt[eta_idx][p_idx] =
            new TH1D("h_dif_z_ismuon_mubdt_" + suffix, tag + ";z_{DiF} [mm]",
                     160, -500., 19500.);
        h_dif_z_ismuon_dllmu_mubdt[eta_idx][p_idx] =
            new TH1D("h_dif_z_ismuon_dllmu_mubdt_" + suffix,
                     tag + ";z_{DiF} [mm]", 160, -500., 19500.);

        h_mu_p_all_nocuts[eta_idx][p_idx] =
            new TH1D("h_mu_p_all_nocuts_" + suffix, tag + ";p_{#mu} [MeV]", 60,
                     1000., 7000.);
        h_mu_p_all[eta_idx][p_idx] = new TH1D(
            "h_mu_p_all_" + suffix, tag + ";p_{#mu} [MeV]", 60, 1000., 7000.);
        h_mu_p_ismuon[eta_idx][p_idx] =
            new TH1D("h_mu_p_ismuon_" + suffix, tag + ";p_{#mu} [MeV]", 60,
                     1000., 7000.);
        h_mu_p_ismuon_dllmu[eta_idx][p_idx] =
            new TH1D("h_mu_p_ismuon_dllmu_" + suffix, tag + ";p_{#mu} [MeV]",
                     60, 1000., 7000.);
        h_mu_p_ismuon_mubdt[eta_idx][p_idx] =
            new TH1D("h_mu_p_ismuon_mubdt_" + suffix, tag + ";p_{#mu} [MeV]",
                     60, 1000., 7000.);
        h_mu_p_ismuon_dllmu_mubdt[eta_idx][p_idx] =
            new TH1D("h_mu_p_ismuon_dllmu_mubdt_" + suffix,
                     tag + ";p_{#mu} [MeV]", 60, 1000., 7000.);

        h_dif_d0m_hadron_norich[eta_idx][p_idx] =
            new TH1D("h_dif_d0m_hadron_norich_" + suffix,
                     tag + ";m(D^{0}) [GeV]", 50, 1.825, 1.910);
        h_dif_d0m_hadron_rich1[eta_idx][p_idx] =
            new TH1D("h_dif_d0m_hadron_rich1_" + suffix,
                     tag + ";m(D^{0}) [GeV]", 50, 1.825, 1.910);
        h_dif_d0m_hadron_rich12[eta_idx][p_idx] =
            new TH1D("h_dif_d0m_hadron_rich12_" + suffix,
                     tag + ";m(D^{0}) [GeV]", 50, 1.825, 1.910);
        h_nondif_d0m[eta_idx][p_idx] =
            new TH1D("h_nondif_d0m_" + suffix, tag + ";m(D^{0}) [GeV]", 50,
                     1.825, 1.910);

        h2_dif_d0m_z_all_nocuts[eta_idx][p_idx] =
            new TH2D("h2_dif_d0m_z_all_nocuts_" + suffix,
                     tag + ";z_{DiF} [mm];m(D^{0}) [GeV]", 100, -600., 19400.,
                     50, 1., 2.);
        h2_dif_d0m_z_all[eta_idx][p_idx] = new TH2D(
            "h2_dif_d0m_z_all_" + suffix, tag + ";z_{DiF} [mm];m(D^{0}) [GeV]",
            100, -600., 19400., 50, 1., 2.);
        h2_dif_d0m_z_ismuon[eta_idx][p_idx] =
            new TH2D("h2_dif_d0m_z_ismuon_" + suffix,
                     tag + ";z_{DiF} [mm];m(D^{0}) [GeV]", 100, -600., 19400.,
                     50, 1., 2.);
        h2_dif_d0m_z_ismuon_dllmu[eta_idx][p_idx] =
            new TH2D("h2_dif_d0m_z_ismuon_dllmu_" + suffix,
                     tag + ";z_{DiF} [mm];m(D^{0}) [GeV]", 100, -600., 19400.,
                     50, 1., 2.);
        h2_dif_d0m_z_ismuon_mubdt[eta_idx][p_idx] =
            new TH2D("h2_dif_d0m_z_ismuon_mubdt_" + suffix,
                     tag + ";z_{DiF} [mm];m(D^{0}) [GeV]", 100, -600., 19400.,
                     50, 1., 2.);
        h2_dif_d0m_z_ismuon_dllmu_mubdt[eta_idx][p_idx] =
            new TH2D("h2_dif_d0m_z_ismuon_dllmu_mubdt_" + suffix,
                     tag + ";z_{DiF} [mm];m(D^{0}) [GeV]", 100, -600., 19400.,
                     50, 1., 2.);

        h2_dif_alpha_z_all_nocuts[eta_idx][p_idx] =
            new TH2D("h2_dif_alpha_z_all_nocuts_" + suffix,
                     tag + ";z_{DiF} [mm];#alpha_{h,#mu}(#circ)", 100, -600.,
                     19400., 50, 0., 1.5);
        h2_dif_alpha_z_all_nomw[eta_idx][p_idx] =
            new TH2D("h2_dif_alpha_z_all_nomw_" + suffix,
                     tag + ";z_{DiF} [mm];#alpha_{h,#mu}(#circ)", 100, -600.,
                     19400., 50, 0., 1.5);
        h2_dif_alpha_z_all[eta_idx][p_idx] =
            new TH2D("h2_dif_alpha_z_all_" + suffix,
                     tag + ";z_{DiF} [mm];#alpha_{h,#mu}(#circ)", 100, -600.,
                     19400., 50, 0., 1.5);
        h2_dif_alpha_z_ismuon[eta_idx][p_idx] =
            new TH2D("h2_dif_alpha_z_ismuon_" + suffix,
                     tag + ";z_{DiF} [mm];#alpha_{h,#mu}(#circ)", 100, -600.,
                     19400., 50, 0., 1.5);
        h2_dif_alpha_z_ismuon_dllmu[eta_idx][p_idx] =
            new TH2D("h2_dif_alpha_z_ismuon_dllmu_" + suffix,
                     tag + ";z_{DiF} [mm];#alpha_{h,#mu}(#circ)", 100, -600.,
                     19400., 50, 0., 1.5);
        h2_dif_alpha_z_ismuon_mubdt[eta_idx][p_idx] =
            new TH2D("h2_dif_alpha_z_ismuon_mubdt_" + suffix,
                     tag + ";z_{DiF} [mm];#alpha_{h,#mu}(#circ)", 100, -600.,
                     19400., 50, 0., 1.5);
        h2_dif_alpha_z_ismuon_dllmu_mubdt[eta_idx][p_idx] =
            new TH2D("h2_dif_alpha_z_ismuon_dllmu_mubdt_" + suffix,
                     tag + ";z_{DiF} [mm];#alpha_{h,#mu}(#circ)", 100, -600.,
                     19400., 50, 0., 1.5);

        h2_dif_dira_z_all_nocuts[eta_idx][p_idx] = new TH2D(
            "h2_dif_dira_z_all_nocuts_" + suffix,
            tag + ";z_{DiF} [mm];DIRA_{PV}", 100, -600., 19400., 50, 0.999, 1.);
        h2_dif_dira_z_all_nomw[eta_idx][p_idx] = new TH2D(
            "h2_dif_dira_z_all_nomw_" + suffix, tag + ";z_{DiF} [mm];DIRA_{PV}",
            100, -600., 19400., 50, 0.999, 1.);
        h2_dif_dira_z_all[eta_idx][p_idx] = new TH2D(
            "h2_dif_dira_z_all_" + suffix, tag + ";z_{DiF} [mm];DIRA_{PV}", 100,
            -600., 19400., 50, 0.999, 1.);
        h2_dif_dira_z_ismuon[eta_idx][p_idx] = new TH2D(
            "h2_dif_dira_z_ismuon_" + suffix, tag + ";z_{DiF} [mm];DIRA_{PV}",
            100, -600., 19400., 50, 0.999, 1.);
        h2_dif_dira_z_ismuon_dllmu[eta_idx][p_idx] = new TH2D(
            "h2_dif_dira_z_ismuon_dllmu_" + suffix,
            tag + ";z_{DiF} [mm];DIRA_{PV}", 100, -600., 19400., 50, 0.999, 1.);
        h2_dif_dira_z_ismuon_mubdt[eta_idx][p_idx] = new TH2D(
            "h2_dif_dira_z_ismuon_mubdt_" + suffix,
            tag + ";z_{DiF} [mm];DIRA_{PV}", 100, -600., 19400., 50, 0.999, 1.);
        h2_dif_dira_z_ismuon_dllmu_mubdt[eta_idx][p_idx] = new TH2D(
            "h2_dif_dira_z_ismuon_dllmu_mubdt_" + suffix,
            tag + ";z_{DiF} [mm];DIRA_{PV}", 100, -600., 19400., 50, 0.999, 1.);

        h2_dif_track_chi2ndof_z_all_nocuts[eta_idx][p_idx] =
            new TH2D("h2_dif_track_chi2ndof_z_all_nocuts_" + suffix,
                     tag + ";z_{DiF} [mm];Track #chi^{2}_{NDOF}", 100, -600.,
                     19400., 50, 0., 5.);
        h2_dif_track_chi2ndof_z_all_nomw[eta_idx][p_idx] =
            new TH2D("h2_dif_track_chi2ndof_z_all_nomw_" + suffix,
                     tag + ";z_{DiF} [mm];Track #chi^{2}_{NDOF}", 100, -600.,
                     19400., 50, 0., 5.);
        h2_dif_track_chi2ndof_z_all[eta_idx][p_idx] =
            new TH2D("h2_dif_track_chi2ndof_z_all_" + suffix,
                     tag + ";z_{DiF} [mm];Track #chi^{2}_{NDOF}", 100, -600.,
                     19400., 50, 0., 5.);
        h2_dif_track_chi2ndof_z_ismuon[eta_idx][p_idx] =
            new TH2D("h2_dif_track_chi2ndof_z_ismuon_" + suffix,
                     tag + ";z_{DiF} [mm];Track #chi^{2}_{NDOF}", 100, -600.,
                     19400., 50, 0., 5.);
        h2_dif_track_chi2ndof_z_ismuon_dllmu[eta_idx][p_idx] =
            new TH2D("h2_dif_track_chi2ndof_z_ismuon_dllmu_" + suffix,
                     tag + ";z_{DiF} [mm];Track #chi^{2}_{NDOF}", 100, -600.,
                     19400., 50, 0., 5.);
        h2_dif_track_chi2ndof_z_ismuon_mubdt[eta_idx][p_idx] =
            new TH2D("h2_dif_track_chi2ndof_z_ismuon_mubdt_" + suffix,
                     tag + ";z_{DiF} [mm];Track #chi^{2}_{NDOF}", 100, -600.,
                     19400., 50, 0., 5.);
        h2_dif_track_chi2ndof_z_ismuon_dllmu_mubdt[eta_idx][p_idx] =
            new TH2D("h2_dif_track_chi2ndof_z_ismuon_dllmu_mubdt_" + suffix,
                     tag + ";z_{DiF} [mm];Track #chi^{2}_{NDOF}", 100, -600.,
                     19400., 50, 0., 5.);

        h2_dif_had_p_z_all_nocuts[eta_idx][p_idx] =
            new TH2D("h2_dif_had_p_z_all_nocuts_" + suffix,
                     tag + ";z_{DiF} [mm];p_{had} [MeV]", 100, -600., 19400.,
                     60, 1000., 7000.);
        h2_dif_had_p_z_all_nomw[eta_idx][p_idx] =
            new TH2D("h2_dif_had_p_z_all_nomw_" + suffix,
                     tag + ";z_{DiF} [mm];p_{had} [MeV]", 100, -600., 19400.,
                     60, 1000., 7000.);
        h2_dif_had_p_z_all[eta_idx][p_idx] = new TH2D(
            "h2_dif_had_p_z_all_" + suffix, tag + ";z_{DiF} [mm];p_{had} [MeV]",
            100, -600., 19400., 60, 1000., 7000.);
        h2_dif_had_p_z_ismuon[eta_idx][p_idx] =
            new TH2D("h2_dif_had_p_z_ismuon_" + suffix,
                     tag + ";z_{DiF} [mm];p_{had} [MeV]", 100, -600., 19400.,
                     60, 1000., 7000.);
        h2_dif_had_p_z_ismuon_dllmu[eta_idx][p_idx] =
            new TH2D("h2_dif_had_p_z_ismuon_dllmu_" + suffix,
                     tag + ";z_{DiF} [mm];p_{had} [MeV]", 100, -600., 19400.,
                     60, 1000., 7000.);
        h2_dif_had_p_z_ismuon_mubdt[eta_idx][p_idx] =
            new TH2D("h2_dif_had_p_z_ismuon_mubdt_" + suffix,
                     tag + ";z_{DiF} [mm];p_{had} [MeV]", 100, -600., 19400.,
                     60, 1000., 7000.);
        h2_dif_had_p_z_ismuon_dllmu_mubdt[eta_idx][p_idx] =
            new TH2D("h2_dif_had_p_z_ismuon_dllmu_mubdt_" + suffix,
                     tag + ";z_{DiF} [mm];p_{had} [MeV]", 100, -600., 19400.,
                     60, 1000., 7000.);

        h2_dif_mu_p_z_all_nocuts[eta_idx][p_idx] =
            new TH2D("h2_dif_mu_p_z_all_nocuts_" + suffix,
                     tag + ";z_{DiF} [mm];p_{#mu} [MeV]", 100, -600., 19400.,
                     60, 1000., 7000.);
        h2_dif_mu_p_z_all_nomw[eta_idx][p_idx] =
            new TH2D("h2_dif_mu_p_z_all_nomw_" + suffix,
                     tag + ";z_{DiF} [mm];p_{#mu} [MeV]", 100, -600., 19400.,
                     60, 1000., 7000.);
        h2_dif_mu_p_z_all[eta_idx][p_idx] = new TH2D(
            "h2_dif_mu_p_z_all_" + suffix, tag + ";z_{DiF} [mm];p_{#mu} [MeV]",
            100, -600., 19400., 60, 1000., 7000.);
        h2_dif_mu_p_z_ismuon[eta_idx][p_idx] =
            new TH2D("h2_dif_mu_p_z_ismuon_" + suffix,
                     tag + ";z_{DiF} [mm];p_{#mu} [MeV]", 100, -600., 19400.,
                     60, 1000., 7000.);
        h2_dif_mu_p_z_ismuon_dllmu[eta_idx][p_idx] =
            new TH2D("h2_dif_mu_p_z_ismuon_dllmu_" + suffix,
                     tag + ";z_{DiF} [mm];p_{#mu} [MeV]", 100, -600., 19400.,
                     60, 1000., 7000.);
        h2_dif_mu_p_z_ismuon_mubdt[eta_idx][p_idx] =
            new TH2D("h2_dif_mu_p_z_ismuon_mubdt_" + suffix,
                     tag + ";z_{DiF} [mm];p_{#mu} [MeV]", 100, -600., 19400.,
                     60, 1000., 7000.);
        h2_dif_mu_p_z_ismuon_dllmu_mubdt[eta_idx][p_idx] =
            new TH2D("h2_dif_mu_p_z_ismuon_dllmu_mubdt_" + suffix,
                     tag + ";z_{DiF} [mm];p_{#mu} [MeV]", 100, -600., 19400.,
                     60, 1000., 7000.);
      }
    }

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
      int dst_id, d0_id, probe_trueid, probe_daughter0_trueid,
          probe_daughter1_trueid, ntracks, dst_vtx_ndof, d0_vtx_ndof,
          probe_mc_mom_nd, probe_mc_mom_id, probe_mc_gd_mom_id, tag_trueid,
          spi_trueid;

      double dst_m, d0_m, probe_p, probe_pz, probe_pt, probe_dllmu, probe_dlle,
          spi_p, spi_pt, tag_p, tag_pz, tag_pt, k_track_chi2ndof,
          pi_track_chi2ndof, spi_track_chi2ndof, k_px, k_py, pi_px, pi_py,
          dst_vtx_chi2, d0_vtx_chi2, k_ghostprob, pi_ghostprob, spi_ghostprob,
          probe_daughter0_true_origz, probe_true_origz, probe_mother_true_px,
          probe_mother_true_py, probe_mother_true_pz, probe_true_px,
          probe_true_py, probe_true_pz, probe_daughter0_true_px,
          probe_daughter0_true_py, probe_daughter0_true_pz,
          probe_daughter1_true_px, probe_daughter1_true_py,
          probe_daughter1_true_pz, d0_dira;

      float probe_mu_ubdt;

      bool probe_ismuon, probe_hasmuon;

      cout << "INFO Setting input branches " << endl;
      ch_mc.SetBranchStatus("*", false);
      ch_mc.SetBranchAddress("dst_M", &dst_m);
      ch_mc.SetBranchAddress("dst_TRUEID", &dst_id);
      ch_mc.SetBranchAddress("dst_ENDVERTEX_CHI2", &dst_vtx_chi2);
      ch_mc.SetBranchAddress("dst_ENDVERTEX_NDOF", &dst_vtx_ndof);
      ch_mc.SetBranchAddress("d0_TRUEID", &d0_id);
      ch_mc.SetBranchAddress("d0_M", &d0_m);
      ch_mc.SetBranchAddress("d0_ENDVERTEX_CHI2", &d0_vtx_chi2);
      ch_mc.SetBranchAddress("d0_ENDVERTEX_NDOF", &d0_vtx_ndof);
      ch_mc.SetBranchAddress("spi_P", &spi_p);
      ch_mc.SetBranchAddress("spi_PT", &spi_pt);
      ch_mc.SetBranchAddress("spi_TRUEID", &spi_trueid);
      ch_mc.SetBranchAddress(tag + "_P", &tag_p);
      ch_mc.SetBranchAddress(tag + "_PT", &tag_pt);
      ch_mc.SetBranchAddress(tag + "_PZ", &tag_pz);
      ch_mc.SetBranchAddress(tag + "_TRUEID", &tag_trueid);
      ch_mc.SetBranchAddress((probe + "_P").c_str(), &probe_p);
      ch_mc.SetBranchAddress((probe + "_PT").c_str(), &probe_pt);
      ch_mc.SetBranchAddress((probe + "_PZ").c_str(), &probe_pz);
      ch_mc.SetBranchAddress((probe + "_isMuon").c_str(), &probe_ismuon);
      ch_mc.SetBranchAddress((probe + "_hasMuon").c_str(), &probe_hasmuon);
      ch_mc.SetBranchAddress((probe + "_PIDmu").c_str(), &probe_dllmu);
      ch_mc.SetBranchAddress((probe + "_PIDe").c_str(), &probe_dlle);
      ch_mc.SetBranchAddress((probe + "_bdt_mu").c_str(), &probe_mu_ubdt);
      ch_mc.SetBranchAddress((probe + "_TRUEID").c_str(), &probe_trueid);
      ch_mc.SetBranchAddress((probe + "_TRUEORIGINVERTEX_Z").c_str(),
                             &probe_true_origz);
      ch_mc.SetBranchAddress((probe + "_DAUGHTER0_ID").c_str(),
                             &probe_daughter0_trueid);
      ch_mc.SetBranchAddress((probe + "_DAUGHTER0_ORIGINVERTEX_Z").c_str(),
                             &probe_daughter0_true_origz);
      ch_mc.SetBranchAddress((probe + "_DAUGHTER1_ID").c_str(),
                             &probe_daughter1_trueid);
      ch_mc.SetBranchAddress("k_PX", &k_px);
      ch_mc.SetBranchAddress("k_PY", &k_py);
      ch_mc.SetBranchAddress("pi_PX", &pi_px);
      ch_mc.SetBranchAddress("pi_PY", &pi_py);
      ch_mc.SetBranchAddress("k_TRACK_CHI2NDOF", &k_track_chi2ndof);
      ch_mc.SetBranchAddress("pi_TRACK_CHI2NDOF", &pi_track_chi2ndof);
      ch_mc.SetBranchAddress("spi_TRACK_CHI2NDOF", &spi_track_chi2ndof);
      ch_mc.SetBranchAddress("k_TRACK_GhostProb", &k_ghostprob);
      ch_mc.SetBranchAddress("pi_TRACK_GhostProb", &pi_ghostprob);
      ch_mc.SetBranchAddress("spi_TRACK_GhostProb", &spi_ghostprob);
      ch_mc.SetBranchAddress("nTracks", &ntracks);
      ch_mc.SetBranchAddress((probe + "_MC_MOTHER_ID").c_str(),
                             &probe_mc_mom_id);
      ch_mc.SetBranchAddress((probe + "_MC_MOTHER_ND").c_str(),
                             &probe_mc_mom_nd);
      ch_mc.SetBranchAddress((probe + "_MC_MOTHER_TRUEPX").c_str(),
                             &probe_mother_true_px);
      ch_mc.SetBranchAddress((probe + "_MC_MOTHER_TRUEPY").c_str(),
                             &probe_mother_true_py);
      ch_mc.SetBranchAddress((probe + "_MC_MOTHER_TRUEPZ").c_str(),
                             &probe_mother_true_pz);
      ch_mc.SetBranchAddress((probe + "_MC_GD_MOTHER_ID").c_str(),
                             &probe_mc_gd_mom_id);
      ch_mc.SetBranchAddress((probe + "_TRUEP_X").c_str(), &probe_true_px);
      ch_mc.SetBranchAddress((probe + "_TRUEP_Y").c_str(), &probe_true_py);
      ch_mc.SetBranchAddress((probe + "_TRUEP_Z").c_str(), &probe_true_pz);
      ch_mc.SetBranchAddress((probe + "_DAUGHTER0_PX").c_str(),
                             &probe_daughter0_true_px);
      ch_mc.SetBranchAddress((probe + "_DAUGHTER0_PY").c_str(),
                             &probe_daughter0_true_py);
      ch_mc.SetBranchAddress((probe + "_DAUGHTER0_PZ").c_str(),
                             &probe_daughter0_true_pz);
      ch_mc.SetBranchAddress((probe + "_DAUGHTER1_PX").c_str(),
                             &probe_daughter1_true_px);
      ch_mc.SetBranchAddress((probe + "_DAUGHTER1_PY").c_str(),
                             &probe_daughter1_true_py);
      ch_mc.SetBranchAddress((probe + "_DAUGHTER1_PZ").c_str(),
                             &probe_daughter1_true_pz);

      ch_mc.SetBranchAddress("d0_DIRA_OWNPV", &d0_dira);

      const int entries_mc = ch_mc.GetEntries();
      cout << "INFO Starting MC event loop over " << entries_mc << " entries"
           << endl;

      int eta_bin = 0, p_bin = 0, ntrks_bin = 0;

      for (int evt = 0; evt < entries_mc; evt++) {
        ch_mc.GetEntry(evt);

        // Offline truth-matching
        if (std::abs(dst_id) != Dst_ID) continue;
        if (std::abs(d0_id) != D0_ID) continue;

        // Deduce hadron ID
        int probe_id = 0, tag_id = 0;
        if (probe == "k") {
          probe_id = K_ID;
          tag_id   = PI_ID;
        } else if (probe == "pi") {
          probe_id = PI_ID;
          tag_id   = K_ID;
        }

        if (std::abs(probe_trueid) != probe_id &&
            std::abs(probe_trueid) != MU_ID)
          continue;
        if (std::abs(tag_trueid) != tag_id && std::abs(tag_trueid) != MU_ID)
          continue;
        if (std::abs(spi_trueid) != PI_ID && std::abs(spi_trueid) != MU_ID)
          continue;

        // Check for decay in flight of probe hadron
        const bool dif =
            check_dif(probe_trueid, probe_daughter0_trueid,
                      probe_daughter1_trueid, probe_mc_mom_nd, probe);

        const bool probe_matched_mu = std::abs(probe_trueid) == MU_ID;

        // if (!probe_matched_mu) continue;  // FIXME Just for tests!
        // if (probe_mc_mom_nd == 2) continue;//

        double alpha     = 2.95;
        double mu_true_p = 0.;
        double h_true_p  = 0.;

        if (dif && !probe_matched_mu) {
          const XYZVector mu_vec3(probe_daughter0_true_px,
                                  probe_daughter0_true_py,
                                  probe_daughter0_true_pz);
          const XYZVector nu_vec3(probe_daughter1_true_px,
                                  probe_daughter1_true_py,
                                  probe_daughter1_true_pz);
          const XYZVector h_vec3 = mu_vec3 + nu_vec3;

          const double cos_alpha =
              h_vec3.Dot(mu_vec3) / std::sqrt(h_vec3.Mag2() * mu_vec3.Mag2());
          alpha = std::acos(cos_alpha) * 360. / TMath::TwoPi();  // degrees

          mu_true_p =
              sqrt_sum_sq(probe_daughter0_true_px, probe_daughter0_true_py,
                          probe_daughter0_true_pz);
          h_true_p = sqrt_sum_sq(probe_true_px, probe_true_py, probe_true_pz);
        }
        if (dif && probe_matched_mu) {
          mu_true_p = sqrt_sum_sq(probe_true_px, probe_true_py, probe_true_pz);
          h_true_p  = sqrt_sum_sq(probe_mother_true_px, probe_mother_true_py,
                                  probe_mother_true_pz);
        }

        const double dif_z =
            probe_matched_mu ? probe_true_origz : probe_daughter0_true_origz;

        const double dm = (dst_m - d0_m) * 0.001;
        d0_m            = d0_m * 0.001;

        if (probe_p < 3000. || probe_p >= 100000.) continue;
        const double probe_eta =
            0.5 * std::log((probe_p + probe_pz) / (probe_p - probe_pz));
        if (probe_eta < 1.7 || probe_eta >= 5.0) continue;

        // Determine kinematical bin
        const int kin_bin = histo_binning.FindBin(probe_p, probe_eta, ntracks);
        histo_binning.GetBinXYZ(kin_bin, p_bin, eta_bin, ntrks_bin);

        const double probe_track_chi2ndof =
            (probe == "k") ? k_track_chi2ndof : pi_track_chi2ndof;

        if (dif) {
          h_mu_p_all_nocuts[eta_bin - 1][p_bin - 1]->Fill(mu_true_p);
          h_dif_z_all_nocuts[eta_bin - 1][p_bin - 1]->Fill(dif_z);
          h2_dif_d0m_z_all_nocuts[eta_bin - 1][p_bin - 1]->Fill(dif_z, d0_m);
          h2_dif_alpha_z_all_nocuts[eta_bin - 1][p_bin - 1]->Fill(dif_z, alpha);
          h2_dif_dira_z_all_nocuts[eta_bin - 1][p_bin - 1]->Fill(dif_z,
                                                                 d0_dira);
          h2_dif_mu_p_z_all_nocuts[eta_bin - 1][p_bin - 1]->Fill(dif_z,
                                                                 mu_true_p);
          h2_dif_had_p_z_all_nocuts[eta_bin - 1][p_bin - 1]->Fill(dif_z,
                                                                  h_true_p);
          h2_dif_track_chi2ndof_z_all_nocuts[eta_bin - 1][p_bin - 1]->Fill(
              dif_z, probe_track_chi2ndof);
        }

        // Conditional cuts
        if (!probe_hasmuon) continue;

        if ((k_ghostprob > 0.4) || (pi_ghostprob > 0.4) ||
            (spi_ghostprob > 0.4))
          continue;

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

        if ((dst_vtx_chi2 / dst_vtx_ndof) > 15 ||
            (d0_vtx_chi2 / d0_vtx_ndof) > 10)
          continue;

        double k_pz, pi_pz;
        if (probe == "k") {
          k_pz  = probe_pz;
          pi_pz = tag_pz;
        } else {
          k_pz  = tag_pz;
          pi_pz = probe_pz;
        }

        const PxPyPzMVector p_k(k_px, k_py, k_pz, K_M * 1000.);
        const PxPyPzMVector p_pi(pi_px, pi_py, pi_pz, PI_M * 1000.);

        if ((std::max(probe_pt, tag_pt) < 1000.) || ((p_k + p_pi).Pt() < 1500.))
          continue;

        const PxPyPzMVector p_k_wm(k_px, k_py, k_pz, PI_M * 1000.);
        const PxPyPzMVector p_pi_wm(pi_px, pi_py, pi_pz, K_M * 1000.);

        // Veto wrong-mass hypotheses
        // pi pi
        const double wm_pipi = (p_k_wm + p_pi).M();
        if (std::abs(wm_pipi - D0_M * 1000.) < 25.) continue;

        // pi K (wrong sign)
        const double wm_pik = (p_k_wm + p_pi_wm).M();
        if (std::abs(wm_pik - D0_M * 1000.) < 25.) continue;

        // K K
        const double wm_kk = (p_k + p_pi_wm).M();
        if (std::abs(wm_kk - D0_M * 1000.) < 25.) continue;

        // MuonUnbiased equivalent to L0+HLT1 TIS. See MuonUnbiased defition
        // at
        // https://gitlab.cern.ch/lhcb-datapkg/WG/PIDCalib/-/blob/master/scriptsR2/makeTuples_pp_2016_reprocessing.py#L71
        // and
        // https://mattermost.web.cern.ch/lhcb/pl/893yre484jggigooti5u3gqb8w

        // Not applying MuonUnbiased and GHOSTPROB conditional cuts in MC
        // since some kinematical bins have very low statistics

        if (probe_p < 3000. || probe_p >= 100000. || ntracks >= 600) continue;
        if (probe_eta < 1.7 || probe_eta >= 5.0) continue;

        // Fill histograms

        if (dif && in_range(DM_min, dm, DM_max)) {
          h2_dif_d0m_z_all[eta_bin - 1][p_bin - 1]->Fill(dif_z, d0_m);
          h2_dif_alpha_z_all_nomw[eta_bin - 1][p_bin - 1]->Fill(dif_z, alpha);
          h2_dif_dira_z_all_nomw[eta_bin - 1][p_bin - 1]->Fill(dif_z, d0_dira);
          h2_dif_track_chi2ndof_z_all_nomw[eta_bin - 1][p_bin - 1]->Fill(
              dif_z, probe_track_chi2ndof);
          h2_dif_had_p_z_all_nomw[eta_bin - 1][p_bin - 1]->Fill(dif_z,
                                                                h_true_p);
          h2_dif_mu_p_z_all_nomw[eta_bin - 1][p_bin - 1]->Fill(dif_z,
                                                               mu_true_p);

          if (in_range(D0_M_min, d0_m, D0_M_max)) {
            h2_dif_alpha_z_all[eta_bin - 1][p_bin - 1]->Fill(dif_z, alpha);
            h2_dif_dira_z_all[eta_bin - 1][p_bin - 1]->Fill(dif_z, d0_dira);
            h2_dif_track_chi2ndof_z_all[eta_bin - 1][p_bin - 1]->Fill(
                dif_z, probe_track_chi2ndof);
            h2_dif_had_p_z_all[eta_bin - 1][p_bin - 1]->Fill(dif_z, h_true_p);
            h2_dif_mu_p_z_all[eta_bin - 1][p_bin - 1]->Fill(dif_z, mu_true_p);
          }
          if (probe_ismuon && probe_dlle < 1.0) {
            h2_dif_d0m_z_ismuon[eta_bin - 1][p_bin - 1]->Fill(dif_z, d0_m);
            if (in_range(D0_M_min, d0_m, D0_M_max)) {
              h2_dif_alpha_z_ismuon[eta_bin - 1][p_bin - 1]->Fill(dif_z, alpha);
              h2_dif_dira_z_ismuon[eta_bin - 1][p_bin - 1]->Fill(dif_z,
                                                                 d0_dira);
              h2_dif_track_chi2ndof_z_ismuon[eta_bin - 1][p_bin - 1]->Fill(
                  dif_z, probe_track_chi2ndof);
              h2_dif_had_p_z_ismuon[eta_bin - 1][p_bin - 1]->Fill(dif_z,
                                                                  h_true_p);
              h2_dif_mu_p_z_ismuon[eta_bin - 1][p_bin - 1]->Fill(dif_z,
                                                                 mu_true_p);
            }
            if (probe_dllmu > 2.0) {
              h2_dif_d0m_z_ismuon_dllmu[eta_bin - 1][p_bin - 1]->Fill(dif_z,
                                                                      d0_m);
              if (in_range(D0_M_min, d0_m, D0_M_max)) {
                h2_dif_alpha_z_ismuon_dllmu[eta_bin - 1][p_bin - 1]->Fill(
                    dif_z, alpha);
                h2_dif_dira_z_ismuon_dllmu[eta_bin - 1][p_bin - 1]->Fill(
                    dif_z, d0_dira);
                h2_dif_track_chi2ndof_z_ismuon_dllmu[eta_bin - 1][p_bin - 1]
                    ->Fill(dif_z, probe_track_chi2ndof);
                h2_dif_had_p_z_ismuon_dllmu[eta_bin - 1][p_bin - 1]->Fill(
                    dif_z, h_true_p);
                h2_dif_mu_p_z_ismuon_dllmu[eta_bin - 1][p_bin - 1]->Fill(
                    dif_z, mu_true_p);
              }
            }
            if (probe_mu_ubdt > 0.25) {
              h2_dif_d0m_z_ismuon_mubdt[eta_bin - 1][p_bin - 1]->Fill(dif_z,
                                                                      d0_m);
              if (in_range(D0_M_min, d0_m, D0_M_max)) {
                h2_dif_alpha_z_ismuon_mubdt[eta_bin - 1][p_bin - 1]->Fill(
                    dif_z, alpha);
                h2_dif_dira_z_ismuon_mubdt[eta_bin - 1][p_bin - 1]->Fill(
                    dif_z, d0_dira);
                h2_dif_track_chi2ndof_z_ismuon_mubdt[eta_bin - 1][p_bin - 1]
                    ->Fill(dif_z, probe_track_chi2ndof);
                h2_dif_had_p_z_ismuon_mubdt[eta_bin - 1][p_bin - 1]->Fill(
                    dif_z, h_true_p);
                h2_dif_mu_p_z_ismuon_mubdt[eta_bin - 1][p_bin - 1]->Fill(
                    dif_z, mu_true_p);
              }
            }
            if (probe_dllmu > 2.0 && probe_mu_ubdt > 0.25) {
              h2_dif_d0m_z_ismuon_dllmu_mubdt[eta_bin - 1][p_bin - 1]->Fill(
                  dif_z, d0_m);
              if (in_range(D0_M_min, d0_m, D0_M_max)) {
                h2_dif_alpha_z_ismuon_dllmu_mubdt[eta_bin - 1][p_bin - 1]->Fill(
                    dif_z, alpha);
                h2_dif_dira_z_ismuon_dllmu_mubdt[eta_bin - 1][p_bin - 1]->Fill(
                    dif_z, d0_dira);
                h2_dif_track_chi2ndof_z_ismuon_dllmu_mubdt[eta_bin - 1][p_bin -
                                                                        1]
                    ->Fill(dif_z, probe_track_chi2ndof);
                h2_dif_had_p_z_ismuon_dllmu_mubdt[eta_bin - 1][p_bin - 1]->Fill(
                    dif_z, h_true_p);
                h2_dif_mu_p_z_ismuon_dllmu_mubdt[eta_bin - 1][p_bin - 1]->Fill(
                    dif_z, mu_true_p);
              }
            }
          }
        }

        const bool in_pidcalib_window =
            in_range(D0_M_min, d0_m, D0_M_max) && in_range(DM_min, dm, DM_max);
        if (!in_pidcalib_window) continue;

        // IMPORTANT: The FAKE_MU pid requirement is !isMuon, but here we flip
        // it and calculate the complementary efficiency so that the K/pi
        // misid case (with less stats) falls in the "passed" sub-sample. The
        // proper efficiency is calculated afterwards as 1 - flipped_eff.
        const bool pid_ok = check_mu_pid(probe_ismuon, probe_dllmu, probe_dlle,
                                         probe_mu_ubdt, sample, run2ang);

        if (pid_ok) {
          // Fill "passed" samples
          if (dif) {
            d0_m_var.setVal(d0_m);
            dm_var.setVal(dm);
            datasets_mc_passed_dif[ntrks_bin - 1][eta_bin - 1][p_bin - 1]
                ->addFast(fit_vars);

          } else {
            d0_m_var.setVal(d0_m);
            dm_var.setVal(dm);
            datasets_mc_passed_nondif[ntrks_bin - 1][eta_bin - 1][p_bin - 1]
                ->addFast(fit_vars);
          }
        } else {
          // Fill "failed" samples
          if (dif) {
            d0_m_var.setVal(d0_m);
            dm_var.setVal(dm);
            datasets_mc_failed_dif[ntrks_bin - 1][eta_bin - 1][p_bin - 1]
                ->addFast(fit_vars);

          } else {
            d0_m_var.setVal(d0_m);
            dm_var.setVal(dm);
            datasets_mc_failed_nondif[ntrks_bin - 1][eta_bin - 1][p_bin - 1]
                ->addFast(fit_vars);
          }
        }

        // Fill histograms
        if (dif) {
          // Fill DiF z distributions
          h_dif_z_all[eta_bin - 1][p_bin - 1]->Fill(dif_z);
          h_mu_p_all[eta_bin - 1][p_bin - 1]->Fill(mu_true_p);
          if (probe_ismuon && probe_dlle < 1.0) {
            h_dif_z_ismuon[eta_bin - 1][p_bin - 1]->Fill(dif_z);
            h_mu_p_ismuon[eta_bin - 1][p_bin - 1]->Fill(mu_true_p);
            if (probe_dllmu > 2.0) {
              h_dif_z_ismuon_dllmu[eta_bin - 1][p_bin - 1]->Fill(dif_z);
              h_mu_p_ismuon_dllmu[eta_bin - 1][p_bin - 1]->Fill(mu_true_p);
            }
            if (probe_mu_ubdt > 0.25) {
              h_dif_z_ismuon_mubdt[eta_bin - 1][p_bin - 1]->Fill(dif_z);
              h_mu_p_ismuon_mubdt[eta_bin - 1][p_bin - 1]->Fill(mu_true_p);
            }
            if (probe_dllmu > 2.0 && probe_mu_ubdt > 0.25) {
              h_dif_z_ismuon_dllmu_mubdt[eta_bin - 1][p_bin - 1]->Fill(dif_z);
              h_mu_p_ismuon_dllmu_mubdt[eta_bin - 1][p_bin - 1]->Fill(
                  mu_true_p);
            }
          }

          // Now fill m(D0) distributions separating according to z DiF region
          if (!after_rich1(dif_z)) {
            h_dif_d0m_hadron_norich[eta_bin - 1][p_bin - 1]->Fill(d0_m);
          } else if (before_t3(dif_z)) {
            h_dif_d0m_hadron_rich1[eta_bin - 1][p_bin - 1]->Fill(d0_m);
          } else {
            h_dif_d0m_hadron_rich12[eta_bin - 1][p_bin - 1]->Fill(d0_m);
          }

        } else {
          h_nondif_d0m[eta_bin - 1][p_bin - 1]->Fill(d0_m);
        }
      }
    }

    // Make plots for probe
    for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
      for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
        for (int ntrks_idx = 0; ntrks_idx < N_BINS_NTRACKS; ntrks_idx++) {
          const TString suffix = TString::Format("%s_%d_%d_%d", probe.c_str(),
                                                 ntrks_idx, eta_idx, p_idx);
          const TString tag    = TString::Format(
              "(%.0f < nTracks < %.0f, %.1f < #eta < %.1f, %.0f < p < %.0f)",
              BINS_NTRACKS[ntrks_idx], BINS_NTRACKS[ntrks_idx + 1],
              BINS_ETA[eta_idx], BINS_ETA[eta_idx + 1], BINS_P[p_idx],
              BINS_P[p_idx + 1]);

          const auto &ds_mc_passed_nondif =
              datasets_mc_passed_nondif[ntrks_idx][eta_idx][p_idx];
          const auto &ds_mc_failed_nondif =
              datasets_mc_failed_nondif[ntrks_idx][eta_idx][p_idx];
          const auto &ds_mc_passed_dif =
              datasets_mc_passed_dif[ntrks_idx][eta_idx][p_idx];
          const auto &ds_mc_failed_dif =
              datasets_mc_failed_dif[ntrks_idx][eta_idx][p_idx];

          RooDataSet ds_mc_passed(
              "ds_mc_passed_" + suffix, "ds_mc_passed_" + suffix, fit_vars,
              Import(*ds_mc_passed_nondif), Import(*ds_mc_passed_dif));
          RooDataSet ds_mc_failed(
              "ds_mc_failed_" + suffix, "ds_mc_failed_" + suffix, fit_vars,
              Import(*ds_mc_failed_nondif), Import(*ds_mc_failed_dif));

          c.cd();

          const TString bin_dir_path = opath + "/" + suffix;
          gSystem->mkdir(bin_dir_path);

          // D0 M
          unique_ptr<TH1D> h_d0_m_failed(
              static_cast<TH1D *>(ds_mc_failed.createHistogram(
                  "h_d0_m_failed", d0_m_var, Binning(40, D0_M_min, D0_M_max))));
          h_d0_m_failed->SetTitle(tag);
          h_d0_m_failed->SetMarkerColor(kBlack);
          h_d0_m_failed->SetLineColor(kBlack);
          h_d0_m_failed->Scale(1. / h_d0_m_failed->GetEntries());
          unique_ptr<TH1D> h_d0_m_passed(
              static_cast<TH1D *>(ds_mc_passed.createHistogram(
                  "h_d0_m_passed", d0_m_var, Binning(40, D0_M_min, D0_M_max))));
          h_d0_m_passed->SetTitle(tag);
          h_d0_m_passed->SetMarkerColor(kRed);
          h_d0_m_passed->SetLineColor(kRed);
          h_d0_m_passed->Scale(1. / h_d0_m_passed->GetEntries());

          h_d0_m_failed->Draw();
          h_d0_m_passed->Draw("SAME");
          c.SaveAs(bin_dir_path + "/h_d0_m_" + suffix + ".pdf");

          // dm
          unique_ptr<TH1D> h_dm_failed(
              static_cast<TH1D *>(ds_mc_failed.createHistogram(
                  "h_dm_failed", dm_var, Binning(40, DM_min, DM_max))));
          h_dm_failed->SetTitle(tag);
          h_dm_failed->SetMarkerColor(kBlack);
          h_dm_failed->SetLineColor(kBlack);
          h_dm_failed->Scale(1. / h_dm_failed->GetEntries());
          unique_ptr<TH1D> h_dm_passed(
              static_cast<TH1D *>(ds_mc_passed.createHistogram(
                  "h_dm_passed", dm_var, Binning(40, DM_min, DM_max))));
          h_dm_passed->SetTitle(tag);
          h_dm_passed->SetMarkerColor(kRed);
          h_dm_passed->SetLineColor(kRed);
          h_dm_passed->Scale(1. / h_dm_passed->GetEntries());

          h_dm_failed->Draw();
          h_dm_passed->Draw("SAME");
          c.SaveAs(bin_dir_path + "/h_dm_" + suffix + ".pdf");

          // Now plot histograms with DiF z distributions
          vector<TLine *> sd_lines;
          vector<TText *> sd_labels;

          const auto &h_all_nocuts   = h_dif_z_all_nocuts[eta_idx][p_idx];
          const auto &h_all          = h_dif_z_all[eta_idx][p_idx];
          const auto &h_ismuon       = h_dif_z_ismuon[eta_idx][p_idx];
          const auto &h_ismuon_dllmu = h_dif_z_ismuon_dllmu[eta_idx][p_idx];
          const auto &h_ismuon_mubdt = h_dif_z_ismuon_mubdt[eta_idx][p_idx];
          const auto &h_ismuon_dllmu_mubdt =
              h_dif_z_ismuon_dllmu_mubdt[eta_idx][p_idx];

          h_all_nocuts->Draw("HIST");
          plot_sd_lines(c, sds, sd_lines, h_all_nocuts);
          plot_sd_labels(c, sds, sd_labels, h_all_nocuts);
          c.SaveAs(bin_dir_path + "/h_dif_z_" + suffix + "_all_nocuts.pdf");

          h_all->Draw("HIST");
          plot_sd_lines(c, sds, sd_lines, h_all);
          plot_sd_labels(c, sds, sd_labels, h_all);
          c.SaveAs(bin_dir_path + "/h_dif_z_" + suffix + "_all.pdf");

          h_ismuon->Draw("HIST");
          plot_sd_lines(c, sds, sd_lines, h_ismuon);
          plot_sd_labels(c, sds, sd_labels, h_ismuon);
          c.SaveAs(bin_dir_path + "/h_dif_z_" + suffix + "_ismuon.pdf");

          h_ismuon_dllmu->Draw("HIST");
          plot_sd_lines(c, sds, sd_lines, h_ismuon_dllmu);
          plot_sd_labels(c, sds, sd_labels, h_ismuon_dllmu);
          c.SaveAs(bin_dir_path + "/h_dif_z_" + suffix + "_ismuon_dllmu.pdf");

          h_ismuon_mubdt->Draw("HIST");
          plot_sd_lines(c, sds, sd_lines, h_ismuon_mubdt);
          plot_sd_labels(c, sds, sd_labels, h_ismuon_mubdt);
          c.SaveAs(bin_dir_path + "/h_dif_z_" + suffix + "_ismuon_mubdt.pdf");

          h_ismuon_dllmu_mubdt->Draw("HIST");
          plot_sd_lines(c, sds, sd_lines, h_ismuon_dllmu_mubdt);
          plot_sd_labels(c, sds, sd_labels, h_ismuon_dllmu_mubdt);
          c.SaveAs(bin_dir_path + "/h_dif_z_" + suffix +
                   "_ismuon_dllmu_mubdt.pdf");

          TH1D *h_ratio_ismu_to_all = static_cast<TH1D *>(
              h_all->Clone("h_ratio_ismu_to_all_" + suffix));
          h_ratio_ismu_to_all->Divide(h_ismuon, h_all, 1.0, 1.0, "B");
          h_ratio_ismu_to_all->Draw("PE");
          plot_sd_lines(c, sds, sd_lines, h_ratio_ismu_to_all);
          plot_sd_labels(c, sds, sd_labels, h_ratio_ismu_to_all);
          c.SaveAs(bin_dir_path + "/h_dif_z_" + suffix +
                   "_ratio_ismu_to_all.pdf");

          TH1D *h_ratio_ismuon_dllmu_to_ismuon = static_cast<TH1D *>(
              h_all->Clone("h_ratio_ismuon_dllmu_to_ismuon_" + suffix));
          h_ratio_ismuon_dllmu_to_ismuon->Divide(h_ismuon_dllmu, h_ismuon, 1.0,
                                                 1.0, "B");
          h_ratio_ismuon_dllmu_to_ismuon->Draw("PE");
          plot_sd_lines(c, sds, sd_lines, h_ratio_ismuon_dllmu_to_ismuon);
          plot_sd_labels(c, sds, sd_labels, h_ratio_ismuon_dllmu_to_ismuon);
          c.SaveAs(bin_dir_path + "/h_dif_z_" + suffix +
                   "_ratio_ismuon_dllmu_to_ismuon.pdf");

          TH1D *h_ratio_ismuon_mubdt_to_ismuon = static_cast<TH1D *>(
              h_all->Clone("h_ratio_ismuon_mubdt_to_ismuon_" + suffix));
          h_ratio_ismuon_mubdt_to_ismuon->Divide(h_ismuon_mubdt, h_ismuon, 1.0,
                                                 1.0, "B");
          h_ratio_ismuon_mubdt_to_ismuon->Draw("PE");
          plot_sd_lines(c, sds, sd_lines, h_ratio_ismuon_mubdt_to_ismuon);
          plot_sd_labels(c, sds, sd_labels, h_ratio_ismuon_mubdt_to_ismuon);
          c.SaveAs(bin_dir_path + "/h_dif_z_" + suffix +
                   "_ratio_ismuon_mubdt_to_ismuon.pdf");

          TH1D *h_ratio_ismuon_dllmu_mubdt_to_ismuon = static_cast<TH1D *>(
              h_all->Clone("h_ratio_ismuon_dllmu_mubdt_to_ismuon_" + suffix));
          h_ratio_ismuon_dllmu_mubdt_to_ismuon->Divide(h_ismuon_dllmu_mubdt,
                                                       h_ismuon, 1.0, 1.0, "B");
          h_ratio_ismuon_dllmu_mubdt_to_ismuon->Draw("PE");
          plot_sd_lines(c, sds, sd_lines, h_ratio_ismuon_dllmu_mubdt_to_ismuon);
          plot_sd_labels(c, sds, sd_labels,
                         h_ratio_ismuon_dllmu_mubdt_to_ismuon);
          c.SaveAs(bin_dir_path + "/h_dif_z_" + suffix +
                   "_ratio_ismuon_dllmu_mubdt_to_ismuon.pdf");

          h_mu_p_all_nocuts[eta_idx][p_idx]->Draw("HIST");
          c.SaveAs(bin_dir_path + "/h_mu_p_" + suffix + "_all_nocuts.pdf");

          h_mu_p_all[eta_idx][p_idx]->Draw("HIST");
          c.SaveAs(bin_dir_path + "/h_mu_p_" + suffix + "_all.pdf");

          h_mu_p_ismuon[eta_idx][p_idx]->Draw("HIST");
          c.SaveAs(bin_dir_path + "/h_mu_p_" + suffix + "_ismuon.pdf");

          h_mu_p_ismuon_dllmu[eta_idx][p_idx]->Draw("HIST");
          c.SaveAs(bin_dir_path + "/h_mu_p_" + suffix + "_ismuon_dllmu.pdf");

          h_mu_p_ismuon_mubdt[eta_idx][p_idx]->Draw("HIST");
          c.SaveAs(bin_dir_path + "/h_mu_p_" + suffix + "_ismuon_mubdt.pdf");

          h_mu_p_ismuon_dllmu_mubdt[eta_idx][p_idx]->Draw("HIST");
          c.SaveAs(bin_dir_path + "/h_mu_p_" + suffix +
                   "_ismuon_dllmu_mubdt.pdf");

          TH1D *h_mu_p_ratio_ismu_to_all =
              static_cast<TH1D *>(h_mu_p_all[eta_idx][p_idx]->Clone(
                  "h_ratio_ismu_to_all_" + suffix));
          h_mu_p_ratio_ismu_to_all->Divide(h_mu_p_ismuon[eta_idx][p_idx],
                                           h_mu_p_all[eta_idx][p_idx], 1.0, 1.0,
                                           "B");
          h_mu_p_ratio_ismu_to_all->Draw("PE");
          c.SaveAs(bin_dir_path + "/h_mu_p_" + suffix +
                   "_ratio_ismu_to_all.pdf");

          TH1D *h_mu_p_ratio_ismuon_dllmu_to_ismuon =
              static_cast<TH1D *>(h_mu_p_all[eta_idx][p_idx]->Clone(
                  "h_ratio_ismuon_dllmu_to_ismuon_" + suffix));
          h_mu_p_ratio_ismuon_dllmu_to_ismuon->Divide(
              h_mu_p_ismuon_dllmu[eta_idx][p_idx],
              h_mu_p_ismuon[eta_idx][p_idx], 1.0, 1.0, "B");
          h_mu_p_ratio_ismuon_dllmu_to_ismuon->Draw("PE");
          c.SaveAs(bin_dir_path + "/h_mu_p_" + suffix +
                   "_ratio_ismuon_dllmu_to_ismuon.pdf");

          TH1D *h_mu_p_ratio_ismuon_mubdt_to_ismuon =
              static_cast<TH1D *>(h_mu_p_all[eta_idx][p_idx]->Clone(
                  "h_ratio_ismuon_mubdt_to_ismuon_" + suffix));
          h_mu_p_ratio_ismuon_mubdt_to_ismuon->Divide(
              h_mu_p_ismuon_mubdt[eta_idx][p_idx],
              h_mu_p_ismuon[eta_idx][p_idx], 1.0, 1.0, "B");
          h_mu_p_ratio_ismuon_mubdt_to_ismuon->Draw("PE");
          c.SaveAs(bin_dir_path + "/h_mu_p_" + suffix +
                   "_ratio_ismuon_mubdt_to_ismuon.pdf");

          TH1D *h_mu_p_ratio_ismuon_dllmu_mubdt_to_ismuon =
              static_cast<TH1D *>(h_mu_p_all[eta_idx][p_idx]->Clone(
                  "h_ratio_ismuon_dllmu_mubdt_to_ismuon_" + suffix));
          h_mu_p_ratio_ismuon_dllmu_mubdt_to_ismuon->Divide(
              h_mu_p_ismuon_dllmu_mubdt[eta_idx][p_idx],
              h_mu_p_ismuon[eta_idx][p_idx], 1.0, 1.0, "B");
          h_mu_p_ratio_ismuon_dllmu_mubdt_to_ismuon->Draw("PE");
          c.SaveAs(bin_dir_path + "/h_mu_p_" + suffix +
                   "_ratio_ismuon_dllmu_mubdt_to_ismuon.pdf");

          h_dif_d0m_hadron_rich12[eta_idx][p_idx]->SetLineColor(kBlue);
          h_dif_d0m_hadron_rich12[eta_idx][p_idx]->Draw("HIST");
          c.SaveAs(bin_dir_path + "/h_dif_d0m_rich12_" + suffix + "_all.pdf");

          h_dif_d0m_hadron_rich1[eta_idx][p_idx]->SetLineColor(kOrange);
          h_dif_d0m_hadron_rich1[eta_idx][p_idx]->Draw("HIST");
          c.SaveAs(bin_dir_path + "/h_dif_d0m_rich1_" + suffix + "_all.pdf");

          h_dif_d0m_hadron_norich[eta_idx][p_idx]->SetLineColor(kRed);
          h_dif_d0m_hadron_norich[eta_idx][p_idx]->Draw("HIST");
          c.SaveAs(bin_dir_path + "/h_dif_d0m_norich_z_" + suffix + "_all.pdf");

          h_nondif_d0m[eta_idx][p_idx]->SetLineColor(kBlack);
          h_nondif_d0m[eta_idx][p_idx]->Draw("HIST");
          c.SaveAs(bin_dir_path + "/h_nondif_d0m_" + suffix + "_all.pdf");

          h_dif_d0m_hadron_rich12[eta_idx][p_idx]->Scale(
              1.0 / h_dif_d0m_hadron_rich12[eta_idx][p_idx]->Integral());
          h_dif_d0m_hadron_rich1[eta_idx][p_idx]->Scale(
              1.0 / h_dif_d0m_hadron_rich1[eta_idx][p_idx]->Integral());
          h_dif_d0m_hadron_norich[eta_idx][p_idx]->Scale(
              1.0 / h_dif_d0m_hadron_norich[eta_idx][p_idx]->Integral());
          h_nondif_d0m[eta_idx][p_idx]->Scale(
              1.0 / h_nondif_d0m[eta_idx][p_idx]->Integral());

          h_nondif_d0m[eta_idx][p_idx]->Draw("PE");
          h_dif_d0m_hadron_rich12[eta_idx][p_idx]->Draw("PE SAME");
          h_dif_d0m_hadron_rich1[eta_idx][p_idx]->Draw("PE SAME");
          h_dif_d0m_hadron_norich[eta_idx][p_idx]->Draw("PE SAME");

          c.SaveAs(bin_dir_path + "/h_dif_d0m_rich_comp_" + suffix +
                   "_all.pdf");

          h2_dif_d0m_z_all_nocuts[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_d0m_z_all_nocuts[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_d0m_z_all_nocuts[eta_idx][p_idx]);
          plot_d0m_lines(c, sd_lines, h2_dif_d0m_z_all_nocuts[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_d0m_z_all_nocuts_" + suffix +
                   ".pdf");

          h2_dif_d0m_z_all[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines, h2_dif_d0m_z_all[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels, h2_dif_d0m_z_all[eta_idx][p_idx]);
          plot_d0m_lines(c, sd_lines, h2_dif_d0m_z_all[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_d0m_z_all_" + suffix + ".pdf");

          h2_dif_d0m_z_ismuon[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines, h2_dif_d0m_z_ismuon[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_d0m_z_ismuon[eta_idx][p_idx]);
          plot_d0m_lines(c, sd_lines, h2_dif_d0m_z_ismuon[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_d0m_z_ismuon_" + suffix + ".pdf");

          h2_dif_d0m_z_ismuon_dllmu[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_d0m_z_ismuon_dllmu[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_d0m_z_ismuon_dllmu[eta_idx][p_idx]);
          plot_d0m_lines(c, sd_lines,
                         h2_dif_d0m_z_ismuon_dllmu[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_d0m_z_ismuon_dllmu_" + suffix +
                   ".pdf");

          h2_dif_d0m_z_ismuon_mubdt[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_d0m_z_ismuon_mubdt[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_d0m_z_ismuon_mubdt[eta_idx][p_idx]);
          plot_d0m_lines(c, sd_lines,
                         h2_dif_d0m_z_ismuon_mubdt[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_d0m_z_ismuon_mubdt_" + suffix +
                   ".pdf");

          h2_dif_d0m_z_ismuon_dllmu_mubdt[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_d0m_z_ismuon_dllmu_mubdt[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_d0m_z_ismuon_dllmu_mubdt[eta_idx][p_idx]);
          plot_d0m_lines(c, sd_lines,
                         h2_dif_d0m_z_ismuon_dllmu_mubdt[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_d0m_z_ismuon_dllmu_mubdt_" + suffix +
                   ".pdf");

          h2_dif_alpha_z_all_nocuts[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_alpha_z_all_nocuts[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_alpha_z_all_nocuts[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_alpha_z_all_nocuts_" + suffix +
                   ".pdf");

          h2_dif_alpha_z_all_nomw[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_alpha_z_all_nomw[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_alpha_z_all_nomw[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_alpha_z_all_nomw_" + suffix +
                   ".pdf");

          h2_dif_alpha_z_all[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines, h2_dif_alpha_z_all[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels, h2_dif_alpha_z_all[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_alpha_z_all_" + suffix + ".pdf");

          h2_dif_alpha_z_ismuon[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_alpha_z_ismuon[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_alpha_z_ismuon[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_alpha_z_ismuon_" + suffix + ".pdf");

          h2_dif_alpha_z_ismuon_dllmu[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_alpha_z_ismuon_dllmu[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_alpha_z_ismuon_dllmu[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_alpha_z_ismuon_dllmu_" + suffix +
                   ".pdf");

          h2_dif_alpha_z_ismuon_mubdt[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_alpha_z_ismuon_mubdt[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_alpha_z_ismuon_mubdt[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_alpha_z_ismuon_mubdt_" + suffix +
                   ".pdf");

          h2_dif_alpha_z_ismuon_dllmu_mubdt[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_alpha_z_ismuon_dllmu_mubdt[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_alpha_z_ismuon_dllmu_mubdt[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_alpha_z_ismuon_dllmu_mubdt_" +
                   suffix + ".pdf");

          h2_dif_dira_z_all_nocuts[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_dira_z_all_nocuts[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_dira_z_all_nocuts[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_dira_z_all_nocuts_" + suffix +
                   ".pdf");

          h2_dif_dira_z_all_nomw[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_dira_z_all_nomw[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_dira_z_all_nomw[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_dira_z_all_nomw_" + suffix + ".pdf");

          h2_dif_dira_z_all[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines, h2_dif_dira_z_all[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels, h2_dif_dira_z_all[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_dira_z_all_" + suffix + ".pdf");

          h2_dif_dira_z_ismuon[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines, h2_dif_dira_z_ismuon[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_dira_z_ismuon[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_dira_z_ismuon_" + suffix + ".pdf");

          h2_dif_dira_z_ismuon_dllmu[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_dira_z_ismuon_dllmu[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_dira_z_ismuon_dllmu[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_dira_z_ismuon_dllmu_" + suffix +
                   ".pdf");

          h2_dif_dira_z_ismuon_mubdt[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_dira_z_ismuon_mubdt[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_dira_z_ismuon_mubdt[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_dira_z_ismuon_mubdt_" + suffix +
                   ".pdf");

          h2_dif_dira_z_ismuon_dllmu_mubdt[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_dira_z_ismuon_dllmu_mubdt[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_dira_z_ismuon_dllmu_mubdt[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_dira_z_ismuon_dllmu_mubdt_" +
                   suffix + ".pdf");

          h2_dif_track_chi2ndof_z_all_nocuts[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_track_chi2ndof_z_all_nocuts[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_track_chi2ndof_z_all_nocuts[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_track_chi2ndof_z_all_nocuts_" +
                   suffix + ".pdf");

          h2_dif_track_chi2ndof_z_all_nomw[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_track_chi2ndof_z_all_nomw[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_track_chi2ndof_z_all_nomw[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_track_chi2ndof_z_all_nomw_" +
                   suffix + ".pdf");

          h2_dif_track_chi2ndof_z_all[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_track_chi2ndof_z_all[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_track_chi2ndof_z_all[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_track_chi2ndof_z_all_" + suffix +
                   ".pdf");

          h2_dif_track_chi2ndof_z_ismuon[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_track_chi2ndof_z_ismuon[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_track_chi2ndof_z_ismuon[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_track_chi2ndof_z_ismuon_" + suffix +
                   ".pdf");

          h2_dif_track_chi2ndof_z_ismuon_dllmu[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_track_chi2ndof_z_ismuon_dllmu[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_track_chi2ndof_z_ismuon_dllmu[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_track_chi2ndof_z_ismuon_dllmu_" +
                   suffix + ".pdf");

          h2_dif_track_chi2ndof_z_ismuon_mubdt[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_track_chi2ndof_z_ismuon_mubdt[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_track_chi2ndof_z_ismuon_mubdt[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_track_chi2ndof_z_ismuon_mubdt_" +
                   suffix + ".pdf");

          h2_dif_track_chi2ndof_z_ismuon_dllmu_mubdt[eta_idx][p_idx]->Draw(
              "COLZ");
          plot_sd_lines(
              c, sds, sd_lines,
              h2_dif_track_chi2ndof_z_ismuon_dllmu_mubdt[eta_idx][p_idx]);
          plot_sd_labels(
              c, sds, sd_labels,
              h2_dif_track_chi2ndof_z_ismuon_dllmu_mubdt[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path +
                   "/h2_dif_track_chi2ndof_z_ismuon_dllmu_mubdt_" + suffix +
                   ".pdf");

          h2_dif_had_p_z_all_nocuts[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_had_p_z_all_nocuts[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_had_p_z_all_nocuts[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_had_p_z_all_nocuts_" + suffix +
                   ".pdf");

          h2_dif_had_p_z_all_nomw[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_had_p_z_all_nomw[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_had_p_z_all_nomw[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_had_p_z_all_nomw_" + suffix +
                   ".pdf");

          h2_dif_had_p_z_all[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines, h2_dif_had_p_z_all[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels, h2_dif_had_p_z_all[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_had_p_z_all_" + suffix + ".pdf");

          h2_dif_had_p_z_ismuon[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_had_p_z_ismuon[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_had_p_z_ismuon[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_had_p_z_ismuon_" + suffix + ".pdf");

          h2_dif_had_p_z_ismuon_dllmu[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_had_p_z_ismuon_dllmu[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_had_p_z_ismuon_dllmu[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_had_p_z_ismuon_dllmu_" + suffix +
                   ".pdf");

          h2_dif_had_p_z_ismuon_mubdt[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_had_p_z_ismuon_mubdt[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_had_p_z_ismuon_mubdt[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_had_p_z_ismuon_mubdt_" + suffix +
                   ".pdf");

          h2_dif_had_p_z_ismuon_dllmu_mubdt[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_had_p_z_ismuon_dllmu_mubdt[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_had_p_z_ismuon_dllmu_mubdt[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_had_p_z_ismuon_dllmu_mubdt_" +
                   suffix + ".pdf");

          h2_dif_mu_p_z_all_nocuts[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_mu_p_z_all_nocuts[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_mu_p_z_all_nocuts[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_mu_p_z_all_nocuts_" + suffix +
                   ".pdf");

          h2_dif_mu_p_z_all_nomw[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_mu_p_z_all_nomw[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_mu_p_z_all_nomw[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_mu_p_z_all_nomw_" + suffix + ".pdf");

          h2_dif_mu_p_z_all[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines, h2_dif_mu_p_z_all[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels, h2_dif_mu_p_z_all[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_mu_p_z_all_" + suffix + ".pdf");

          h2_dif_mu_p_z_ismuon[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines, h2_dif_mu_p_z_ismuon[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_mu_p_z_ismuon[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_mu_p_z_ismuon_" + suffix + ".pdf");

          h2_dif_mu_p_z_ismuon_dllmu[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_mu_p_z_ismuon_dllmu[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_mu_p_z_ismuon_dllmu[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_mu_p_z_ismuon_dllmu_" + suffix +
                   ".pdf");

          h2_dif_mu_p_z_ismuon_mubdt[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_mu_p_z_ismuon_mubdt[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_mu_p_z_ismuon_mubdt[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_mu_p_z_ismuon_mubdt_" + suffix +
                   ".pdf");

          h2_dif_mu_p_z_ismuon_dllmu_mubdt[eta_idx][p_idx]->Draw("COLZ");
          plot_sd_lines(c, sds, sd_lines,
                        h2_dif_mu_p_z_ismuon_dllmu_mubdt[eta_idx][p_idx]);
          plot_sd_labels(c, sds, sd_labels,
                         h2_dif_mu_p_z_ismuon_dllmu_mubdt[eta_idx][p_idx]);
          c.SaveAs(bin_dir_path + "/h2_dif_mu_p_z_ismuon_dllmu_mubdt_" +
                   suffix + ".pdf");
        }
      }
    }
  }
}