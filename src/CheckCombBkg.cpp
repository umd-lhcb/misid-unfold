// Author: Lucas Meyer Garcia
// License: BSD 2-clause
//
// Description: Plot m(D0) distributions from wrong-sign samples
// [D*+ -> D0 (-> K+ pi+) pi+]cc built from the NoBias stream to check how the
// PIDCalib kinematic cuts affect the shape of combinatorial background bin by
// bin.

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <Math/Vector4D.h>
#include "TCanvas.h"
#include "TChain.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TH2D.h"
#include "TString.h"

#include "RooArgSet.h"
#include "RooBinning.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooFitResult.h"
#include "RooHelpers.h"
#include "RooKeysPdf.h"  // Customized version in root-curated
#include "RooPlot.h"
#include "RooPowerLaw.h"  // Custom PDF defined in root-curated
#include "RooProdPdf.h"
#include "RooRealVar.h"

#include <yaml-cpp/yaml.h>
#include <cxxopts.hpp>

#include "misid.h"
#include "utils.h"

using std::cout, std::endl, std::right, std::left, std::setw, std::setprecision,
    std::fixed, std::to_string;
using std::string, std::unique_ptr, std::unordered_map, std::vector, std::map;

using ROOT::Math::PxPyPzMVector;

int main(int argc, char **argv) {
  cxxopts::Options argOpts(
      "CheckCombBkg",
      "Check shape of comb bkg in pidcalib samples using WS D* samples.");

  // clang-format off
  argOpts.add_options()
    ("h,help", "Print help")
    ("d,debug", "Enable debug mode",
     cxxopts::value<bool>()->default_value("false"))
    ("c,config", "Specify input YAML config file",
     cxxopts::value<string>())
    ("o,output", "Specify output folder",
     cxxopts::value<string>()->default_value("gen/"))
    ;
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  const auto ymlFile   = parsedArgs["config"].as<string>();
  const auto ymlConfig = YAML::LoadFile(ymlFile)["misid_corrections"];

  const string opath = parsedArgs["output"].as<string>();
  cout << "INFO Output will be saved in " << opath << endl;

  int    ntracks, dst_vtx_ndof, d0_vtx_ndof;
  double dst_m, d0_m, probe_p, probe_pz, tag_pz, k_px, k_py, pi_px, pi_py,
      k_track_chi2ndof, pi_track_chi2ndof, spi_track_chi2ndof, probe_mudll,
      k_ghostprob, pi_ghostprob, spi_ghostprob, dst_vtx_chi2, d0_vtx_chi2;
  bool probe_ismuon, probe_hasmuon;

  const vector<string> particles = {"k", "pi"};

  TH2D histo_binning("histo_binning", ";p;eta", N_BINS_P, BINS_P, N_BINS_ETA,
                     BINS_ETA);

  // Calculate histogram ranges with additional bins
  // Limits are the same for all samples

  constexpr double extended_d0_m_max = 1.920;
  constexpr double extended_d0_m_min = 1.810;
  constexpr double extended_dm_max   = DM_max;
  constexpr double extended_dm_min   = DM_min;

  RooBinning bins_histos_d0_m(55, extended_d0_m_min, extended_d0_m_max,
                              "bins_histos_d0_m");
  RooBinning bins_histos_dm(60, extended_dm_min, extended_dm_max,
                            "bins_histos_dm");

  // Observables
  RooRealVar d0_m_var("d0_m_var", "d0_m_var", extended_d0_m_min,
                      extended_d0_m_max, "GeV");
  RooRealVar dm_var("dm_var", "dm_var", extended_dm_min, extended_dm_max,
                    "GeV");

  RooArgSet fit_vars(d0_m_var, dm_var);

  RooRealVar     k_comb("k_comb", "k_comb", 0., -10., 10., "");
  RooExponential d0_comb("d0_comb", "d0_comb", d0_m_var, k_comb);

  RooConstVar pi_m_pdg("pi_m_pdg", "pi_m_pdg", PI_M);  // GeV

  RooRealVar  c_comb("c_comb", "c_comb", 0., 1., "");
  RooPowerLaw dm_comb("dm_comb", "dm_comb", dm_var, pi_m_pdg, c_comb);

  TCanvas c_single("c_single", "c_single", 960, 720);
  TCanvas c_double("c_double", "c_double", 1920, 720);
  c_double.Divide(2, 1);

  for (auto &probe : particles) {
    const TString tag = (probe == "pi") ? "k" : "pi";

    array<array<RooDataSet *, N_BINS_P>, N_BINS_ETA> datasets;
    array<array<RooDataSet *, N_BINS_P>, N_BINS_ETA> datasets_wmcut;
    array<array<RooDataSet *, N_BINS_P>, N_BINS_ETA> datasets_wmcut_offline;
    array<array<RooDataSet *, N_BINS_P>, N_BINS_ETA> datasets_wmcut_mupid;

    // Initialize MC datasets
    cout << "INFO Initializing signal MC datahists" << endl;
    for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
      for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
        TString suffix =
            TString::Format("%s_%d_%d", probe.c_str(), eta_idx, p_idx);
        datasets[eta_idx][p_idx] =
            new RooDataSet("ds_" + suffix, "ds_" + suffix, fit_vars);
        datasets_wmcut[eta_idx][p_idx] = new RooDataSet(
            "ds_wmcut_" + suffix, "ds_wmcut_" + suffix, fit_vars);
        datasets_wmcut_offline[eta_idx][p_idx] =
            new RooDataSet("ds_wmcut_offline_" + suffix,
                           "ds_wmcut_offline_" + suffix, fit_vars);
        datasets_wmcut_mupid[eta_idx][p_idx] = new RooDataSet(
            "ds_wmcut_mupid_" + suffix, "ds_wmcut_mupid_" + suffix, fit_vars);
      }
    }

    // Open and loop over data files
    const auto    mc_path = ymlConfig["mc_ntps"]["nobias"][probe].as<string>();
    const TString tdir_wmcut =
        (probe == "k") ? "TupleDstANNK_OS_WMcut" : "TupleDstANNPi_OS_WMcut";
    TChain ch_wmcut(tdir_wmcut + "/DecayTree");

    cout << "INFO Opening files: " << mc_path << endl;
    ch_wmcut.Add(mc_path.c_str());

    cout << "INFO Opened files:" << endl;
    print_files(ch_wmcut);

    ch_wmcut.SetBranchStatus("*", false);
    ch_wmcut.SetBranchAddress("dst_M", &dst_m);
    ch_wmcut.SetBranchAddress("dst_ENDVERTEX_CHI2", &dst_vtx_chi2);
    ch_wmcut.SetBranchAddress("dst_ENDVERTEX_NDOF", &dst_vtx_ndof);
    ch_wmcut.SetBranchAddress("d0_M", &d0_m);
    ch_wmcut.SetBranchAddress("d0_ENDVERTEX_CHI2", &d0_vtx_chi2);
    ch_wmcut.SetBranchAddress("d0_ENDVERTEX_NDOF", &d0_vtx_ndof);
    ch_wmcut.SetBranchAddress((probe + "_P").c_str(), &probe_p);
    ch_wmcut.SetBranchAddress((probe + "_PZ").c_str(), &probe_pz);
    ch_wmcut.SetBranchAddress((probe + "_isMuon").c_str(), &probe_ismuon);
    ch_wmcut.SetBranchAddress((probe + "_hasMuon").c_str(), &probe_hasmuon);
    ch_wmcut.SetBranchAddress((probe + "_PIDmu").c_str(), &probe_mudll);
    ch_wmcut.SetBranchAddress("k_TRACK_CHI2NDOF", &k_track_chi2ndof);
    ch_wmcut.SetBranchAddress("pi_TRACK_CHI2NDOF", &pi_track_chi2ndof);
    ch_wmcut.SetBranchAddress("spi_TRACK_CHI2NDOF", &spi_track_chi2ndof);
    ch_wmcut.SetBranchAddress("k_TRACK_GhostProb", &k_ghostprob);
    ch_wmcut.SetBranchAddress("pi_TRACK_GhostProb", &pi_ghostprob);
    ch_wmcut.SetBranchAddress("spi_TRACK_GhostProb", &spi_ghostprob);
    ch_wmcut.SetBranchAddress("nTracks", &ntracks);

    const int entries_wmcut = ch_wmcut.GetEntries();
    cout << "INFO Starting MC event loop over " << entries_wmcut << " entries"
         << endl;

    int eta_bin = 0, p_bin = 0, dummy = 0;

    int count_a = 0, count_b = 0, count_c = 0, count_d = 0, count_e = 0,
        count_f = 0;

    for (int evt = 0; evt < entries_wmcut; evt++) {
      ch_wmcut.GetEntry(evt);

      // Conditional cuts
      if (!probe_hasmuon) continue;
      count_a++;

      if ((k_ghostprob > 0.4) || (pi_ghostprob > 0.4) || (spi_ghostprob > 0.4))
        continue;
      count_b++;

      if ((k_track_chi2ndof > 3.) || (pi_track_chi2ndof > 3.) ||
          (spi_track_chi2ndof > 3.))
        continue;
      count_c++;

      if ((dst_vtx_chi2 / dst_vtx_ndof) > 15 ||
          (d0_vtx_chi2 / d0_vtx_ndof) > 10)
        continue;
      count_d++;

      if (probe_p < 3000. || probe_p >= 100000. || ntracks >= 600) continue;
      count_e++;

      const double probe_eta =
          0.5 * std::log((probe_p + probe_pz) / (probe_p - probe_pz));
      if (probe_eta < 1.7 || probe_eta >= 5.0) continue;

      // Determine kinematical bin
      const int kin_bin = histo_binning.FindBin(probe_p, probe_eta, ntracks);
      histo_binning.GetBinXYZ(kin_bin, p_bin, eta_bin, dummy);

      // Fill histograms
      const double dm = (dst_m - d0_m) * 0.001;
      d0_m            = d0_m * 0.001;

      const bool in_extended_window =
          in_var_range(d0_m_var, d0_m) && in_var_range(dm_var, dm);
      if (!in_extended_window) continue;

      count_f++;
      d0_m_var.setVal(d0_m);
      dm_var.setVal(dm);
      datasets_wmcut[eta_bin - 1][p_bin - 1]->addFast(fit_vars);
      if (probe_ismuon && probe_mudll > 2.)
        datasets_wmcut_mupid[eta_bin - 1][p_bin - 1]->addFast(fit_vars);
    }

    cout << "INFO cutflow:\n";
    cout << " - Total MC events:              " << entries_wmcut << "\n";
    cout << " - After a: " << count_a << "(" << count_a * 100. / entries_wmcut
         << "%)\n";
    cout << " - After b: " << count_b << "(" << count_b * 100. / count_a
         << "%)\n";
    cout << " - After c: " << count_c << "(" << count_c * 100. / count_b
         << "%)\n";
    cout << " - After d: " << count_d << "(" << count_d * 100. / count_c
         << "%)\n";
    cout << " - After e: " << count_e << "(" << count_e * 100. / count_d
         << "%)\n";
    cout << " - After f: " << count_f << "(" << count_f * 100. / count_e
         << "%)\n";
    cout << " Overall: " << "(" << count_f * 100. / entries_wmcut << "%)\n"
         << endl;

    const TString tdir =
        (probe == "k") ? "TupleDstANNK_OS" : "TupleDstANNPi_OS";
    TChain ch(tdir + "/DecayTree");

    cout << "INFO Opening files: " << mc_path << endl;
    ch.Add(mc_path.c_str());

    cout << "INFO Opened files:" << endl;
    print_files(ch);

    ch.SetBranchStatus("*", false);
    ch.SetBranchAddress("dst_M", &dst_m);
    ch.SetBranchAddress("dst_ENDVERTEX_CHI2", &dst_vtx_chi2);
    ch.SetBranchAddress("dst_ENDVERTEX_NDOF", &dst_vtx_ndof);
    ch.SetBranchAddress("d0_M", &d0_m);
    ch.SetBranchAddress("d0_ENDVERTEX_CHI2", &d0_vtx_chi2);
    ch.SetBranchAddress("d0_ENDVERTEX_NDOF", &d0_vtx_ndof);
    ch.SetBranchAddress((probe + "_P").c_str(), &probe_p);
    ch.SetBranchAddress((probe + "_PZ").c_str(), &probe_pz);
    ch.SetBranchAddress((probe + "_hasMuon").c_str(), &probe_hasmuon);
    ch.SetBranchAddress(tag + "_PZ", &tag_pz);
    ch.SetBranchAddress("k_PX", &k_px);
    ch.SetBranchAddress("k_PY", &k_py);
    ch.SetBranchAddress("pi_PX", &pi_px);
    ch.SetBranchAddress("pi_PY", &pi_py);
    ch.SetBranchAddress("k_TRACK_CHI2NDOF", &k_track_chi2ndof);
    ch.SetBranchAddress("pi_TRACK_CHI2NDOF", &pi_track_chi2ndof);
    ch.SetBranchAddress("spi_TRACK_CHI2NDOF", &spi_track_chi2ndof);
    ch.SetBranchAddress("k_TRACK_GhostProb", &k_ghostprob);
    ch.SetBranchAddress("pi_TRACK_GhostProb", &pi_ghostprob);
    ch.SetBranchAddress("spi_TRACK_GhostProb", &spi_ghostprob);
    ch.SetBranchAddress("nTracks", &ntracks);

    const int entries = ch.GetEntries();
    cout << "INFO Starting event loop over " << entries << " entries" << endl;

    count_a = 0;
    count_b = 0;
    count_c = 0;
    count_d = 0;
    count_e = 0;
    count_f = 0;

    for (int evt = 0; evt < entries; evt++) {
      ch.GetEntry(evt);

      // Conditional cuts
      if (!probe_hasmuon) continue;
      count_a++;

      if ((k_ghostprob > 0.4) || (pi_ghostprob > 0.4) || (spi_ghostprob > 0.4))
        continue;
      count_b++;

      if ((k_track_chi2ndof > 3.) || (pi_track_chi2ndof > 3.) ||
          (spi_track_chi2ndof > 3.))
        continue;
      count_c++;

      if ((dst_vtx_chi2 / dst_vtx_ndof) > 15 ||
          (d0_vtx_chi2 / d0_vtx_ndof) > 10)
        continue;
      count_d++;

      if (probe_p < 3000. || probe_p >= 100000. || ntracks >= 600) continue;
      count_e++;

      const double probe_eta =
          0.5 * std::log((probe_p + probe_pz) / (probe_p - probe_pz));
      if (probe_eta < 1.7 || probe_eta >= 5.0) continue;

      // Determine kinematical bin
      const int kin_bin = histo_binning.FindBin(probe_p, probe_eta, ntracks);
      histo_binning.GetBinXYZ(kin_bin, p_bin, eta_bin, dummy);

      // Fill histograms
      const double dm = (dst_m - d0_m) * 0.001;
      d0_m            = d0_m * 0.001;

      const bool in_extended_window =
          in_var_range(d0_m_var, d0_m) && in_var_range(dm_var, dm);
      if (!in_extended_window) continue;

      count_f++;
      d0_m_var.setVal(d0_m);
      dm_var.setVal(dm);
      datasets[eta_bin - 1][p_bin - 1]->addFast(fit_vars);

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

      datasets_wmcut_offline[eta_bin - 1][p_bin - 1]->addFast(fit_vars);
    }

    cout << "INFO cutflow:\n";
    cout << " - Total MC events:              " << entries << "\n";
    cout << " - After a: " << count_a << "(" << count_a * 100. / entries
         << "%)\n";
    cout << " - After b: " << count_b << "(" << count_b * 100. / count_a
         << "%)\n";
    cout << " - After c: " << count_c << "(" << count_c * 100. / count_b
         << "%)\n";
    cout << " - After d: " << count_d << "(" << count_d * 100. / count_c
         << "%)\n";
    cout << " - After e: " << count_e << "(" << count_e * 100. / count_d
         << "%)\n";
    cout << " - After f: " << count_f << "(" << count_f * 100. / count_e
         << "%)\n";
    cout << " Overall: " << "(" << count_f * 100. / entries << "%)\n" << endl;

    d0_m_var.setRange("pidcalib", D0_M_min, D0_M_max);
    dm_var.setRange("pidcalib", DM_min, DM_max);

    for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
      for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
        const TString tag = TString::Format(
            "(%.1f < #eta < %.1f, %.0f < p < %.0f)", BINS_ETA[eta_idx],
            BINS_ETA[eta_idx + 1], BINS_P[p_idx], BINS_P[p_idx + 1]);
        const TString suffix =
            TString::Format("%s_%d_%d", probe.c_str(), eta_idx, p_idx);

        const auto &dataset_wmcut_mupid = datasets_wmcut_mupid[eta_idx][p_idx];
        const auto &dataset_wmcut       = datasets_wmcut[eta_idx][p_idx];
        const auto &dataset_wmcut_offline =
            datasets_wmcut_offline[eta_idx][p_idx];
        const auto &dataset = datasets[eta_idx][p_idx];

        d0_m_var.setRange("fitRange", 1.825, 1.910);
        dm_var.setRange("fitRange", 0.141, 0.153);

        // Cut regions where combinatorial bkg deviates from exponential
        // distribution
        if (probe == "k") {
          if (p_idx == 5) {
            d0_m_var.setRange("fitRange", 1.825, 1.900);
          }
          if (p_idx == 4) {
            d0_m_var.setRange("fitRange", 1.825, 1.902);
          }
          if (p_idx == 3 && eta_idx == 0) {
            d0_m_var.setRange("fitRange", 1.825, 1.906);
          }
        }
        if (probe == "pi") {
          if (p_idx == 2) {
            d0_m_var.setRange("fitRange", 1.825, 1.906);
          }
          if (p_idx == 1) {
            d0_m_var.setRange("fitRange", 1.825, 1.904);
          }
          if (p_idx == 0) {
            d0_m_var.setRange("fitRange", 1.825, 1.902);
          }
        }

        // Fit sample without WM cut
        d0_comb.fitTo(*dataset, Strategy(2), Range("fitRange"), PrintLevel(0));

        unique_ptr<RooPlot> frame_d0(
            d0_m_var.frame(Title("D0 M - No WM cut " + tag)));

        dataset->plotOn(frame_d0.get(), Binning(bins_histos_d0_m));
        d0_comb.plotOn(frame_d0.get(), LineWidth(2), LineColor(kBlue),
                       Range(""), NormRange("fitRange"));
        d0_comb.plotOn(frame_d0.get(), LineWidth(2), LineColor(kMagenta),
                       VLines(), Range("pidcalib"), NormRange("fitRange"));
        d0_comb.plotOn(frame_d0.get(), LineWidth(2), LineColor(kRed), VLines(),
                       LineStyle(kDashed), Range("fitRange"));

        // Fit sample with WM cut
        d0_comb.fitTo(*dataset_wmcut, Strategy(2), Range("fitRange"),
                      PrintLevel(0));

        unique_ptr<RooPlot> frame_d0_wmcut(
            d0_m_var.frame(Title("D0 M - With WM cut " + tag)));

        dataset_wmcut->plotOn(frame_d0_wmcut.get(), Binning(bins_histos_d0_m));
        d0_comb.plotOn(frame_d0_wmcut.get(), LineWidth(2), LineColor(kBlue),
                       Range(""), NormRange("fitRange"));
        d0_comb.plotOn(frame_d0_wmcut.get(), LineWidth(2), LineColor(kMagenta),
                       VLines(), Range("pidcalib"), NormRange("fitRange"));
        d0_comb.plotOn(frame_d0_wmcut.get(), LineWidth(2), LineColor(kRed),
                       VLines(), LineStyle(kDashed), Range("fitRange"));

        c_double.cd(1);
        frame_d0->Draw();
        c_double.cd(2);
        frame_d0_wmcut->Draw();

        c_double.SaveAs(opath + "/fits/d0m_wmcut_" + suffix + ".pdf");

        const double k_total     = k_comb.getVal();
        const double k_total_err = k_comb.getError();

        // Compare DaVinci and offline WM cuts

        unique_ptr<RooPlot> frame_d0_wm_comp(
            d0_m_var.frame(Title("D0 M - WM cut comparison " + tag)));

        dataset_wmcut->plotOn(frame_d0_wm_comp.get(),
                              Binning(bins_histos_d0_m));
        dataset_wmcut_offline->plotOn(frame_d0_wm_comp.get(),
                                      Binning(bins_histos_d0_m),
                                      MarkerColor(kRed));

        c_single.cd();
        frame_d0_wm_comp->Draw();
        c_single.SaveAs(opath + "/fits/d0m_wmcut_comp_" + suffix + ".pdf");

        TH1D *h_d0_wm_dv = static_cast<TH1D *>(dataset_wmcut->createHistogram(
            "d0_m_var", Binning(bins_histos_d0_m)));
        h_d0_wm_dv->SetName("h_d0_wm_dv_" + suffix);
        TH1D *h_d0_wm_offline =
            static_cast<TH1D *>(dataset_wmcut_offline->createHistogram(
                "d0_m_var", Binning(bins_histos_d0_m)));
        h_d0_wm_offline->SetName("h_d0_wm_offline_" + suffix);
        TGraphAsymmErrors tg(h_d0_wm_dv, h_d0_wm_offline, "pois");
        tg.Draw("AP");
        c_single.SaveAs(opath + "/fits/d0m_wmcut_ratio_" + suffix + ".pdf");

        // Fit sample with WM cut and mu pid
        // Rearrange plotting here so tha tmagenta retains the slope from
        // "total" fit

        unique_ptr<RooPlot> frame_d0_wmcut_mupid(
            d0_m_var.frame(Title("D0 M - With WM cut and mu PID " + tag)));
        dataset_wmcut_mupid->plotOn(frame_d0_wmcut_mupid.get(),
                                    Binning(bins_histos_d0_m));
        d0_comb.plotOn(frame_d0_wmcut_mupid.get(), LineWidth(2),
                       LineColor(kMagenta), VLines(), Range("pidcalib"),
                       NormRange("fitRange"));

        d0_comb.fitTo(*dataset_wmcut_mupid, Strategy(2), Range("fitRange"),
                      PrintLevel(0));

        d0_comb.plotOn(frame_d0_wmcut_mupid.get(), LineWidth(2),
                       LineColor(kBlue), Range(""), NormRange("fitRange"));
        d0_comb.plotOn(frame_d0_wmcut_mupid.get(), LineWidth(2),
                       LineColor(kRed), VLines(), LineStyle(kDashed),
                       Range("fitRange"));

        c_double.cd(1);
        frame_d0_wmcut->Draw();
        c_double.cd(2);
        frame_d0_wmcut_mupid->Draw();

        c_double.SaveAs(opath + "/fits/d0m_mupid_" + suffix + ".pdf");

        const double k_passed     = k_comb.getVal();
        const double k_passed_err = k_comb.getError();

        const double diff = k_passed - k_total;
        const double diff_err =
            std::sqrt(k_passed_err * k_passed_err + k_total_err * k_total_err);

        cout << "\nINFO k_passed - k_total = (" << k_passed << " +- "
             << k_passed_err << ") - (" << k_total << " +- " << k_total_err
             << ") = (" << diff << " +- " << diff_err << ")" << endl;
        cout << "INFO Pull = " << diff / diff_err << endl << endl;

        ////////
        // dm //
        ////////

        // Fit sample without WM cut
        dm_comb.fitTo(*dataset, Strategy(2), Range("fitRange"), PrintLevel(0));

        unique_ptr<RooPlot> frame_dm(
            dm_var.frame(Title("dM - No WM cut " + tag)));

        dataset->plotOn(frame_dm.get(), Binning(bins_histos_dm));
        dm_comb.plotOn(frame_dm.get(), LineWidth(2), LineColor(kBlue));

        // Fit sample with WM cut
        dm_comb.fitTo(*dataset_wmcut, Strategy(2), Range("fitRange"),
                      PrintLevel(0));

        unique_ptr<RooPlot> frame_dm_wmcut(
            dm_var.frame(Title("dM - With WM cut " + tag)));

        dataset_wmcut->plotOn(frame_dm_wmcut.get(), Binning(bins_histos_dm));
        dm_comb.plotOn(frame_dm_wmcut.get(), LineWidth(2), LineColor(kBlue));

        c_double.cd(1);
        frame_dm->Draw();
        c_double.cd(2);
        frame_dm_wmcut->Draw();

        c_double.SaveAs(opath + "/fits/dm_wmcut_" + suffix + ".pdf");

        const double c_total     = c_comb.getVal();
        const double c_total_err = c_comb.getError();

        // Fit sample with WM cut and mu pid
        // Rearrange plotting here so tha tmagenta retains the slope from
        // "total" fit

        unique_ptr<RooPlot> frame_dm_wmcut_mupid(
            dm_var.frame(Title("D0 M - With WM cut and mu PID " + tag)));
        dataset_wmcut_mupid->plotOn(frame_dm_wmcut_mupid.get(),
                                    Binning(bins_histos_dm));
        dm_comb.plotOn(frame_dm_wmcut_mupid.get(), LineWidth(2),
                       LineColor(kMagenta));

        dm_comb.fitTo(*dataset_wmcut_mupid, Strategy(2), Range("fitRange"),
                      PrintLevel(0));

        dm_comb.plotOn(frame_dm_wmcut_mupid.get(), LineWidth(2),
                       LineColor(kBlue));

        c_double.cd(1);
        frame_dm_wmcut->Draw();
        c_double.cd(2);
        frame_dm_wmcut_mupid->Draw();

        c_double.SaveAs(opath + "/fits/dm_mupid_" + suffix + ".pdf");

        const double c_passed     = c_comb.getVal();
        const double c_passed_err = c_comb.getError();

        const double diff_dm = c_passed - c_total;
        const double diff_dm_err =
            std::sqrt(c_passed_err * c_passed_err + c_total_err * c_total_err);

        cout << "\nINFO c_passed - c_total = (" << c_passed << " +- "
             << c_passed_err << ") - (" << c_total << " +- " << c_total_err
             << ") = (" << diff_dm << " +- " << diff_dm_err << ")" << endl;
        cout << "INFO Pull = " << diff_dm / diff_dm_err << endl << endl;
      }
    }
  }
}