// Author: Lucas Meyer Garcia
// License: BSD 2-clause
//
// Description: Calculate efficiency of PIDCalib selection requirements for its
// D*+ -> D0(-> K-pi+) pi+ sample

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <Math/Vector4D.h>
#include "TChain.h"
#include "TEfficiency.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include "RooArgSet.h"
#include "RooDataSet.h"

#include "misid.h"
#include "utils.h"

#include <yaml-cpp/yaml.h>
#include <cxxopts.hpp>

using ROOT::Math::PxPyPzMVector;

int main(int argc, char **argv) {
  cxxopts::Options argOpts("PIDCalibSelEffs",
                           "Calculate selection efficiencies for PIDCalib's "
                           "K/pi calibration sample");

  // clang-format off
  argOpts.add_options()
    ("h,help", "Print help")
    ("d,debug", "Enable debug mode",
     cxxopts::value<bool>()->default_value("false"))
    ("c,config", "Specify input YAML config file",
     cxxopts::value<string>())
    ("o,output", "Specify output folder",
     cxxopts::value<string>()->default_value("gen/"))
    ("f,file", "Specify output YAML file base name",
     cxxopts::value<string>()->default_value("pidcalib_sel_effs.yml"))
    ("p,particles", "Specify probed particle",
     cxxopts::value<vector<string>>()->default_value("k,pi"))
    ("r,run2ang", "Use Run2Ang PID cuts",
     cxxopts::value<bool>()->default_value("false"))
    ;
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  const auto particles = parsedArgs["particles"].as<vector<string>>();
  const auto ymlFile   = parsedArgs["config"].as<string>();
  const auto ymlConfig = YAML::LoadFile(ymlFile)["pidcalib_sel_effs"];
  const auto run2ang   = parsedArgs["run2ang"].as<bool>();

  if (run2ang) {
    cout << "INFO Using Run 2 Angular Analysis cuts" << endl;
  } else {
    cout << "INFO Using Run 1 RDx Analysis cuts" << endl;
  }

  const string opath = parsedArgs["output"].as<string>() + "/";
  cout << "INFO Output will be saved in " << opath << endl;

  const vector<string> years_mc = {"2016", "2017", "2018"};

  // Define histogram to easily determine kinematical bins
  TH3D histo_binning("histo_binning", ";p;#eta;nTracks", N_BINS_P, BINS_P,
                     N_BINS_ETA, BINS_ETA, N_BINS_NTRACKS, BINS_NTRACKS);

  constexpr bool fix_mass_eff = true;

  // Create structure to hold yaml data
  map<string,
      map<string,
          map<string,
              map<string, map<int, map<int, map<int, map<string, int>>>>>>>>
      ymlContent;

  const auto outFile = opath + parsedArgs["file"].as<string>();
  cout << "INFO Output will be saved in " << outFile << endl;

  for (const auto &[sample, sample_name] : SAMPLES) {
    cout << "INFO Producing " << sample_name << " selection efficiencies"
         << endl;

    for (const auto &probe : particles) {
      cout << "INFO Selecting " << probe << endl;
      const TString tag = (probe == "pi") ? "k" : "pi";

      // Histograms to store counts
      TH2D sel_counts_passed_dif_uprich("sel_counts_passed_dif_uprich",
                                        ";p;#eta", N_BINS_P, BINS_P, N_BINS_ETA,
                                        BINS_ETA);
      TH2D sel_counts_passed_dif_dwrich("sel_counts_passed_dif_dwrich",
                                        ";p;#eta", N_BINS_P, BINS_P, N_BINS_ETA,
                                        BINS_ETA);

      TH2D sel_counts_failed_dif("sel_counts_failed_dif", ";p;#eta", N_BINS_P,
                                 BINS_P, N_BINS_ETA, BINS_ETA);
      TH3D sel_counts_nondif("sel_counts_nondif", ";p;#eta;nTracks", N_BINS_P,
                             BINS_P, N_BINS_ETA, BINS_ETA, N_BINS_NTRACKS,
                             BINS_NTRACKS);

      TH2D total_counts_passed_dif_uprich("total_counts_passed_dif_uprich",
                                          ";p;#eta", N_BINS_P, BINS_P,
                                          N_BINS_ETA, BINS_ETA);
      TH2D total_counts_passed_dif_dwrich("total_counts_passed_dif_dwrich",
                                          ";p;#eta", N_BINS_P, BINS_P,
                                          N_BINS_ETA, BINS_ETA);

      TH2D total_counts_failed_dif("total_counts_failed_dif", ";p;#eta",
                                   N_BINS_P, BINS_P, N_BINS_ETA, BINS_ETA);
      TH3D total_counts_nondif("total_counts_nondif", ";p;#eta;nTracks",
                               N_BINS_P, BINS_P, N_BINS_ETA, BINS_ETA,
                               N_BINS_NTRACKS, BINS_NTRACKS);

      // Histograms to monitor corresponding mass distributions
      gStyle->SetOptStat("eou");
      array<array<array<TH2D *, N_BINS_P>, N_BINS_ETA>, N_BINS_NTRACKS>
          d0m_dm_passed_dif_uprich;

      array<array<array<TH2D *, N_BINS_P>, N_BINS_ETA>, N_BINS_NTRACKS>
          d0m_dm_passed_dif_dwrich;

      array<array<array<TH2D *, N_BINS_P>, N_BINS_ETA>, N_BINS_NTRACKS>
          d0m_dm_passed_nondif;

      array<array<array<TH2D *, N_BINS_P>, N_BINS_ETA>, N_BINS_NTRACKS>
          d0m_dm_failed_dif;

      array<array<array<TH2D *, N_BINS_P>, N_BINS_ETA>, N_BINS_NTRACKS>
          d0m_dm_failed_nondif;

      for (int ntrks_idx = 0; ntrks_idx < N_BINS_NTRACKS; ntrks_idx++) {
        for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
          for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
            const TString suffix = TString::Format("%s_%d_%d_%d", probe.c_str(),
                                                   ntrks_idx, eta_idx, p_idx);
            const TString tag    = TString::Format(
                "(%.0f < nTracks < %.0f, %.1f < #eta < %.1f, %.0f < p < %.0f)",
                BINS_NTRACKS[ntrks_idx], BINS_NTRACKS[ntrks_idx + 1],
                BINS_ETA[eta_idx], BINS_ETA[eta_idx + 1], BINS_P[p_idx],
                BINS_P[p_idx + 1]);

            d0m_dm_passed_dif_uprich[ntrks_idx][eta_idx][p_idx] = new TH2D(
                "d0m_dm_passed_dif_uprich_" + suffix,
                "Passed DiF - before RICH1 " + tag + ";" + d0_m_x_dm_label, 150,
                0.5, 3.5, 100, 0.1, 0.4);
            d0m_dm_passed_dif_dwrich[ntrks_idx][eta_idx][p_idx] = new TH2D(
                "d0m_dm_passed_dif_dwrich_" + suffix,
                "Passed DiF - after RICH1 " + tag + ";" + d0_m_x_dm_label, 150,
                0.5, 3.5, 100, 0.1, 0.4);
            d0m_dm_passed_nondif[ntrks_idx][eta_idx][p_idx] =
                new TH2D("d0m_dm_passed_nondif_" + suffix,
                         "Passed non-DiF " + tag + ";" + d0_m_x_dm_label, 150,
                         0.5, 3.5, 100, 0.1, 0.4);

            d0m_dm_failed_dif[ntrks_idx][eta_idx][p_idx] =
                new TH2D("d0m_dm_failed_dif_" + suffix,
                         "Failed DiF " + tag + ";" + d0_m_x_dm_label, 150, 0.5,
                         3.5, 100, 0.1, 0.4);
            d0m_dm_failed_nondif[ntrks_idx][eta_idx][p_idx] =
                new TH2D("d0m_dm_failed_nondif_" + suffix,
                         "Failed non-DiF " + tag + ";" + d0_m_x_dm_label, 150,
                         0.5, 3.5, 100, 0.1, 0.4);
          }
        }
      }

      for (auto year : years_mc) {
        cout << "INFO Reading " << year << " ntuples" << endl;

        // Open MC files
        const auto mc_path = ymlConfig["ntps"][year][probe].as<string>();
        TChain     ch_mc("tree");
        ch_mc.Add(mc_path.c_str());

        cout << "INFO Opened MC files:" << endl;
        print_files(ch_mc);

        // Define variable to access input ntuples
        int dst_trueid, d0_trueid, probe_trueid, probe_daughter0_trueid,
            probe_daughter1_trueid, probe_mc_mom_nd, ntracks, dst_vtx_ndof,
            d0_vtx_ndof;

        double dst_m, d0_m, probe_p, probe_pz, probe_pt, tag_p, probe_dllmu,
            probe_dlle, tag_pz, tag_pt, k_track_chi2ndof, pi_track_chi2ndof,
            spi_track_chi2ndof, k_px, k_py, pi_px, pi_py, dst_vtx_chi2,
            d0_vtx_chi2, k_ghostprob, pi_ghostprob, spi_ghostprob,
            probe_true_origz, probe_daughter0_true_origz, spi_p, spi_pt;

        float probe_mu_ubdt;

        bool probe_ismuon, probe_hasmuon;

        cout << "INFO Setting input branches " << endl;
        ch_mc.SetBranchStatus("*", false);
        ch_mc.SetBranchAddress("dst_M", &dst_m);
        ch_mc.SetBranchAddress("dst_TRUEID", &dst_trueid);
        ch_mc.SetBranchAddress("dst_ENDVERTEX_CHI2", &dst_vtx_chi2);
        ch_mc.SetBranchAddress("dst_ENDVERTEX_NDOF", &dst_vtx_ndof);
        ch_mc.SetBranchAddress("d0_M", &d0_m);
        ch_mc.SetBranchAddress("d0_TRUEID", &d0_trueid);
        ch_mc.SetBranchAddress("d0_ENDVERTEX_CHI2", &d0_vtx_chi2);
        ch_mc.SetBranchAddress("d0_ENDVERTEX_NDOF", &d0_vtx_ndof);
        // ch_mc.SetBranchAddress("d0_DIRA_OWNPV", &d0_dira);
        // ch_mc.SetBranchAddress("spi_TRUEID", &spi_trueid);
        // ch_mc.SetBranchAddress("spi_MC_MOTHER_ID", &spi_mom_trueid);
        // ch_mc.SetBranchAddress("spi_TRUEENDVERTEX_Z", &spi_trueendvtx_z);
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
        ch_mc.SetBranchAddress((probe + "_TRUEORIGINVERTEX_Z").c_str(),
                               &probe_true_origz);
        ch_mc.SetBranchAddress((probe + "_DAUGHTER0_ID").c_str(),
                               &probe_daughter0_trueid);
        ch_mc.SetBranchAddress((probe + "_DAUGHTER0_ORIGINVERTEX_Z").c_str(),
                               &probe_daughter0_true_origz);
        ch_mc.SetBranchAddress((probe + "_DAUGHTER1_ID").c_str(),
                               &probe_daughter1_trueid);
        ch_mc.SetBranchAddress((probe + "_MC_MOTHER_ND").c_str(),
                               &probe_mc_mom_nd);
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

        const int entries_mc = ch_mc.GetEntries();
        cout << "INFO Starting MC event loop over " << entries_mc << " entries"
             << endl;

        int eta_bin = 0, p_bin = 0, ntrks_bin = 0;

        for (int evt = 0; evt < entries_mc; evt++) {
          ch_mc.GetEntry(evt);

          if (std::abs(dst_trueid) != Dst_ID) continue;
          if (std::abs(d0_trueid) != D0_ID) continue;

          // Conditional cuts
          if (!probe_hasmuon) continue;

          if ((k_ghostprob > 0.4) || (pi_ghostprob > 0.4) ||
              (spi_ghostprob > 0.4))
            continue;

          if ((k_track_chi2ndof > 3.) || (pi_track_chi2ndof > 3.) ||
              (spi_track_chi2ndof > 3.))
            continue;

          // Add a loose cut on D0 vertex as conditional since in RDx
          // the misidentified mu track must still satisfy vertex quality
          // requirements with the D0. This matches the HLT2 level cut (which is
          // tightened at stripping level to 6).
          if ((d0_vtx_chi2 / d0_vtx_ndof) > 15.) continue;

          // Range cuts
          if (probe_p < 3000. || probe_p >= 100000. || ntracks >= 600) continue;

          const double probe_eta =
              0.5 * std::log((probe_p + probe_pz) / (probe_p - probe_pz));
          if (probe_eta < 1.7 || probe_eta >= 5.0) continue;

          // Determine kinematical bin
          const int kin_bin =
              histo_binning.FindBin(probe_p, probe_eta, ntracks);
          histo_binning.GetBinXYZ(kin_bin, p_bin, eta_bin, ntrks_bin);

          // Check for decay in flight of probe hadron
          const bool dif =
              check_dif(probe_trueid, probe_daughter0_trueid,
                        probe_daughter1_trueid, probe_mc_mom_nd, probe);

          // IMPORTANT: The FAKE_MU pid requirement is !isMuon, but here we flip
          // it and calculate the complementary efficiency so that the K/pi
          // misid case (with less stats) falls in the "passed" sub-sample. The
          // proper efficiency is calculated afterwards as 1 - flipped_eff.
          const bool pid_ok =
              check_mu_pid(probe_ismuon, probe_dllmu, probe_dlle, probe_mu_ubdt,
                           sample, run2ang);

          const double dm = (dst_m - d0_m) * 0.001;
          d0_m            = d0_m * 0.001;

          // Define wrong mass hypotheses
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

          const double wm_pipi = (p_k_wm + p_pi).M();
          const double wm_pik  = (p_k_wm + p_pi_wm).M();
          const double wm_kk   = (p_k + p_pi_wm).M();

          const double dif_z = std::abs(probe_trueid) == MU_ID
                                   ? probe_true_origz
                                   : probe_daughter0_true_origz;

          // Fill total counts
          if (fix_mass_eff) {
            if (pid_ok) {
              if (dif && (before_t3(dif_z) || !run2ang)) {
                if (!after_rich1(dif_z)) {
                  total_counts_passed_dif_uprich.Fill(probe_p, probe_eta);
                  d0m_dm_passed_dif_uprich[ntrks_bin - 1][eta_bin - 1]
                                          [p_bin - 1]
                                              ->Fill(d0_m, dm);
                } else {
                  total_counts_passed_dif_dwrich.Fill(probe_p, probe_eta);
                  d0m_dm_passed_dif_dwrich[ntrks_bin - 1][eta_bin - 1]
                                          [p_bin - 1]
                                              ->Fill(d0_m, dm);
                }
              } else {
                total_counts_nondif.Fill(probe_p, probe_eta, ntracks);
                d0m_dm_passed_nondif[ntrks_bin - 1][eta_bin - 1][p_bin - 1]
                    ->Fill(d0_m, dm);
              }
            } else {
              if (dif) {
                total_counts_failed_dif.Fill(probe_p, probe_eta);
                d0m_dm_failed_dif[ntrks_bin - 1][eta_bin - 1][p_bin - 1]->Fill(
                    d0_m, dm);
              } else {
                total_counts_nondif.Fill(probe_p, probe_eta, ntracks);
                d0m_dm_failed_nondif[ntrks_bin - 1][eta_bin - 1][p_bin - 1]
                    ->Fill(d0_m, dm);
              }
            }
          } else {
            // In previous code, all the cuts below are applied to the
            // denominator
            if ((tag_p < 2000.) || (tag_pt < 250.)) continue;
            if ((probe_p < 2000.) || (probe_pt < 250.)) continue;
            if ((spi_p < 1000.) || (spi_pt < 100.)) continue;

            if ((dst_vtx_chi2 / dst_vtx_ndof) > 15 ||
                (d0_vtx_chi2 / d0_vtx_ndof) > 10)
              continue;

            if ((std::max(probe_pt, tag_pt) < 1000.) ||
                ((p_k + p_pi).Pt() < 1500.))
              continue;

            // Veto wrong-mass hypotheses
            // pi pi
            if (std::abs(wm_pipi - D0_M * 1000.) < 25.) continue;
            // pi K (wrong sign)
            if (std::abs(wm_pik - D0_M * 1000.) < 25.) continue;
            // K K
            if (std::abs(wm_kk - D0_M * 1000.) < 25.) continue;
            if (pid_ok) {
              if (dif && (before_t3(dif_z) || !run2ang)) {
                if (!after_rich1(dif_z)) {
                  total_counts_passed_dif_uprich.Fill(probe_p, probe_eta);
                  d0m_dm_passed_dif_uprich[ntrks_bin - 1][eta_bin - 1]
                                          [p_bin - 1]
                                              ->Fill(d0_m, dm);
                } else {
                  total_counts_passed_dif_dwrich.Fill(probe_p, probe_eta);
                  d0m_dm_passed_dif_dwrich[ntrks_bin - 1][eta_bin - 1]
                                          [p_bin - 1]
                                              ->Fill(d0_m, dm);
                }
              } else {
                total_counts_nondif.Fill(probe_p, probe_eta, ntracks);
                d0m_dm_passed_nondif[ntrks_bin - 1][eta_bin - 1][p_bin - 1]
                    ->Fill(d0_m, dm);
              }
            } else {
              if (dif) {
                total_counts_failed_dif.Fill(probe_p, probe_eta);
                d0m_dm_failed_dif[ntrks_bin - 1][eta_bin - 1][p_bin - 1]->Fill(
                    d0_m, dm);
              } else {
                total_counts_nondif.Fill(probe_p, probe_eta, ntracks);
                d0m_dm_failed_nondif[ntrks_bin - 1][eta_bin - 1][p_bin - 1]
                    ->Fill(d0_m, dm);
              }
            }
          }

          // Apply relevant PIDCalib cuts

          const bool in_fit_window =
              in_range(D0_M_min, d0_m,
                       get_d0m_upper_limit(probe, p_bin - 1, eta_bin - 1)) &&
              in_range(DM_min, dm, DM_max);

          if (!in_fit_window) continue;

          // Fill passed counts
          if (fix_mass_eff) {
            // if ((tag_p < 2000.) || (tag_pt < 250.)) continue;
            if (probe_pt < 250.) continue;
            // if ((spi_p < 1000.) || (spi_pt < 100.)) continue;

            if ((dst_vtx_chi2 / dst_vtx_ndof) > 15 ||
                (d0_vtx_chi2 / d0_vtx_ndof) > 10)
              continue;

            if ((std::max(probe_pt, tag_pt) < 1000.) ||
                ((p_k + p_pi).Pt() < 1500.))
              continue;

            // if (d0_dira < 0.9999) continue;

            // Veto wrong-mass hypotheses
            // pi pi
            if (std::abs(wm_pipi - D0_M * 1000.) < 25.) continue;
            // pi K (wrong sign)
            if (std::abs(wm_pik - D0_M * 1000.) < 25.) continue;
            // K K
            if (std::abs(wm_kk - D0_M * 1000.) < 25.) continue;

            if (pid_ok) {
              if (dif && (before_t3(dif_z) || !run2ang)) {
                if (!after_rich1(dif_z)) {
                  sel_counts_passed_dif_uprich.Fill(probe_p, probe_eta);
                } else {
                  sel_counts_passed_dif_dwrich.Fill(probe_p, probe_eta);
                }
              } else {
                sel_counts_nondif.Fill(probe_p, probe_eta, ntracks);
              }
            } else {
              if (dif) {
                sel_counts_failed_dif.Fill(probe_p, probe_eta);
              } else {
                sel_counts_nondif.Fill(probe_p, probe_eta, ntracks);
              }
            }
          } else {
            if (pid_ok) {
              if (dif && (before_t3(dif_z) || !run2ang)) {
                if (!after_rich1(dif_z)) {
                  sel_counts_passed_dif_uprich.Fill(probe_p, probe_eta);
                } else {
                  sel_counts_passed_dif_dwrich.Fill(probe_p, probe_eta);
                }
              } else {
                sel_counts_nondif.Fill(probe_p, probe_eta, ntracks);
              }
            } else {
              if (dif) {
                sel_counts_failed_dif.Fill(probe_p, probe_eta);
              } else {
                sel_counts_nondif.Fill(probe_p, probe_eta, ntracks);
              }
            }
          }
        }
      }

      // Plot mass histograms
      TCanvas c_six("c_six", "c_six", 1440, 720);
      c_six.Divide(3, 2);

      for (int ntrks_idx = 0; ntrks_idx < N_BINS_NTRACKS; ntrks_idx++) {
        for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
          for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
            c_six.cd(1);
            gPad->SetLogz();
            d0m_dm_passed_dif_uprich[ntrks_idx][eta_idx][p_idx]->SetMinimum(1.);
            d0m_dm_passed_dif_uprich[ntrks_idx][eta_idx][p_idx]->Draw("COLZ");

            c_six.cd(2);
            gPad->SetLogz();
            d0m_dm_passed_dif_dwrich[ntrks_idx][eta_idx][p_idx]->SetMinimum(1.);
            d0m_dm_passed_dif_dwrich[ntrks_idx][eta_idx][p_idx]->Draw("COLZ");

            c_six.cd(3);
            gPad->SetLogz();
            d0m_dm_passed_nondif[ntrks_idx][eta_idx][p_idx]->SetMinimum(1.);
            d0m_dm_passed_nondif[ntrks_idx][eta_idx][p_idx]->Draw("COLZ");

            c_six.cd(5);
            gPad->SetLogz();
            d0m_dm_failed_dif[ntrks_idx][eta_idx][p_idx]->SetMinimum(1.);
            d0m_dm_failed_dif[ntrks_idx][eta_idx][p_idx]->Draw("COLZ");

            c_six.cd(6);
            gPad->SetLogz();
            d0m_dm_failed_nondif[ntrks_idx][eta_idx][p_idx]->SetMinimum(1.);
            d0m_dm_failed_nondif[ntrks_idx][eta_idx][p_idx]->Draw("COLZ");

            const TString suffix =
                TString::Format("%s_%s_%d_%d_%d", probe.c_str(),
                                sample_name.c_str(), ntrks_idx, eta_idx, p_idx);
            c_six.SaveAs(opath + "/h2_d0m_dm_" + suffix + ".pdf");

            delete d0m_dm_passed_dif_uprich[ntrks_idx][eta_idx][p_idx];
            delete d0m_dm_passed_dif_dwrich[ntrks_idx][eta_idx][p_idx];
            delete d0m_dm_passed_nondif[ntrks_idx][eta_idx][p_idx];
            delete d0m_dm_failed_dif[ntrks_idx][eta_idx][p_idx];
            delete d0m_dm_failed_nondif[ntrks_idx][eta_idx][p_idx];
          }
        }
      }

      // Move counts from histos to structure
      for (int ntrks_idx = 0; ntrks_idx < N_BINS_NTRACKS; ntrks_idx++) {
        for (int eta_idx = 0; eta_idx < N_BINS_ETA; eta_idx++) {
          for (int p_idx = 0; p_idx < N_BINS_P; p_idx++) {
            const double p_bin     = p_idx + 1;
            const double eta_bin   = eta_idx + 1;
            const double ntrks_bin = ntrks_idx + 1;
            cout << "\nINFO Bin " << ntrks_idx << " " << eta_idx << " " << p_idx
                 << endl;

            // Passed DiF before RICH
            const int sel_count_passed_dif_uprich =
                sel_counts_passed_dif_uprich.GetBinContent(p_bin, eta_bin);
            const int total_count_passed_dif_uprich =
                total_counts_passed_dif_uprich.GetBinContent(p_bin, eta_bin);

            ymlContent[sample_name][probe]["passed"]["dif_uprich"][p_idx]
                      [eta_idx][ntrks_idx]["sel"] = sel_count_passed_dif_uprich;
            ymlContent[sample_name][probe]["passed"]["dif_uprich"][p_idx]
                      [eta_idx][ntrks_idx]["total"] =
                          total_count_passed_dif_uprich;

            if (total_count_passed_dif_uprich > 0) {
              const double eff_passed_dif_uprich =
                  sel_count_passed_dif_uprich * 1.0 /
                  total_count_passed_dif_uprich;
              const double eff_passed_dif_uprich_unc_hi =
                  TEfficiency::Wilson(total_count_passed_dif_uprich,
                                      sel_count_passed_dif_uprich, ONE_SIGMA,
                                      true) -
                  eff_passed_dif_uprich;
              const double eff_passed_dif_uprich_unc_lo =
                  eff_passed_dif_uprich -
                  TEfficiency::Wilson(total_count_passed_dif_uprich,
                                      sel_count_passed_dif_uprich, ONE_SIGMA,
                                      false);

              cout << "INFO - Passed dif up:  " << eff_passed_dif_uprich * 100.
                   << " + " << eff_passed_dif_uprich_unc_hi * 100. << " - "
                   << eff_passed_dif_uprich_unc_lo * 100. << "% ("
                   << sel_count_passed_dif_uprich << " / "
                   << total_count_passed_dif_uprich << ")" << endl;
            } else {
              cout << "INFO - Passed dif up:  INVALID ("
                   << sel_count_passed_dif_uprich << " / "
                   << total_count_passed_dif_uprich << ")" << endl;
            }

            // Passed DiF after RICH
            const int sel_count_passed_dif_dwrich =
                sel_counts_passed_dif_dwrich.GetBinContent(p_bin, eta_bin);
            const int total_count_passed_dif_dwrich =
                total_counts_passed_dif_dwrich.GetBinContent(p_bin, eta_bin);

            ymlContent[sample_name][probe]["passed"]["dif_dwrich"][p_idx]
                      [eta_idx][ntrks_idx]["sel"] = sel_count_passed_dif_dwrich;
            ymlContent[sample_name][probe]["passed"]["dif_dwrich"][p_idx]
                      [eta_idx][ntrks_idx]["total"] =
                          total_count_passed_dif_dwrich;

            if (total_count_passed_dif_dwrich > 0) {
              const double eff_passed_dif_dwrich =
                  sel_count_passed_dif_dwrich * 1.0 /
                  total_count_passed_dif_dwrich;
              const double eff_passed_dif_dwrich_unc_hi =
                  TEfficiency::Wilson(total_count_passed_dif_dwrich,
                                      sel_count_passed_dif_dwrich, ONE_SIGMA,
                                      true) -
                  eff_passed_dif_dwrich;
              const double eff_passed_dif_dwrich_unc_lo =
                  eff_passed_dif_dwrich -
                  TEfficiency::Wilson(total_count_passed_dif_dwrich,
                                      sel_count_passed_dif_dwrich, ONE_SIGMA,
                                      false);

              cout << "INFO - Passed dif dw:  " << eff_passed_dif_dwrich * 100.
                   << " + " << eff_passed_dif_dwrich_unc_hi * 100. << " - "
                   << eff_passed_dif_dwrich_unc_lo * 100. << "% ("
                   << sel_count_passed_dif_dwrich << " / "
                   << total_count_passed_dif_dwrich << ")" << endl;
            } else {
              cout << "INFO - Passed dif dw:  INVALID ("
                   << sel_count_passed_dif_dwrich << " / "
                   << total_count_passed_dif_dwrich << ")" << endl;
            }

            // Passed non-DiF
            const int sel_count_passed_nondif =
                sel_counts_nondif.GetBinContent(p_bin, eta_bin, ntrks_bin);
            const int total_count_passed_nondif =
                total_counts_nondif.GetBinContent(p_bin, eta_bin, ntrks_bin);

            ymlContent[sample_name][probe]["passed"]["nondif"][p_idx][eta_idx]
                      [ntrks_idx]["sel"] = sel_count_passed_nondif;
            ymlContent[sample_name][probe]["passed"]["nondif"][p_idx][eta_idx]
                      [ntrks_idx]["total"] = total_count_passed_nondif;

            if (total_count_passed_nondif > 0) {
              const double eff_passed_nondif =
                  sel_count_passed_nondif * 1.0 / total_count_passed_nondif;
              const double eff_passed_nondif_unc_hi =
                  TEfficiency::Wilson(total_count_passed_nondif,
                                      sel_count_passed_nondif, ONE_SIGMA,
                                      true) -
                  eff_passed_nondif;
              const double eff_passed_nondif_unc_lo =
                  eff_passed_nondif -
                  TEfficiency::Wilson(total_count_passed_nondif,
                                      sel_count_passed_nondif, ONE_SIGMA,
                                      false);

              cout << "INFO - Passed non-dif: " << eff_passed_nondif * 100.
                   << " + " << eff_passed_nondif_unc_hi * 100. << " - "
                   << eff_passed_nondif_unc_lo * 100. << "% ("
                   << sel_count_passed_nondif << " / "
                   << total_count_passed_nondif << ")" << endl;
            } else {
              cout << "INFO - Passed non-dif: INVALID ("
                   << sel_count_passed_nondif << " / "
                   << total_count_passed_nondif << ")" << endl;
            }

            // Failed DiF
            const int sel_count_failed_dif =
                sel_counts_failed_dif.GetBinContent(p_bin, eta_bin);
            const int total_count_failed_dif =
                total_counts_failed_dif.GetBinContent(p_bin, eta_bin);

            ymlContent[sample_name][probe]["failed"]["dif"][p_idx][eta_idx]
                      [ntrks_idx]["sel"] = sel_count_failed_dif;
            ymlContent[sample_name][probe]["failed"]["dif"][p_idx][eta_idx]
                      [ntrks_idx]["total"] = total_count_failed_dif;

            if (total_count_failed_dif > 0) {
              const double eff_failed_dif =
                  sel_count_failed_dif * 1.0 / total_count_failed_dif;
              const double eff_failed_dif_unc_hi =
                  TEfficiency::Wilson(total_count_failed_dif,
                                      sel_count_failed_dif, ONE_SIGMA, true) -
                  eff_failed_dif;
              const double eff_failed_dif_unc_lo =
                  eff_failed_dif - TEfficiency::Wilson(total_count_failed_dif,
                                                       sel_count_failed_dif,
                                                       ONE_SIGMA, false);

              cout << "INFO - Failed dif:     " << eff_failed_dif * 100.
                   << " + " << eff_failed_dif_unc_hi * 100. << " - "
                   << eff_failed_dif_unc_lo * 100. << "% ("
                   << sel_count_failed_dif << " / " << total_count_failed_dif
                   << ")" << endl;
            } else {
              cout << "INFO - Failed dif:     INVALID (" << sel_count_failed_dif
                   << " / " << total_count_failed_dif << ")" << endl;
            }

            // Failed non-DiF
            const int sel_count_failed_nondif =
                sel_counts_nondif.GetBinContent(p_bin, eta_bin, ntrks_bin);
            const int total_count_failed_nondif =
                total_counts_nondif.GetBinContent(p_bin, eta_bin, ntrks_bin);

            ymlContent[sample_name][probe]["failed"]["nondif"][p_idx][eta_idx]
                      [ntrks_idx]["sel"] = sel_count_failed_nondif;

            ymlContent[sample_name][probe]["failed"]["nondif"][p_idx][eta_idx]
                      [ntrks_idx]["total"] = total_count_failed_nondif;

            if (total_count_failed_nondif > 0) {
              const double eff_failed_nondif =
                  sel_count_failed_nondif * 1.0 / total_count_failed_nondif;
              const double eff_failed_nondif_unc_hi =
                  TEfficiency::Wilson(total_count_failed_nondif,
                                      sel_count_failed_nondif, ONE_SIGMA,
                                      true) -
                  eff_failed_nondif;
              const double eff_failed_nondif_unc_lo =
                  eff_failed_nondif -
                  TEfficiency::Wilson(total_count_failed_nondif,
                                      sel_count_failed_nondif, ONE_SIGMA,
                                      false);

              cout << "INFO - Failed non-dif: " << eff_failed_nondif * 100.
                   << " + " << eff_failed_nondif_unc_hi * 100. << " - "
                   << eff_failed_nondif_unc_lo * 100. << "% ("
                   << sel_count_failed_nondif << " / "
                   << total_count_failed_nondif << ")" << endl;
            } else {
              cout << "INFO - Failed non-dif: INVALID ("
                   << sel_count_failed_nondif << " / "
                   << total_count_failed_nondif << ")" << endl;
            }
          }
        }
      }
    }
  }
  // Write yaml file for this year
  YAML::Emitter ymlOut;
  ymlOut << ymlContent;

  std::ofstream ofout(outFile);
  ofout << ymlOut.c_str();
  ofout.close();
}