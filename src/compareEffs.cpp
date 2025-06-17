#include <iostream>
#include <string>

#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH3.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TRatioPlot.h"
#include "TString.h"
#include "TStyle.h"

#include <cxxopts.hpp>

using std::cout, std::endl;
using std::string, std::vector;

int main(int argc, char** argv) {
  cxxopts::Options argOpts("compareEffs", "Compare efficiency histograms");

  // clang-format off
  argOpts.add_options()
    ("h,help", "print help")
    ("p,particles", "specify probed particle",
     cxxopts::value<vector<string>>()->default_value("pi,k"))
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
  const auto    particles = parsedArgs["particles"].as<vector<string>>();
  const TString opath     = parsedArgs["output"].as<string>();
  const auto    vmu       = parsedArgs["vmu"].as<bool>();
  const auto    fake_mu   = parsedArgs["fake_mu"].as<bool>();

  // // VMU
  // const TString path_ref_k =
  //     "/home/lmeyerga/misid-unfold/histos/ctrl_sample/"
  //     "rdx-25_04_08_11_31-true_to_tag_glacier-2016/kTrueToMuTag_nom.root";
  // const TString path_ref_pi =
  //     "/home/lmeyerga/misid-unfold/histos/ctrl_sample/"
  //     "rdx-25_04_08_11_31-true_to_tag_glacier-2016/piTrueToMuTag_nom.root";
  // const TString path_new_k =
  //     "/home/lmeyerga/misid-unfold/histos/"
  //     "rdx-25_06_16_09_00-misid-mc-corrections-vmu/"
  //     "kTrueToMuTag_nom.root";
  // const TString path_new_pi =
  //     "/home/lmeyerga/misid-unfold/histos/"
  //     "rdx-25_06_16_09_00-misid-mc-corrections-vmu/"
  //     "piTrueToMuTag_nom.root";

  // ISO+CTRL
  const TString path_ref_k =
      "/home/lmeyerga/misid-unfold/histos/default/"
      "rdx-25_04_08_06_47-true_to_tag_glacier-2016/kTrueToMuTag_nom.root";
  const TString path_ref_pi =
      "/home/lmeyerga/misid-unfold/histos/default/"
      "rdx-25_04_08_06_47-true_to_tag_glacier-2016/piTrueToMuTag_nom.root";
  const TString path_new_k =
      "/home/lmeyerga/misid-unfold/histos/rdx-25_06_16_08_58-misid-mc-corrections/"
      "kTrueToMuTag_nom.root";
  const TString path_new_pi =
      "/home/lmeyerga/misid-unfold/histos/rdx-25_06_16_08_58-misid-mc-corrections/"
      "piTrueToMuTag_nom.root";

  // // FAKE_MU
  // const TString path_ref_k =
  //     "/home/lmeyerga/misid-unfold/histos/default/"
  //     "rdx-25_04_08_06_47-true_to_tag_glacier-2016/kTrueToMuTag_denom.root";
  // const TString path_ref_pi =
  //     "/home/lmeyerga/misid-unfold/histos/default/"
  //     "rdx-25_04_08_06_47-true_to_tag_glacier-2016/piTrueToMuTag_denom.root";
  // const TString path_new_k =
  //     "/home/lmeyerga/misid-unfold/histos/"
  //     "rdx-25_06_16_08_59-misid-mc-corrections-fake_mu/"
  //     "kTrueToMuTag_denom.root";
  // const TString path_new_pi =
  //     "/home/lmeyerga/misid-unfold/histos/"
  //     "rdx-25_06_16_08_59-misid-mc-corrections-fake_mu/"
  //     "piTrueToMuTag_denom.root";

  TFile ifile_ref_k(path_ref_k, "READ");
  TFile ifile_ref_pi(path_ref_pi, "READ");
  TFile ifile_new_k(path_new_k, "READ");
  TFile ifile_new_pi(path_new_pi, "READ");

  TH3D* histo_ref_k;
  TH3D* histo_ref_pi;
  TH3D* histo_new_k;
  TH3D* histo_new_pi;

  ifile_ref_k.GetObject("eff", histo_ref_k);
  ifile_ref_pi.GetObject("eff", histo_ref_pi);
  ifile_new_k.GetObject("eff", histo_new_k);
  ifile_new_pi.GetObject("eff", histo_new_pi);

  // cout << "DEBUG Printing histo_new_k contents:" << endl;
  // histo_new_k->Print("all");
  // cout << endl;

  // c.cd();
  // histo_ref_k->Draw("BOX2");
  // c.SaveAs(opath + "/th3_k_ref.pdf");
  // histo_ref_pi->Draw("BOX2");
  // c.SaveAs(opath + "/th3_pi_ref.pdf");
  // histo_new_k->Draw("BOX2");
  // c.SaveAs(opath + "/th3_k_new.pdf");
  // histo_new_pi->Draw("BOX2");
  // c.SaveAs(opath + "/th3_pi_new.pdf");

  const int n_p_bins       = histo_ref_k->GetNbinsX();
  const int n_eta_bins     = histo_ref_k->GetNbinsY();
  const int n_ntracks_bins = histo_ref_k->GetNbinsZ();

  for (int ntrks_idx = 0; ntrks_idx < n_ntracks_bins; ntrks_idx++) {
    for (int eta_idx = 0; eta_idx < n_eta_bins; eta_idx++) {
      TString suffix = TString::Format("%d_%d", ntrks_idx, eta_idx);

      // I tried using histo_ref_k->ProjectionX("_px", eta_idx + 1, eta_idx + 1,
      // ntrks_idx + 1, ntrks_idx + 1) to get a 1D slice of the 3D histogram,
      // but for some reason it produces wrong results. To be sure, let's do it
      // manually
      TH1D histo_ref_k_proj("histo_ref_k_proj" + suffix, "", n_p_bins,
                            histo_ref_k->GetXaxis()->GetXbins()->GetArray());
      TH1D histo_ref_pi_proj("histo_ref_pi_proj" + suffix, "", n_p_bins,
                             histo_ref_pi->GetXaxis()->GetXbins()->GetArray());
      TH1D histo_new_k_proj("histo_new_k_proj" + suffix, "", n_p_bins,
                            histo_new_k->GetXaxis()->GetXbins()->GetArray());
      TH1D histo_new_pi_proj("histo_new_pi_proj" + suffix, "", n_p_bins,
                             histo_new_pi->GetXaxis()->GetXbins()->GetArray());

      for (int p_idx = 0; p_idx < n_p_bins; p_idx++) {
        const int bin_idx =
            histo_ref_k->GetBin(p_idx + 1, eta_idx + 1, ntrks_idx + 1);

        histo_ref_k_proj.SetBinContent(p_idx + 1,
                                       histo_ref_k->GetBinContent(bin_idx));
        histo_ref_pi_proj.SetBinContent(p_idx + 1,
                                        histo_ref_pi->GetBinContent(bin_idx));
        histo_new_k_proj.SetBinContent(p_idx + 1,
                                       histo_new_k->GetBinContent(bin_idx));
        histo_new_pi_proj.SetBinContent(p_idx + 1,
                                        histo_new_pi->GetBinContent(bin_idx));

        histo_ref_k_proj.SetBinError(p_idx + 1,
                                     histo_ref_k->GetBinError(bin_idx));
        histo_ref_pi_proj.SetBinError(p_idx + 1,
                                      histo_ref_pi->GetBinError(bin_idx));
        histo_new_k_proj.SetBinError(p_idx + 1,
                                     histo_new_k->GetBinError(bin_idx));
        histo_new_pi_proj.SetBinError(p_idx + 1,
                                      histo_new_pi->GetBinError(bin_idx));
      }

      histo_ref_k_proj.SetLineColor(kBlack);
      histo_ref_pi_proj.SetLineColor(kBlack);
      histo_new_k_proj.SetLineColor(kRed);
      histo_new_pi_proj.SetLineColor(kRed);

      // cout << "DEBUG Printing histo_new_k_slice " << suffix
      //      << " contents:" << endl;
      // histo_new_k_slice->Print("all");
      // cout << endl;

      // cout << "DEBUG Printing histo_new_k_slice FIX " << suffix
      //      << " contents:" << endl;
      // histo_new_k_proj.Print("all");
      // cout << endl;

      // cout << "DEBUG Printing histo_ref_k_slice " << suffix
      //      << " contents:" << endl;
      // histo_ref_k_slice->Print("all");
      // cout << endl;

      // cout << "DEBUG Printing histo_ref_k_slice FIX " << suffix
      //      << " contents:" << endl;
      // histo_ref_k_proj.Print("all");
      // cout << endl;

      // c.cd();
      // histo_ref_k_proj.Draw("");
      // c.SaveAs(opath + "/histo_k_ref_" + suffix + ".pdf");
      // histo_ref_pi_proj.Draw("");
      // c.SaveAs(opath + "/histo_pi_ref_" + suffix + ".pdf");
      // histo_new_k_proj.Draw("");
      // c.SaveAs(opath + "/histo_k_new_" + suffix + ".pdf");
      // histo_new_pi_proj.Draw("");
      // c.SaveAs(opath + "/histo_pi_new_" + suffix + ".pdf");

      TCanvas c("c", "c", 640, 480);

      TGraphErrors tge_ref_k(&histo_ref_k_proj);
      tge_ref_k.SetName("tge_ref_k");
      tge_ref_k.SetLineColor(kBlack);
      tge_ref_k.SetMarkerColor(kAzure);
      tge_ref_k.SetMarkerStyle(8);
      TGraphErrors tge_ref_pi(&histo_ref_pi_proj);
      tge_ref_pi.SetName("tge_ref_pi");
      tge_ref_pi.SetLineColor(kBlack);
      tge_ref_pi.SetMarkerColor(kAzure);
      tge_ref_pi.SetMarkerStyle(8);
      TGraphErrors tge_new_k(&histo_new_k_proj);
      tge_new_k.SetName("tge_new_k");
      tge_new_k.SetLineColor(kBlack);
      tge_new_k.SetMarkerColor(kRed);
      tge_new_k.SetMarkerStyle(8);
      TGraphErrors tge_new_pi(&histo_new_pi_proj);
      tge_new_pi.SetName("tge_new_pi");
      tge_new_pi.SetLineColor(kBlack);
      tge_new_pi.SetMarkerColor(kRed);
      tge_new_pi.SetMarkerStyle(8);

      tge_ref_k.SetTitle("K ref");
      tge_ref_pi.SetTitle("#pi ref");
      tge_new_k.SetTitle("K new");
      tge_new_pi.SetTitle("#pi new");

      TMultiGraph mg_k("mg_k", "");
      TMultiGraph mg_pi("mg_pi", "");

      mg_k.Add(&tge_new_k);
      mg_k.Add(&tge_ref_k);
      mg_pi.Add(&tge_new_pi);
      mg_pi.Add(&tge_ref_pi);

      const auto    axis_eta     = histo_ref_k->GetYaxis();
      const auto    axis_ntracks = histo_ref_k->GetZaxis();
      const TString tag =
          TString::Format("(%.0f < nTracks < %.0f, %.1f < #eta < %.1f)",
                          axis_ntracks->GetBinLowEdge(ntrks_idx + 1),
                          axis_ntracks->GetBinUpEdge(ntrks_idx + 1),
                          axis_eta->GetBinLowEdge(eta_idx + 1),
                          axis_eta->GetBinUpEdge(eta_idx + 1));

      mg_k.SetTitle("K " + tag);
      mg_k.Draw("AP");
      c.BuildLegend();
      c.SaveAs(opath + "/comparison_k_" + suffix + ".pdf");

      mg_pi.SetTitle("#pi " + tag);
      mg_pi.Draw("AP");
      c.BuildLegend();
      c.SaveAs(opath + "/comparison_pi_" + suffix + ".pdf");

      c.cd();
      TRatioPlot rp_k(&histo_new_k_proj, &histo_ref_k_proj);
      rp_k.Draw();
      double max_k = 1.1 * std::max(histo_new_k_proj.GetMaximum(),
                                    histo_ref_k_proj.GetMaximum());
      rp_k.GetUpperRefYaxis()->SetRangeUser(0, max_k);
      rp_k.GetLowerRefYaxis()->SetTitle("New/Ref");
      rp_k.GetUpperPad()->cd();
      TLegend legend_k(0.3, 0.7, 0.7, 0.85);
      legend_k.AddEntry(histo_new_k_proj.GetName(), "New", "l");
      legend_k.AddEntry(histo_ref_k_proj.GetName(), "Ref", "le");
      legend_k.Draw();
      c.SaveAs(opath + "/ratio_k_" + suffix + ".pdf");

      c.cd();
      TRatioPlot rp_pi(&histo_new_pi_proj, &histo_ref_pi_proj);
      rp_pi.Draw();
      double max_pi = 1.1 * std::max(histo_new_pi_proj.GetMaximum(),
                                     histo_ref_pi_proj.GetMaximum());
      rp_pi.GetUpperRefYaxis()->SetRangeUser(0, max_pi);
      rp_pi.GetLowerRefYaxis()->SetTitle("New/Ref");
      rp_pi.GetUpperPad()->cd();
      TLegend legend_pi(0.3, 0.7, 0.7, 0.85);
      legend_pi.AddEntry(histo_new_pi_proj.GetName(), "New", "l");
      legend_pi.AddEntry(histo_ref_pi_proj.GetName(), "Ref", "le");
      legend_pi.Draw();
      c.SaveAs(opath + "/ratio_pi_" + suffix + ".pdf");
    }
  }

  ifile_ref_k.Close();
  ifile_ref_pi.Close();
  ifile_new_k.Close();
  ifile_new_pi.Close();
}