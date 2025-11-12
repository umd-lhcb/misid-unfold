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
#include "TRatioPlot.h"
#include "TString.h"

#include <cxxopts.hpp>

using std::cout, std::endl;
using std::string, std::vector;

int main(int argc, char** argv) {
  cxxopts::Options argOpts("compareEffs", "Compare efficiency histograms");

  // clang-format off
  argOpts.add_options()
    ("h,help", "print help")
    ("o,output", "specify output folder",
     cxxopts::value<string>())
    ("r,refPath", "specify path to old efficiencies",
     cxxopts::value<string>())
    ("n,newPath", "specify path to new efficiencies",
     cxxopts::value<string>())
    ;
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  // Read arguments
  const TString opath   = parsedArgs["output"].as<string>();
  const TString pathRef = parsedArgs["refPath"].as<string>();
  const TString pathNew = parsedArgs["newPath"].as<string>();

  // Get efficiencies to be compared
  TFile ifile_ref(pathRef, "READ");
  TFile ifile_new(pathNew, "READ");

  TH1* histo_ref = nullptr;
  TH1* histo_new = nullptr;

  ifile_ref.GetObject("eff", histo_ref);

  if (!histo_ref) {
    cout << "FATAL Could not find histogram 'eff' in " << pathRef << endl;
    exit(1);
  }

  ifile_new.GetObject("eff", histo_new);

  if (!histo_new) {
    cout << "FATAL Could not find histogram 'eff' in " << pathNew << endl;
    exit(1);
  }

  const int n_p_bins       = histo_ref->GetNbinsX();
  const int n_eta_bins     = histo_ref->GetNbinsY();
  const int n_ntracks_bins = histo_ref->GetNbinsZ();

  TCanvas c("c", "c", 640, 480);
  TCanvas c_ratio("c_ratio", "c_ratio", 640, 640);

  for (int ntrks_idx = 0; ntrks_idx < n_ntracks_bins; ntrks_idx++) {
    for (int eta_idx = 0; eta_idx < n_eta_bins; eta_idx++) {
      c.cd();

      TString suffix = TString::Format("%d_%d", ntrks_idx, eta_idx);

      // I tried using histo_ref_k->ProjectionX("_px", eta_idx + 1, eta_idx + 1,
      // ntrks_idx + 1, ntrks_idx + 1) to get a 1D slice of the 3D histogram,
      // but for some reason it produces wrong results. To be sure, let's do it
      // manually
      TH1D histo_ref_proj("histo_ref_proj" + suffix, "", n_p_bins,
                          histo_ref->GetXaxis()->GetXbins()->GetArray());
      TH1D histo_new_proj("histo_new_proj" + suffix, "", n_p_bins,
                          histo_new->GetXaxis()->GetXbins()->GetArray());

      for (int p_idx = 0; p_idx < n_p_bins; p_idx++) {
        const int bin_idx =
            histo_ref->GetBin(p_idx + 1, eta_idx + 1, ntrks_idx + 1);

        histo_ref_proj.SetBinContent(p_idx + 1,
                                     histo_ref->GetBinContent(bin_idx));
        histo_new_proj.SetBinContent(p_idx + 1,
                                     histo_new->GetBinContent(bin_idx));

        histo_ref_proj.SetBinError(p_idx + 1, histo_ref->GetBinError(bin_idx));
        histo_new_proj.SetBinError(p_idx + 1, histo_new->GetBinError(bin_idx));
      }

      histo_ref_proj.SetLineColor(kBlack);
      histo_new_proj.SetLineColor(kRed);

      TGraphErrors tge_ref(&histo_ref_proj);
      tge_ref.SetName("tge_ref");
      tge_ref.SetLineColor(kBlack);
      tge_ref.SetMarkerColor(kBlack);
      tge_ref.SetMarkerStyle(8);

      TGraphErrors tge_new(&histo_new_proj);
      tge_new.SetName("tge_new");
      tge_new.SetLineColor(kRed);
      tge_new.SetMarkerColor(kRed);
      tge_new.SetMarkerStyle(8);

      tge_ref.SetTitle("ref");
      tge_new.SetTitle("new");

      TMultiGraph mg("mg", "");

      mg.Add(&tge_ref);
      mg.Add(&tge_new);

      const auto    axis_eta     = histo_ref->GetYaxis();
      const auto    axis_ntracks = histo_ref->GetZaxis();
      const TString tag =
          TString::Format("(%.0f < nTracks < %.0f, %.1f < #eta < %.1f)",
                          axis_ntracks->GetBinLowEdge(ntrks_idx + 1),
                          axis_ntracks->GetBinUpEdge(ntrks_idx + 1),
                          axis_eta->GetBinLowEdge(eta_idx + 1),
                          axis_eta->GetBinUpEdge(eta_idx + 1));

      mg.SetTitle(tag);
      mg.Draw("AP");
      c.BuildLegend();
      c.SaveAs(opath + "/comparison_" + suffix + ".pdf");

      c_ratio.cd();

      TRatioPlot rp(&histo_new_proj, &histo_ref_proj);
      rp.Draw();
      const double max = 1.1 * std::max(histo_new_proj.GetMaximum(),
                                        histo_ref_proj.GetMaximum());
      rp.GetUpperRefYaxis()->SetRangeUser(0, max);
      rp.GetLowerRefYaxis()->SetTitle("New/Ref");
      rp.GetLowerRefGraph()->SetMarkerStyle(kFullDotLarge);
      rp.GetLowerRefGraph()->SetMarkerColor(kBlue);
      rp.GetLowerRefGraph()->SetLineColor(kBlue);
      rp.GetUpperPad()->cd();
      TLegend legend(0.3, 0.7, 0.7, 0.85);
      legend.AddEntry(histo_new_proj.GetName(), "New", "l");
      legend.AddEntry(histo_ref_proj.GetName(), "Ref", "le");
      legend.Draw();

      c_ratio.SaveAs(opath + "/ratio_" + suffix + ".pdf");
    }
  }

  ifile_ref.Close();
  ifile_new.Close();
}
