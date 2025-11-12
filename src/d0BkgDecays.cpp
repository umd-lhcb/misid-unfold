#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TEfficiency.h"
#include "TH1D.h"
#include "TString.h"

#include <yaml-cpp/yaml.h>
#include <cxxopts.hpp>

#include "utils.h"

using std::cout, std::endl, std::right, std::left, std::setw, std::to_string;
using std::string, std::unique_ptr, std::unordered_map, std::vector, std::map;

///////////////
// Constants //
///////////////

constexpr int Dst_ID = 413;

constexpr double ONE_SIGMA = 0.682689492137086;

constexpr double DM_min   = 0.141;  // GeV
constexpr double DM_max   = 0.153;  // GeV
constexpr double D0_M_min = 1.825;  // GeV
constexpr double D0_M_max = 1.910;  // GeV

const std::unordered_map<int, std::string> id_to_string{
    {1000070140, "N14"},
    {1000040080, "Be8"},
    {1000030070, "Li7"},
    {1000020040, "alpha"},
    {-1000020040, "anti-alpha"},
    {1000010030, "tritium"},
    {-1000010030, "anti-tritium"},
    {1000010020, "deuteron"},
    {-1000010020, "anti-deuteron"},
    {20213, "a_1+"},
    {-20213, "a_1-"},
    {10323, "K_1+"},
    {-10323, "K_1-"},
    {10313, "K_10"},
    {-10313, "anti-K_10"},
    {3112, "Sigma-"},
    {-3112, "anti-Sigma+"},
    {2212, "p+"},
    {-2212, "anti-p-"},
    {2112, "n0"},
    {-2112, "anti-n0"},
    {413, "D*+"},
    {-413, "D*-"},
    {421, "D0"},
    {-421, "anti-D0"},
    {333, "phi"},
    {-333, "phi"},
    {331, "eta'"},
    {-331, "eta'"},
    {325, "K_2*+"},
    {-325, "K_2*-"},
    {323, "K*+"},
    {-323, "K*-"},
    {321, "K+"},
    {-321, "K-"},
    {313, "K*0"},
    {-313, "anti-K*0"},
    {311, "K0"},
    {-311, "anti-K0"},
    {310, "K_S0"},
    {-310, "K_S0"},
    {223, "omega"},
    {-223, "omega"},
    {221, "eta"},
    {-221, "eta"},
    {213, "rho+"},
    {-213, "rho-"},
    {211, "pi+"},
    {-211, "pi-"},
    {130, "K_L0"},
    {-130, "K_L0"},
    {113, "rho0"},
    {-113, "rho0"},
    {111, "pi0"},
    {-111, "pi0"},
    {22, "gamma"},
    {-22, "gamma"},
    {14, "nu_mu"},
    {-14, "anti-nu_mu"},
    {13, "mu-"},
    {-13, "mu+"},
    {12, "nu_e"},
    {-12, "anti-nu_e"},
    {11, "e-"},
    {-11, "e+"}};

constexpr double BINS_P[] = {3e+3, 6e+3, 10e+3, 15.6e+3, 27e+3, 60e+3, 100e+3};
constexpr int    N_BINS_P = sizeof(BINS_P) / sizeof(double) - 1;

// Comparator function
bool comp_pdg_ids(const int &a, const int &b) {
  if (abs(a) != abs(b)) {
    return abs(a) > abs(b);
  } else {
    return a > b;
  }
}

int main(int argc, char **argv) {
  cxxopts::Options argOpts(
      "GetMisIDCorrections",
      "Calculate misid corrections (CrystalBall core -> core + tail).");

  // clang-format off
  argOpts.add_options()
    ("h,help", "Print help")
    ("d,debug", "Enable debug mode",
     cxxopts::value<bool>()->default_value("false"))
    ("c,config", "Specify input YAML config file",
     cxxopts::value<string>())
    ("o,output", "Specify output folder",
     cxxopts::value<string>()->default_value("gen/"))
    ("f,file", "Specify output YAML config file",
     cxxopts::value<string>())
    ;
  // clang-format on

  auto parsedArgs = argOpts.parse(argc, argv);
  if (parsedArgs.count("help")) {
    cout << argOpts.help() << endl;
    return 0;
  }

  const auto ymlFile   = parsedArgs["config"].as<string>();
  const auto ymlConfig = YAML::LoadFile(ymlFile)["misid_corrections"];

  const string opath =
      parsedArgs["output"].as<string>() + "/" + parsedArgs["file"].as<string>();
  cout << "INFO Output will be saved in " << opath << endl;

  int dst_trueid, d0_trueid, d0_dauther0_id, d0_dauther1_id, d0_dauther2_id,
    d0_dauther3_id, d0_dauther4_id, ntracks;
  double dst_m, d0_m, probe_p, probe_pz, k_track_chi2ndof, pi_track_chi2ndof,
      spi_track_chi2ndof, probe_mudll, probe_edll;
  float probe_mu_ubdt, nshared;
  bool  probe_ismuon, probe_hasmuon;

  TH1D histo_binning("histo_binning", ";p", N_BINS_P, BINS_P);

  const vector<string> samples = {"iso_ctrl", "vmu", "fake_mu"};

  const vector<string> particles = {"k", "pi"};

  map<string, map<string, map<string, map<int, map<string, int>>>>> ymlContent;

  for (auto &sample : samples) {
    cout << "INFO Producing " << sample << " constraints" << endl;

    unordered_map<string, vector<map<string, int>>> counts_particles_passed;
    unordered_map<string, vector<map<string, int>>> counts_particles_failed;

    for (auto &probe : particles) {
      counts_particles_passed.emplace(probe, vector<map<string, int>>(6));
      counts_particles_failed.emplace(probe, vector<map<string, int>>(6));

      auto &counts_passed = counts_particles_passed[probe];
      auto &counts_failed = counts_particles_failed[probe];

      // Open and loop over MC files
      const auto mc_path = ymlConfig["mc_ntps"]["mb"][probe].as<string>();
      TChain     ch_mc("tree");

      cout << "INFO Opening MC files: " << mc_path << endl;
      ch_mc.Add(mc_path.c_str());

      cout << "INFO Opened MC files:" << endl;
      print_files(ch_mc);

      ch_mc.SetBranchStatus("*", false);
      ch_mc.SetBranchAddress("dst_TRUEID", &dst_trueid);
      ch_mc.SetBranchAddress("dst_M", &dst_m);
      ch_mc.SetBranchAddress("d0_TRUEID", &d0_trueid);
      ch_mc.SetBranchAddress("d0_M", &d0_m);
      ch_mc.SetBranchAddress("d0_DAUGHTER0_ID", &d0_dauther0_id);
      ch_mc.SetBranchAddress("d0_DAUGHTER1_ID", &d0_dauther1_id);
      ch_mc.SetBranchAddress("d0_DAUGHTER2_ID", &d0_dauther2_id);
      ch_mc.SetBranchAddress("d0_DAUGHTER3_ID", &d0_dauther3_id);
      ch_mc.SetBranchAddress("d0_DAUGHTER4_ID", &d0_dauther4_id);
      ch_mc.SetBranchAddress((probe + "_P").c_str(), &probe_p);
      ch_mc.SetBranchAddress((probe + "_PZ").c_str(), &probe_pz);
      ch_mc.SetBranchAddress((probe + "_isMuon").c_str(), &probe_ismuon);
      ch_mc.SetBranchAddress((probe + "_hasMuon").c_str(), &probe_hasmuon);
      ch_mc.SetBranchAddress((probe + "_PIDmu").c_str(), &probe_mudll);
      ch_mc.SetBranchAddress((probe + "_PIDe").c_str(), &probe_edll);
      ch_mc.SetBranchAddress("k_TRACK_CHI2NDOF", &k_track_chi2ndof);
      ch_mc.SetBranchAddress("pi_TRACK_CHI2NDOF", &pi_track_chi2ndof);
      ch_mc.SetBranchAddress("spi_TRACK_CHI2NDOF", &spi_track_chi2ndof);
      ch_mc.SetBranchAddress((probe + "_bdt_mu").c_str(), &probe_mu_ubdt);
      ch_mc.SetBranchAddress("nTracks", &ntracks);
      ch_mc.SetBranchAddress("MuonNShared", &nshared);
      
      const int entries_mc = ch_mc.GetEntries();
      cout << "INFO Starting MC event loop over " << entries_mc << " entries"
           << endl;
      for (int evt = 0; evt < entries_mc; evt++) {
        ch_mc.GetEntry(evt);

        if (abs(dst_trueid) != Dst_ID) {
          continue;
        }

        // Conditional cuts
        if (!probe_hasmuon) continue;

        if ((probe_p < 3000.) || (probe_p > 100000.) || (ntracks >= 600)) {
          continue;
        }

        const double probe_eta =
            0.5 * log((probe_p + probe_pz) / (probe_p - probe_pz));
        if (probe_eta < 1.7 || probe_eta >= 5.0) continue;

        if ((k_track_chi2ndof > 3.) || (pi_track_chi2ndof > 3.) ||
            (spi_track_chi2ndof > 3.))
          continue;

        // Determine kinematical bin
        const int p_bin = histo_binning.FindBin(probe_p);

        const bool reduced_fit_range = (probe == "pi") && (p_bin <= 2);

        const double dm = (dst_m - d0_m) * 0.001;
        d0_m            = d0_m * 0.001;

        const bool in_fit_window = reduced_fit_range
                                       ? in_range(D0_M_min, d0_m, 1.900) &&
                                             in_range(DM_min, dm, DM_max)
                                       : in_range(D0_M_min, d0_m, D0_M_max) &&
                                             in_range(DM_min, dm, DM_max);
        if (!in_fit_window) continue;

        bool pid_ok;
        if (sample == "fake_mu") {
          // For FAKE_MU, we calculate the complementary efficiency so that
          // the "passed" sample always corresponds to the K/pi misid case
          pid_ok = probe_ismuon;
        } else if (sample == "vmu") {
          pid_ok = probe_ismuon && probe_mudll > 2.0 && probe_edll < -1.0 &&
                   probe_mu_ubdt < 0.65 && nshared==1;
          // pid_ok = probe_ismuon && probe_mudll > 2.0 && probe_edll < 1.0 &&
          //          probe_mu_ubdt < 0.25;
        } else {
          pid_ok = probe_ismuon && probe_mudll > 2.0 && probe_edll < -1.0 &&
                   probe_mu_ubdt > 0.65 && nshared==1;
        }

        const int sign = (d0_trueid >= 0) ? 1 : -1;

        const int p_idx = p_bin - 1;

        vector<int> daughter_ids{d0_dauther0_id * sign, d0_dauther1_id * sign,
                                 d0_dauther2_id * sign, d0_dauther3_id * sign,
                                 d0_dauther4_id * sign};
        sort(daughter_ids.begin(), daughter_ids.end(), comp_pdg_ids);
        string decay_string = (id_to_string.count(d0_trueid * sign) > 0)
                                  ? id_to_string.at(d0_trueid * sign) + " ->"
                                  : to_string(d0_trueid * sign) + " ->";
        for (const auto &p_id : daughter_ids) {
          if (p_id != 0) {
            const string part_string = (id_to_string.count(p_id) > 0)
                                           ? id_to_string.at(p_id)
                                           : to_string(p_id);
            decay_string += " " + part_string;
          }
        }

        if (pid_ok) {
          counts_passed[p_idx][decay_string]++;
        } else {
          counts_failed[p_idx][decay_string]++;
        }
      }
    }

    // Probably not the most efficient solution, but avoids duplicating the loop
    // below
    const unordered_map<string, unordered_map<string, vector<map<string, int>>>>
        counts_particles_all{{"passed", counts_particles_passed},
                             {"failed", counts_particles_failed}};

    for (auto &counts_particles_subsample : counts_particles_all) {
      const auto &subsample        = counts_particles_subsample.first;
      const auto &counts_particles = counts_particles_subsample.second;
      for (auto &particle_counts : counts_particles) {
        const auto &particle = particle_counts.first;
        const auto &counts   = particle_counts.second;
        cout << "\nINFO Printing " << subsample << " counts for " << particle
             << endl;
        for (unsigned p_bin = 0; p_bin < counts.size(); p_bin++) {
          cout << " - p bin " << p_bin << endl;
          const auto &count_decs = counts[p_bin];
          int         bin_total  = 0;
          for (const auto &dec_count : count_decs) {
            bin_total += dec_count.second;
          }
          for (const auto &dec_count : count_decs) {
            const auto &decay = dec_count.first;
            const auto &count = dec_count.second;
            cout << " -- " << left << setw(30) << decay << ": " << right
                 << setw(5) << count << " / " << left << setw(5) << bin_total
                 << " = " << count * 100. / bin_total << "%" << endl;
          }
          const int bin_signal = count_decs.at("D0 -> K- pi+");
          const int bin_bkg    = bin_total - bin_signal;

          // Bkg / Sig
          cout << " --- bkg / sig: " << bin_bkg << " / " << bin_signal << " = "
               << bin_bkg * 100. / bin_signal << "%" << endl;

          // Bkg / Total
          const double f_bkg = (bin_bkg * 1.0) / bin_total;
          const double f_bkg_unc_hi =
              TEfficiency::Wilson(bin_total, bin_bkg, ONE_SIGMA, true) - f_bkg;
          const double f_bkg_unc_lo =
              f_bkg - TEfficiency::Wilson(bin_total, bin_bkg, ONE_SIGMA, false);
          cout << " --- bkg / total: " << bin_bkg << " / " << bin_total
               << " = (" << f_bkg * 100. << " + " << f_bkg_unc_hi * 100.
               << " - " << f_bkg_unc_lo * 100. << ")%" << endl;

          // Sig / Total
          const double f_sig = (bin_signal * 1.0) / bin_total;
          const double f_sig_unc_hi =
              TEfficiency::Wilson(bin_total, bin_signal, ONE_SIGMA, true) -
              f_sig;
          const double f_sig_unc_lo =
              f_sig -
              TEfficiency::Wilson(bin_total, bin_signal, ONE_SIGMA, false);
          cout << " --- sig / total: " << bin_signal << " / " << bin_total
               << " = (" << f_sig * 100. << " + " << f_sig_unc_hi * 100.
               << " - " << f_sig_unc_lo * 100. << ")% (" << bin_signal << ", "
               << bin_total << ")\n"
               << endl;

          ymlContent[sample][particle][subsample][p_bin]["signal"] = bin_signal;
          ymlContent[sample][particle][subsample][p_bin]["total"]  = bin_total;
        }
      }
    }
  }

  YAML::Emitter ymlOut;
  ymlOut << ymlContent;

  std::ofstream ofout(opath);
  ofout << ymlOut.c_str();
  ofout.close();

  return 0;
}
