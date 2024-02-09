//#if !defined(__CINT__) || defined(__MAKECINT__)
#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <TROOT.h>
#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TEventList.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TRandom3.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <memory>
#include <map>
#include <tuple>
#include <array>

using namespace std;

// Include files

const bool skipToys = true;
const Int_t numToys = (skipToys) ? 1 : 100;
const Double_t nts = numToys * numToys;

/** @class MisidWeightor MisidWeightor.h BcProcessor/MisidWeightor.h
 *  
 *
 *  @author Jack Wimberley
 *  @date   2015-08-11
 */

const TString year = "16";
const TString pidfile_name = "/home/ejiang/rjpsi_misid/myPid_up.root";
const bool verbose = false;
const bool very_verbose = false;

const TString CutString = "Jpsi_L0DiMuonDecision_TOS==1 && (Jpsi_Hlt1DiMuonHighMassDecision_TOS==1 || Jpsi_Hlt1DiMuonLowMassDecision_TOS==1) && Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS==1 && Bc_ISOLATION_BDT<0.2 && Bc_MM<=6400 && Bc_MM>=3203 && Bc_ENDVERTEX_CHI2<25 && Bc_DOCA<0.15 && BachMu_IPCHI2_OWNPV>4.8 && MuM_IPCHI2_OWNPV>4 && MuP_IPCHI2_OWNPV>4 && Jpsi_MM<=3150 && Jpsi_MM>=3040 && BachMu_P>3000 && BachMu_P<100000 && nTracks<600 && BachMu_PT>750 && BachMu_isMuon==0";

const TString PiCutString = "BachMu_ProbNNpi>0.1&&BachMu_PIDK<1.0&&BachMu_PIDmu<2.0&&BachMu_PIDp<4.0&&!(BachMu_ProbNNe>0.1)&&!(BachMu_ProbNNp>0.1)&&!(BachMu_PIDmu>2)&&!(BachMu_ProbNNk>0.1)";
const TString KCutString = "BachMu_ProbNNk>0.1&&!(BachMu_ProbNNe>0.1)&&!(BachMu_PIDmu>2)";
const TString PCutString = "BachMu_ProbNNp>0.1&&BachMu_PIDK<1&&!(BachMu_ProbNNk>0.1)&&!(BachMu_ProbNNe>0.1)&&!(BachMu_PIDmu>2)";
const TString MuCutString = "BachMu_PIDmu>2&&!(BachMu_ProbNNe>0.1)";
const TString ECutString = "BachMu_ProbNNe>0.1";
const Int_t numClasses = 6;
const Int_t num = numClasses-1;

const TString particle_name[numClasses] = {"PI","K","P","MU","E","U"};

TH3F* loadRootObject(TFile& file, const TString& name) {
  std::cout << "LOADING OBJECT " << name << " FROM FILE " << file.GetName() << ":\n";
  TH3F* p = static_cast<TH3F*>(file.Get(name));
  p->SetDirectory(0);
  return p;
}

using CacheType = std::tuple<Double_t,Double_t, std::array<Double_t,numClasses> >;
using CacheEntry = std::pair<Int_t, CacheType >;

class MisIDFunction {
public: 

  /// Standard constructor
  MisIDFunction(int cuttype = 0)
    : cache(new std::map<Int_t,CacheType >()),
      pidfile(pidfile_name)
  {

    HistoMatrix[0][0] = loadRootObject(pidfile, "H2H_Pi_TO_Pi");
    HistoMatrix[1][0] = loadRootObject(pidfile, "H2H_Pi_TO_K");
    HistoMatrix[2][0] = loadRootObject(pidfile, "H2H_Pi_TO_P");
    HistoMatrix[3][0] = loadRootObject(pidfile, "H2H_Pi_TO_Mu");
    HistoMatrix[4][0] = loadRootObject(pidfile, "H2H_Pi_TO_e");
    HistoMatrix[5][0] = loadRootObject(pidfile, "H2H_Pi_TO_Ghost");
    
    HistoMatrix[0][1] = loadRootObject(pidfile, "H2H_K_TO_Pi");
    HistoMatrix[1][1] = loadRootObject(pidfile, "H2H_K_TO_K");
    HistoMatrix[2][1] = loadRootObject(pidfile, "H2H_K_TO_P");
    HistoMatrix[3][1] = loadRootObject(pidfile, "H2H_K_TO_Mu");
    HistoMatrix[4][1] = loadRootObject(pidfile, "H2H_K_TO_e");
    HistoMatrix[5][1] = loadRootObject(pidfile, "H2H_K_TO_Ghost");
    
    HistoMatrix[0][2] = loadRootObject(pidfile, "H2H_P_TO_Pi");
    HistoMatrix[1][2] = loadRootObject(pidfile, "H2H_P_TO_K");
    HistoMatrix[2][2] = loadRootObject(pidfile, "H2H_P_TO_P");
    HistoMatrix[3][2] = loadRootObject(pidfile, "H2H_P_TO_Mu");
    HistoMatrix[4][2] = loadRootObject(pidfile, "H2H_P_TO_e");
    HistoMatrix[5][2] = loadRootObject(pidfile, "H2H_P_TO_Ghost");
    
    HistoMatrix[0][3] = loadRootObject(pidfile, "H2H_Mu_TO_Pi");
    HistoMatrix[1][3] = loadRootObject(pidfile, "H2H_Mu_TO_K");
    HistoMatrix[2][3] = loadRootObject(pidfile, "H2H_Mu_TO_P");
    HistoMatrix[3][3] = loadRootObject(pidfile, "H2H_Mu_TO_Mu");
    HistoMatrix[4][3] = loadRootObject(pidfile, "H2H_Mu_TO_e");
    HistoMatrix[5][3] = loadRootObject(pidfile, "H2H_Mu_TO_Ghost");
    
    HistoMatrix[0][4] = loadRootObject(pidfile, "H2H_e_TO_Pi");
    HistoMatrix[1][4] = loadRootObject(pidfile, "H2H_e_TO_K");
    HistoMatrix[2][4] = loadRootObject(pidfile, "H2H_e_TO_P");
    HistoMatrix[3][4] = loadRootObject(pidfile, "H2H_e_TO_Mu");
    HistoMatrix[4][4] = loadRootObject(pidfile, "H2H_e_TO_e");
    HistoMatrix[5][4] = loadRootObject(pidfile, "H2H_e_TO_Ghost");

    HistoMatrix[0][5] = loadRootObject(pidfile, "H2H_Ghost_TO_Pi");
    HistoMatrix[1][5] = loadRootObject(pidfile, "H2H_Ghost_TO_K");
    HistoMatrix[2][5] = loadRootObject(pidfile, "H2H_Ghost_TO_P");
    HistoMatrix[3][5] = loadRootObject(pidfile, "H2H_Ghost_TO_Mu");
    HistoMatrix[4][5] = loadRootObject(pidfile, "H2H_Ghost_TO_e");
    HistoMatrix[5][5] = loadRootObject(pidfile, "H2H_Ghost_TO_Ghost");
    
    if (cuttype==0) { // nn
      HistoVector[0] = loadRootObject(pidfile,"EFFNN0_Pi_TO_MU");
      HistoVector[1] = loadRootObject(pidfile,"EFFNN0_K_TO_MU");
      HistoVector[2] = loadRootObject(pidfile,"EFFNN0_P_TO_MU");
      HistoVector[3] = loadRootObject(pidfile,"EFFNN0_Mu_TO_MU");
      HistoVector[4] = loadRootObject(pidfile,"EFFNN0_e_TO_MU");
      HistoVector[5] = loadRootObject(pidfile,"EFFNN0_Ghost_TO_MU");
    } else if (cuttype==1) { //dll
      HistoVector[0] = loadRootObject(pidfile,"EFFDLL0_Pi_TO_MU");
      HistoVector[1] = loadRootObject(pidfile,"EFFDLL0_K_TO_MU");
      HistoVector[2] = loadRootObject(pidfile,"EFFDLL0_P_TO_MU");
      HistoVector[3] = loadRootObject(pidfile,"EFFDLL0_Mu_TO_MU");
      HistoVector[4] = loadRootObject(pidfile,"EFFDLL0_e_TO_MU");
      HistoVector[5] = loadRootObject(pidfile,"EFFDLL0_Ghost_TO_MU");
    } else { // dllk
      HistoVector[0] = loadRootObject(pidfile,"EFFDLLK_Pi_TO_MU");
      HistoVector[1] = loadRootObject(pidfile,"EFFDLLK_K_TO_MU");
      HistoVector[2] = loadRootObject(pidfile,"EFFDLLK_P_TO_MU");
      HistoVector[3] = loadRootObject(pidfile,"EFFDLLK_Mu_TO_MU");
      HistoVector[4] = loadRootObject(pidfile,"EFFDLLK_e_TO_MU");
      HistoVector[5] = loadRootObject(pidfile,"EFFDLLK_Ghost_TO_MU");
    }
    
  };

  /// Destructor
  ~MisIDFunction() {
    pidfile.Close();
    delete cache;
  };

  Double_t operator() (Double_t *x, Double_t *ps = 0) {

    // RETRIEVE VAL
    Int_t PidBin = x[0];
    Int_t HadronHypo = ps[0];
    Bool_t returnError = ps[1];
    Bool_t useNikhefPriors = ps[2];
    Bool_t returnHadronType = ps[3];
    Int_t HadronReturnType = ps[4];
    // std::cout << ((returnError) ? "ERROR" : "VALUE") << " / ";
    // std::cout << ((useNikhefPriors) ? "HYBRID" : "UMD") << std::endl;
    // std::cout << HadronHypo << " - " << useNikhefPriors << " - " << PidBin << std::endl;
    Int_t marker = HadronHypo + 10*useNikhefPriors + 100*PidBin;
    auto it = cache->find(marker);
    if (it == cache->end()) {
    
      // CALCULATE VALUE
      std::cout << "Calculating value for bin " << PidBin << " ... " << std::endl;
      Double_t avgVal = 0.0;
      Double_t avgValSq = 0.0;
      Double_t avgValErr = 0.0;

      Double_t val_umd[numClasses] = {0};
      Double_t valsq_umd[numClasses] = {0};
      Double_t valerr_umd[numClasses] = {0};
      Double_t valbreak_umd[numClasses][numClasses] = {0}; // first index is hat category; second is weight for unfolded type
      Double_t val_hyb[numClasses] = {0};
      Double_t valsq_hyb[numClasses] = {0};
      Double_t valerr_hyb[numClasses] = {0};
      Double_t valbreak_hyb[numClasses][numClasses] = {0}; // first index is hat category; second is weight for unfolded type

      if (!HistoMatrix[0][0]->IsBinUnderflow(PidBin) && !HistoMatrix[0][0]->IsBinOverflow(PidBin)) {
	std::cout<<"not over or underflow"<<std::endl;
	TH1D temp("temp","temp",1,0,2);
	TString cut;

	Double_t piNum = 0, kNum = 0, pNum = 0, muNum = 0, eNum = 0, uNum = 0;
	Species sp = (*MisIDFunction::countCache)[PidBin];
	if (num >= 1)
	  piNum = sp.piNum;
	if (num >= 2)
	  kNum = sp.kNum;
	if (num >= 3)
	  pNum = sp.pNum;
	if (num >= 4)
	  muNum = sp.muNum;
	if (num >= 5)
	  eNum = sp.eNum;
	uNum = sp.uNum;

	std::cout << "PI NUM = " << piNum << std::endl;
	std::cout << "K NUM  = " << kNum << std::endl;
	std::cout << "P NUM  = " << pNum << std::endl;
	std::cout << "MU NUM = " << muNum << std::endl;
	std::cout << "E NUM  = " << eNum << std::endl;
	std::cout << "UN NUM = " << uNum << std::endl;

    
	Double_t totNum = piNum + kNum + pNum + muNum + eNum + uNum;

	if (totNum > 0) {

	  Double_t measured_vector_total[numClasses] = {piNum,kNum,pNum,muNum,eNum,uNum};
	  Double_t measured_vector[numClasses];
	  for (int i = 0; i < numClasses; i++) {
	    measured_vector[i] = measured_vector_total[i];
	  }
      
	  // TOY STUDY LOOP

	  TRandom3 r(6451);
      
	  static int counter = 0;
	  TH1D *measured_hist, *response_truth;
	  TH2D *response_hist;
	  measured_hist = new TH1D(TString::Format("measured_hist_%d",counter),TString::Format("measured_hist_%d",counter),numClasses,0,numClasses);
	  response_truth = new TH1D(TString::Format("response_truth_%d",counter),TString::Format("response_truth_%d",counter),numClasses,0,numClasses);
	  response_hist = new TH2D(TString::Format("response_hist_%d",counter),TString::Format("response_hist_%d",counter),numClasses,0,numClasses,numClasses,0,numClasses);
	  ++counter;
      
	  for (int toy = 0; toy < numToys; toy++) {

	    if (verbose) std::cout << "\n\n TOY #" << toy << ":\n\n";
	    Double_t response_matrix[numClasses][numClasses];
	    for (int j = 0; j < numClasses; j++) {
	      for (int i = 0; i < numClasses; i++) {
		double val = HistoMatrix[i][j]->GetBinContent(PidBin);
		double err = HistoMatrix[i][j]->GetBinError(PidBin);
		double sample = -1;
		if (val < 0) {
		  sample = 0.0;
		} else {
		  while (sample < 0) {
		    if (verbose or very_verbose or skipToys)
		      sample = val;
		    else
		      sample = r.Gaus(val,err);
		  }
		}
		response_matrix[i][j] = sample;
	      }
	    }
#ifdef GHOST_ALTERNATE
	    //static_assert(false);
	    for (int j = 0; j < num; j++) {
	      response_matrix[j][num] = 0.0;
	    }
#endif

	    if (verbose) std::cout << "H -> H RESPONSE MATRIX:" << std::endl;
	    for (int row = 0; row < numClasses; row++) {
	      measured_hist->SetBinContent(row+1,measured_vector[row]);
	      response_truth->SetBinContent(row+1,1);
	      for (int col = 0; col < numClasses; col++) {
		response_hist->SetBinContent(row+1,col+1,response_matrix[row][col]);
		if (verbose) std::cout << std::setw(10) << std::setprecision(4) << response_hist->GetBinContent(row+1,col+1) << " | ";
	      }
	      if (verbose) std::cout << std::endl;
	    }
      
	    RooUnfoldResponse response(0,response_truth,response_hist);
	    int numIters = 4;
	    RooUnfoldBayes unfold(&response,measured_hist,numIters);
	    TH1D* reco = static_cast<TH1D*>(unfold.Hreco());
	    Double_t estNum = 0.0;
	    Double_t response_vector[numClasses];
	    for (int i = 0; i < numClasses; i++) {
	      double val = reco->GetBinContent(i+1);
	      double err = reco->GetBinError(i+1);
	      estNum += val;
	      response_vector[i] = val;
	      if (toy == 0) {
		sum_measured_vector[i] += measured_vector[i];
		sum_response_vector[i] += response_vector[i];
		sum_response_vector_var[i] += err*err;
	      }
	      if (verbose) std::cout << particle_name[i] << ": ";
	      if (verbose) std::cout << measured_vector[i] << " -> " << response_vector[i] << " +/-" << reco->GetBinError(i+1) << std::endl;
	    }

	    // NOW DO CASE BY CASE BASIS
	    Double_t prior_probs_umd[numClasses];
	    Double_t prior_probs_hyb[numClasses];
	    if (verbose) std::cout << "\n A PRIORI PROBS:\n";
	    for (int i = 0; i < numClasses; i++) {
	      prior_probs_umd[i] = 1.0/(numClasses);
	      prior_probs_hyb[i] = response_vector[i] / estNum;
	      if (verbose) std::cout << particle_name[i] << ": " << prior_probs_umd[i] << " OR " << prior_probs_hyb[i] << std::endl;
	    }
	    Double_t posterior_probs_umd[numClasses];
	    Double_t posterior_probs_hyb[numClasses];
	    if (verbose) std::cout << "\n A POSTERIORI PROBS:\n";
	    for (int i = 0; i < numClasses; i++) {
	      posterior_probs_umd[i] = 0.0;
	      posterior_probs_hyb[i] = 0.0;
	      for (int j = 0; j < numClasses; j++) {
		posterior_probs_umd[i] += response_matrix[i][j]*prior_probs_umd[j];
		posterior_probs_hyb[i] += response_matrix[i][j]*prior_probs_hyb[j];
	      }
	      if (verbose) std::cout << particle_name[i] << ": " << posterior_probs_umd[i] << " OR " << posterior_probs_hyb[i] << std::endl;
	    }

	    Double_t inverted_probs_umd[numClasses][numClasses];
	    Double_t inverted_probs_hyb[numClasses][numClasses];
	    for (int i = 0; i < numClasses; i++) {
	      for (int j = 0; j < numClasses; j++) {
		inverted_probs_umd[i][j] = response_matrix[j][i]*prior_probs_umd[i]/posterior_probs_umd[j];
		inverted_probs_hyb[i][j] = response_matrix[j][i]*prior_probs_hyb[i]/posterior_probs_hyb[j];
	      }
	    }
	    if (verbose) {
	      std::cout << "\n\n INVERTED PROBS --UMD :\n";
	      for (int i = 0; i < numClasses; i++) {
		for (int j = 0; j < numClasses; j++) {
		  std::cout << std::setw(8) << inverted_probs_umd[i][j] << " | ";
		}
		std::cout << std::endl;
	      }
	      std::cout << "\n\n INVERTED PROBS --HYBRID :\n";
	      for (int i = 0; i < numClasses; i++) {
		for (int j = 0; j < numClasses; j++) {
		  std::cout << std::setw(8) << inverted_probs_hyb[i][j] << " | ";
		}
		std::cout << std::endl;
	      }
	    }
        
        
	    if (very_verbose) std::cout << "\nFAKE RATES:\n";
	    for (int toy2 = 0; toy2 < numToys; toy2++) {
	      if (very_verbose) std::cout << "\n\n NESTED TOY #" << toy2 << ":\n\n";
	      Double_t eff_vector[numClasses];
	      for (int i = 0; i < numClasses; i++) {

		if (i == 3) { // muons
		  eff_vector[i] = 0.0;
		} else {

		  double val, err;
#ifdef GHOST_ALTERNATE
		  // grab pion
		  Int_t z = (i == num) ? 0 : i;
		  val = HistoVector[z]->GetBinContent(PidBin);
		  err = HistoVector[z]->GetBinError(PidBin);
#else
		  val = HistoVector[i]->GetBinContent(PidBin);
		  err = HistoVector[i]->GetBinError(PidBin);
		  ///val = 0.90;
		  ///err = 0.01;
#endif
		  double sample = -1;
		  if (val < 0) {
		    sample = 0.0;
		  } else {
		    while (sample < 0) {
		      if (verbose or very_verbose or skipToys)
			sample = val;
		      else
			sample = r.Gaus(val,err);
		    }
		  }
		  eff_vector[i] = sample;
		}
		if (very_verbose) std::cout << particle_name[i] << " FAKE RATE : ";
		if (very_verbose) std::cout << eff_vector[i] << std::endl;
	      }

	      Double_t avg = 0.0;
	      for (int j = 0; j < numClasses; j++) {
		avg += response_vector[j]*eff_vector[j] / totNum;
	      }

	      if (very_verbose) std::cout << std::endl;

	      Double_t valMatrix_umd[numClasses][numClasses] = {0};
	      Double_t valMatrix_hyb[numClasses][numClasses] = {0};
	      Double_t valVector_umd[numClasses] = {0};
	      Double_t valVector_hyb[numClasses] = {0};
	      for (int j = 0; j < numClasses; j++) {
		for (int i = 0; i < numClasses; i++) {
		  valMatrix_umd[j][i] = inverted_probs_umd[i][j] * eff_vector[i];
		  valMatrix_hyb[j][i] = inverted_probs_hyb[i][j] * eff_vector[i];
		  valVector_umd[j] += valMatrix_umd[j][i];
		  valVector_hyb[j] += valMatrix_hyb[j][i];
		}
		if (very_verbose)
		  std::cout << particle_name[j] << ": " << valVector_umd[j] << " OR " << valVector_hyb[j] << std::endl;
	      }

	      avgVal += avg;
	      avgValSq += avg*avg;

	      for (int i = 0; i < numClasses; ++i) {
		val_umd[i] += valVector_umd[i];
		val_hyb[i] += valVector_hyb[i];
		valsq_umd[i] += valVector_umd[i]*valVector_umd[i];
		valsq_hyb[i] += valVector_hyb[i]*valVector_hyb[i];
		for (int j = 0; j < numClasses; ++j) {
		  valbreak_umd[i][j] += valMatrix_umd[i][j];
		  valbreak_hyb[i][j] += valMatrix_hyb[i][j];
		}
	      }
	    }

	    if (very_verbose) std::cin.get();
	    if (very_verbose) exit(1);
	  }

	  delete measured_hist;
	  delete response_truth;
	  delete response_hist;

	  avgVal /= numToys*numToys;
	  avgValSq /= numToys*numToys;
	  avgValErr = sqrt(avgValSq-avgVal*avgVal);

	  for (int i = 0; i < numClasses; ++i) {
	    val_umd[i] /= nts;
	    valsq_umd[i] /= nts;
	    valerr_umd[i] = sqrt(valsq_umd[i]-val_umd[i]*val_umd[i]);
	    val_hyb[i] /= nts;
	    valsq_hyb[i] /= nts;
	    valerr_hyb[i] = sqrt(valsq_hyb[i]-val_hyb[i]*val_hyb[i]);
	    for (int j = 0; j < numClasses; ++j) {
	      valbreak_umd[i][j] /= nts;
	      valbreak_hyb[i][j] /= nts;
	    }
	  }

	  std::cout << std::endl;
	  std::cout << "AVG: " << avgVal << " +/- " << avgValErr << std::endl;
	  for (int i = 0; i < numClasses; ++i) {
	    std::cout << particle_name[i] << ":  " << val_umd[i]  << " +/- " << valerr_umd[i]
		      << " OR " << val_hyb[i] << " +/- " << valerr_hyb[i] << std::endl;
	  }
	  if (val_umd[numClasses-1] == 1.0) {
	    std::cout << "WARNING! " << std::endl;
	    std::cin.get();
	  }
	}
      }

      // STORE VAL
      for (int i = 0; i < numClasses; ++i) {
	std::array<Double_t,numClasses> break_umd {};
	std::array<Double_t,numClasses> break_hyb {};
	for (int j = 0; j < numClasses; ++j) {
	  break_umd[j] = valbreak_umd[i][j];
	  break_hyb[j] = valbreak_hyb[i][j];
	}
	cache->insert(CacheEntry(i+100*PidBin,std::make_tuple(val_umd[i],valerr_umd[i],break_umd)));
	cache->insert(CacheEntry(i+10+100*PidBin,std::make_tuple(val_hyb[i],valerr_hyb[i],break_hyb)));
      }
      std::array<Double_t,numClasses> aempty {};
      cache->insert(CacheEntry(numClasses+100*PidBin,std::make_tuple(avgVal,avgValErr,aempty))); // nikhef
      cache->insert(CacheEntry(numClasses+10+100*PidBin,std::make_tuple(avgVal,avgValErr,aempty))); // nikhef
    }

    // Now its been inserted
    if (verbose) std::cin.get();
    if (verbose) exit(1);
    it = cache->find(marker);
    Double_t val, err;
    std::array<Double_t,numClasses> breakdown {};
    std::tie(val,err,breakdown) = it->second;

    if (returnError)
      return err;
    else if (returnHadronType)
      return breakdown[HadronReturnType];
    else
      return val;

  };
  
  struct Species {
    Double_t piNum;
    Double_t kNum;
    Double_t pNum;
    Double_t muNum;
    Double_t eNum;
    Double_t uNum;
  };
  
  void createCountCache(TTree* tree);
  static std::map<Int_t,Species>* countCache;
  static TEventList listPi;
  static TEventList listK;
  static TEventList listP;
  static TEventList listMu;
  static TEventList listE;
  static TEventList listNone;
  
  void printAggregateResults() {
    double sum_measured = 0;
    double sum_response = 0;
    std::cout << "TOTAL PARTICLE NUMBERS, RAW AND UNFOLDED:\n";
    for (int i = 0; i < numClasses; i++) {
      std::cout << particle_name[i] << ": ";
      std::cout << sum_measured_vector[i] << " -> "
		<< sum_response_vector[i] << " +/-" << sqrt(sum_response_vector_var[i]) << std::endl;
      sum_measured += sum_measured_vector[i];
      sum_response += sum_response_vector[i];
    }
    std::cout << "TOTAL PARTICLE NUMBERS, RAW AND UNFOLDED:\n";
    for (int i = 0; i < numClasses; i++) {
      std::cout << particle_name[i] << ": ";
      std::cout << sum_measured_vector[i]/sum_measured << " -> "
		<< sum_response_vector[i]/sum_response << " +/-" << sqrt(sum_response_vector_var[i])/sum_response << std::endl;
    }
  }

private:
  
  std::map<Int_t,CacheType >* cache;
  TFile pidfile;
  bool useNikhefPriors;
  TH3F* HistoMatrix[numClasses][numClasses];
  TH3F* HistoVector[numClasses];

  double sum_measured_vector[numClasses] = {0,0,0,0,0,0};
  double sum_response_vector[numClasses] = {0,0,0,0,0,0};
  double sum_response_vector_var[numClasses] = {0,0,0,0,0,0};

};

std::map<Int_t,MisIDFunction::Species>* MisIDFunction::countCache = nullptr;
TEventList MisIDFunction::listPi("listPi","listPi");
TEventList MisIDFunction::listK("listK","listK");
TEventList MisIDFunction::listP("listP","listP");
TEventList MisIDFunction::listMu("listMu","listMu");
TEventList MisIDFunction::listE("listE","listE");
TEventList MisIDFunction::listNone("listNone","listNone");

void MisIDFunction::createCountCache(TTree* tree) {

  MisIDFunction::countCache = new std::map<Int_t,Species>();
  
  // Readers would be nice for uniformity, but could introduce errors
  // TTreeReader reader(tree);
  // TTreeReaderValue<Int_t> myPidBin(reader,"PidBin");
  // TTreeReaderValue<Int_t> myHadronHypo(reader,"HadronHypo"); // built in HadronHypo is out of date; this will be worth fixing
  // TTreeReaderValue<Bool_t> myHlt2(reader,"Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS");
  // TTreeReaderValue<Double_t> myBDT(reader,"Bc_ISOLATION_BDT");
  // TTreeReaderValue<Double_t> myBDT(reader,"Jpsi_MM");

  //TTreeFormula* PidBinFormula = new TTreeFormula("PidBin","PidBin",tree);
  TTreeFormula* PFormula = new TTreeFormula("BachMu_P","BachMu_P",tree);
  TTreeFormula* EtaFormula = new TTreeFormula("BachMu_ETA","TMath::ATanH(BachMu_PZ/BachMu_P)",tree);
  TTreeFormula* nTracksFormula = new TTreeFormula("nTracks","nTracks",tree);
  TTreeFormula* Cut = new TTreeFormula("Cut",CutString,tree);
  TTreeFormula* PiCut = new TTreeFormula("PiCut",PiCutString,tree);
  TTreeFormula* KCut  = new TTreeFormula("KCut", KCutString, tree);
  TTreeFormula* PCut  = new TTreeFormula("PCut", PCutString, tree);
  TTreeFormula* MuCut = new TTreeFormula("MuCut",MuCutString,tree);
  TTreeFormula* ECut  = new TTreeFormula("ECut", ECutString, tree);

  Long64_t nentries = tree->GetEntries();
  for (Long64_t entry = 0; entry < nentries; ++entry) {
    if (entry % 5000 == 0) {
      std::cout << "ON ENTRY " << entry << " OUT OF " << nentries << std::endl;
    }
    tree->GetEntry(entry);
    if (Cut->EvalInstance()) {
      Double_t picut = PiCut->EvalInstance();
      Double_t kcut = KCut->EvalInstance();
      Double_t pcut = PCut->EvalInstance();
      Double_t mucut = MuCut->EvalInstance();
      Double_t ecut = ECut->EvalInstance();
      Species sp = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      Int_t hh = -1;
      if (PiCut->EvalInstance()) {
	hh = 0;
	MisIDFunction::listPi.Enter(entry);
	sp.piNum = 1.0;
      } else if (KCut->EvalInstance()) {
	hh = 1;
	MisIDFunction::listK.Enter(entry);
	sp.kNum = 1.0;
      } else if (PCut->EvalInstance()) {
	hh = 2;
	MisIDFunction::listP.Enter(entry);
	sp.pNum = 1.0;
      } else if (MuCut->EvalInstance()) {
	hh = 3;
	MisIDFunction::listMu.Enter(entry);
	sp.muNum = 1.0;
      } else if (ECut->EvalInstance()) {
	hh = 4;
	MisIDFunction::listE.Enter(entry);
	sp.eNum = 1.0;
      } else {
	hh = 5;
	MisIDFunction::listNone.Enter(entry);
	sp.uNum = 1.0;
      }
      //Int_t PidBin = static_cast<Int_t>(PidBinFormula->EvalInstance());
      Double_t p = PFormula->EvalInstance();
      Double_t eta = EtaFormula->EvalInstance();
      Double_t nt = nTracksFormula->EvalInstance();
      Int_t PidBin = HistoMatrix[0][0]->FindFixBin(p,eta,nt);
      if(PidBin==62){std::cout<<"HELP, THIS IS 62"<<std::endl;}
      //if (PidBin != 296) continue;
      auto it = MisIDFunction::countCache->find(PidBin);
      if (it == MisIDFunction::countCache->end()) {
	MisIDFunction::countCache->insert(std::pair<Int_t,Species>(PidBin,sp));
      } else {
	if (hh == 0)
	  it->second.piNum += 1.0;
	else if (hh == 1)
	  it->second.kNum += 1.0;
	else if (hh == 2)
	  it->second.pNum += 1.0;
	else if (hh == 3)
	  it->second.muNum += 1.0;
	else if (hh == 4)
	  it->second.eNum += 1.0;
	else
	  it->second.uNum += 1.0;
      }
    }
  }
}

void MisIDWeightor() {
  TChain* oldtree;
  TString oldfilename, newfilename;

  if(year=="16"){
    oldfilename = "/home/ejiang/tuples16/2016_Data_MD_misID_0_skimmed.root";
    newfilename = "/home/ejiang/tuples16/2016_Data_MD_misID_0_skimmed_unfold.root";
    oldtree = new TChain("DecayTree");
  }
  if(year=="17"){
    oldfilename = "/home/ejiang/tuples17/2017_Data_MU_misID_skimmed.root";
    newfilename = "/home/ejiang/tuples17/2017_Data_MU_misID_skimmed_unfold.root";
    oldtree = new TChain("JpsiRecTuple/DecayTree");
  }
  else if(year=="18"){
    oldfilename = "/home/ejiang/tuples18/2018_Data_MD_misID_0_skimmed.root";
    newfilename = "/home/ejiang/tuples18/2018_Data_MD_misID_0_skimmed_unfold.root";
    oldtree = new TChain("JpsiRecTuple/DecayTree");
  }
  // TString oldfilename = "/home/ejiang/tuples17/2017_Data_MU_misID_skimmed.root";
  oldtree->Add(oldfilename.Data());

  // TString newfilename = "/home/ejiang/tuples17/2017_Data_MU_misID_skimmed_unfold.root";
  std::cin.get();

  MisIDFunction* fobj = new MisIDFunction(0);
  fobj->createCountCache(oldtree);

  TF1 PIDWfunc("PIDWfunc",fobj,0,10000,5,"MisIDFunction");
  PIDWfunc.SetParameter(1,kFALSE);
  PIDWfunc.SetParameter(2,kFALSE);
  PIDWfunc.SetParameter(3,kFALSE);

  TF1 PIDWErrfunc("PIDWErrfunc",fobj,0,10000,5,"MisIDFunction");
  PIDWErrfunc.SetParameter(1,kTRUE);
  PIDWErrfunc.SetParameter(2,kFALSE);
  PIDWErrfunc.SetParameter(3,kFALSE);

  TF1* PIDWfunc_breakdown[numClasses];
  for (Int_t i = 0; i < numClasses; ++i) {
    PIDWfunc_breakdown[i] = new TF1("PIDWfunc_"+particle_name[i],fobj,0,10000,5,"MisIDFunction");
    PIDWfunc_breakdown[i]->SetParameter(1,kFALSE);
    PIDWfunc_breakdown[i]->SetParameter(2,kFALSE);
    PIDWfunc_breakdown[i]->SetParameter(3,kTRUE);
    PIDWfunc_breakdown[i]->SetParameter(4,i);
  }

  //Get old file, old tree and set top branch address
  if (kTRUE) {
    Int_t PidBin = 0;
    Double_t BachMu_P;
    Double_t BachMu_PZ;
    Double_t BachMu_ETA;
    Double_t nTracks;
    Int_t HadronHypo = 0;

    // umd, nominal
    Double_t PIDW = 0.0;
    Double_t PIDWVar = 0.0;
    Double_t PIDW_PI = 0.0;
    Double_t PIDW_K = 0.0;
    Double_t PIDW_P = 0.0;
    Double_t PIDW_E = 0.0;
    Double_t PIDW_MU = 0.0;
    Double_t PIDW_U = 0.0;

    //oldtree->SetBranchAddress("PidBin",&PidBin);
    oldtree->SetBranchAddress("BachMu_P",&BachMu_P);
    oldtree->SetBranchAddress("BachMu_PZ",&BachMu_PZ);
    oldtree->SetBranchAddress("nTracks",&nTracks);
    oldtree->SetBranchAddress("HadronHypo",&HadronHypo);
    oldtree->SetBranchAddress("PIDW",&PIDW);
    oldtree->SetBranchAddress("PIDWVar",&PIDWVar);
    ///    oldtree->SetBranchAddress("PIDW_DLL",&PIDW_DLL);
    ///    oldtree->SetBranchAddress("PIDWVar_DLL",&PIDWVar_DLL);
    
    //Create a new file + a clone of old tree in new file
    //TString newfilename = dir + "Bc_MisID_Data_PID_DISJ_MICRO.root";
    TFile *newfile = new TFile(newfilename.Data(),"recreate");
    TTree *newtree = oldtree->CloneTree(0);
    newtree->Branch("PIDW",&PIDW,"PIDW/D");
    newtree->Branch("PIDWVar",&PIDWVar,"PIDWVar/D");
    newtree->Branch("PidBin",&PidBin,"PidBin/I");
    newtree->Branch("PIDW_PI",&PIDW_PI,"PIDW_PI/D");
    newtree->Branch("PIDW_K",&PIDW_K,"PIDW_K/D");
    newtree->Branch("PIDW_P",&PIDW_P,"PIDW_P/D");
    newtree->Branch("PIDW_E",&PIDW_E,"PIDW_E/D");
    newtree->Branch("PIDW_MU",&PIDW_MU,"PIDW_MU/D");
    newtree->Branch("PIDW_U",&PIDW_U,"PIDW_U/D");

    Int_t MISIDTYPE;
    TBranch* b_MISIDTYPE = newtree->Branch("MISIDTYPE",&MISIDTYPE,"MISIDTYPE/I");
    
    auto initialize = [&] (Int_t n) {
      PIDWfunc.SetParameter(0,n);
      PIDWErrfunc.SetParameter(0,n);
      for (int i = 0; i < numClasses; ++i) {
        PIDWfunc_breakdown[i]->SetParameter(0,n);
      }
    };

    TFile pidf(pidfile_name,"READ");
    TH3F* bingetter = loadRootObject(pidf, "H2H_Pi_TO_Pi");
    auto getPidBin = [&] () {
      PidBin = bingetter->FindFixBin(BachMu_P,TMath::ATanH(BachMu_PZ/BachMu_P),nTracks);
    };
    
    auto compute = [&] () {
      //std::cout << "RUNNING UMD: VALUE" << std::endl;
      PIDW = PIDWfunc(PidBin);
      PIDWVar = PIDWErrfunc(PidBin);
      PIDW_PI = PIDWfunc_breakdown[0]->operator()(PidBin);
      PIDW_K = PIDWfunc_breakdown[1]->operator()(PidBin);
      PIDW_P = PIDWfunc_breakdown[2]->operator()(PidBin);
      PIDW_E = PIDWfunc_breakdown[3]->operator()(PidBin);
      PIDW_MU = PIDWfunc_breakdown[4]->operator()(PidBin);
      PIDW_U = PIDWfunc_breakdown[5]->operator()(PidBin);
      ///cout << PIDW <<","<<PIDWVar<<","<<PIDW_PI<<","<<PIDW_K<<","<<PIDW_P<<","<<PIDW_E<<","<<PIDW_MU<<","<<PIDW_U<<endl;
    };

    Int_t N_entries = oldtree->GetEntries();
    cout << "oldtree GetEntries = " << N_entries << endl;
    for (Int_t i = 0; i < N_entries; i++) {
      oldtree->GetEntry(i);
      getPidBin();
     
      PIDW = 0.0;
      PIDWVar = 0.0;
      
      // PIONS
      if (num >= 1) {
	if (MisIDFunction::listPi.Contains(i)) {
	  initialize(0);
	  compute();
	  HadronHypo = 0;
	  MISIDTYPE = 0;
	}
      }
       
      // KAONS
      if (num >= 2) {
	if (MisIDFunction::listK.Contains(i)) {
	  initialize(1);
	  compute();
	  HadronHypo = 1;
	  MISIDTYPE = 0;
	}
      }

      // PROTONS
      if (num >= 3) {
	if (MisIDFunction::listP.Contains(i)) {
	  initialize(2);
	  compute();
	  HadronHypo = 2;
	  MISIDTYPE = 0;
	}
      }
      
      // MUONS
      if (num >= 4) {
	if (MisIDFunction::listMu.Contains(i)) {
	  initialize(3);
	  compute();
	  HadronHypo = 3;
	  MISIDTYPE = 0;
	}
      }

      // ELECTRONS
      if (num >= 5) {
	if (MisIDFunction::listE.Contains(i)) {
	  initialize(4);
	  compute();
	  HadronHypo = 4;
	  MISIDTYPE = 0;
	}
      }

      // GHOSTS
      if (MisIDFunction::listNone.Contains(i)) {
	initialize(5);
	compute();
	HadronHypo = 5;
	MISIDTYPE = 1;
      }
      
      // FILL
      newtree->Fill();
      if (i % 10000 == 0) {cout << "Trying to fill tree!" << endl;}
    }

    // PRINT FRACTIONS
    fobj->printAggregateResults();
    newtree->AutoSave();
    delete newfile;
  }
}

#ifndef __CINT__
int main () {
#ifdef GHOST_ALTERNATE
  std::cout << "OH NO!" << std::endl;
  return 0;
#endif
  std::cout << "STARTING" << std::endl;
  MisIDWeightor();
  return 0;
}  // Main program when run stand-alone
#endif
