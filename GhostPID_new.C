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

void GhostPID_new()
{
    TChain chain("JpsiRecTuple/DecayTree");
    chain.Add("/afs/cern.ch/user/z/ziyang/my_eos/RJpsi/tuples/2016_MisID_MC_Jpsi_MD.root");
    chain.Add("/afs/cern.ch/user/z/ziyang/my_eos/RJpsi/tuples/2016_MisID_MC_Jpsi_MU.root");
    double xbins[7]={3e3,6e3,10e3,15.6e3,27e3,60e3,100e3};
    double etabins[3]={1.7,3.6,5};
    double trkbins[3]={0,200,600};

    TH3F *num=new TH3F("H2H_Ghost_TO_Mu","ghost muIDefficiency",6,xbins,2,etabins,2,trkbins);
    TH3F *num2=new TH3F("EFFNN0_Ghost_TO_Mu","ghost muIDefficiency",6,xbins,2,etabins,2,trkbins);
    TH3F *pi=new TH3F("H2H_Ghost_TO_Pi","ghost piIDefficiency",6,xbins,2,etabins,2,trkbins);
    TH3F *k=new TH3F("H2H_Ghost_TO_K","ghost kIDefficiency",6,xbins,2,etabins,2,trkbins);
    TH3F *p=new TH3F("H2H_Ghost_TO_P","ghost pIDefficiency",6,xbins,2,etabins,2,trkbins);
    TH3F *e=new TH3F("H2H_Ghost_TO_e","ghost eIDefficiency",6,xbins,2,etabins,2,trkbins);
    TH3F *g=new TH3F("H2H_Ghost_TO_Ghost","ghost gIDefficiency",6,xbins,2,etabins,2,trkbins);
    TH3F *total=new TH3F("denominator","denominator",6,xbins,2,etabins,2,trkbins);

    num->Sumw2();
    num2->Sumw2();
    pi->Sumw2();
    k->Sumw2();
    p->Sumw2();
    e->Sumw2();
    g->Sumw2();

    //g->g
    chain.Draw("nTracks:BachMu_ETA:BachMu_P>>H2H_Ghost_TO_Ghost","BachMu_TRUEID==0 && BachMu_ProbNNghost>0.2 && BachMu_isMuon==0 && Bc_BKGCAT==60 && Jpsi_BKGCAT==0");
    //g->e
    chain.Draw("nTracks:BachMu_ETA:BachMu_P>>H2H_Ghost_TO_e","BachMu_TRUEID==0 && !(BachMu_ProbNNghost>0.2) && BachMu_ProbNNe>0.1 && BachMu_isMuon==0 && Bc_BKGCAT==60 && Jpsi_BKGCAT==0");
    //g->mu
    chain.Draw("nTracks:BachMu_ETA:BachMu_P>>H2H_Ghost_TO_Mu","BachMu_TRUEID==0 && !(BachMu_ProbNNghost>0.2) && BachMu_PIDmu>2 && !(BachMu_ProbNNe>0.1) && BachMu_isMuon==1 && Bc_BKGCAT==60 && Jpsi_BKGCAT==0");
    chain.Draw("nTracks:BachMu_ETA:BachMu_P>>EFFNN0_Ghost_TO_Mu","BachMu_TRUEID==0 && !(BachMu_ProbNNghost>0.2) && BachMu_isMuon==1.0 && BachMu_ProbNNmu>0.5 && Bc_BKGCAT==60 && Jpsi_BKGCAT==0");
    //g->k
    chain.Draw("nTracks:BachMu_ETA:BachMu_P>>H2H_Ghost_TO_K","BachMu_TRUEID==0 && !(BachMu_ProbNNghost>0.2) && !(BachMu_PIDmu>2) && !(BachMu_ProbNNe>0.1) && BachMu_ProbNNk>0.1 && (BachMu_PIDK-BachMu_PIDp>2) && BachMu_isMuon==0 && Bc_BKGCAT==60 && Jpsi_BKGCAT==0");
    //g->p
    chain.Draw("nTracks:BachMu_ETA:BachMu_P>>H2H_Ghost_TO_P","BachMu_TRUEID==0 && !(BachMu_ProbNNghost>0.2) && !(BachMu_ProbNNe>0.1) && !(BachMu_PIDmu>2) && !(BachMu_ProbNNk>0.1) && !(BachMu_PIDK-BachMu_PIDp>2) && BachMu_ProbNNp>0.1 && BachMu_isMuon==0 && Bc_BKGCAT==60 && Jpsi_BKGCAT==0");
    //g->pi
    chain.Draw("nTracks:BachMu_ETA:BachMu_P>>H2H_Ghost_TO_Pi","BachMu_TRUEID==0 && !(BachMu_ProbNNghost>0.2) && !(BachMu_ProbNNe>0.1) && !(BachMu_PIDmu>2) && !(BachMu_ProbNNk>0.1) && !(BachMu_ProbNNp>0.1 && BachMu_PIDK-BachMu_PIDp>2) && BachMu_isMuon==0 && Bc_BKGCAT==60 && Jpsi_BKGCAT==0");
    //g total
    chain.Draw("nTracks:BachMu_ETA:BachMu_P>>denominator","BachMu_TRUEID==0 && Bc_BKGCAT==60 && Jpsi_BKGCAT==0");

    num->Divide(num,total,1,1,"B");
    num2->Divide(num2,total,1,1,"B");
    pi->Divide(pi,total,1,1,"B");
    k->Divide(k,total,1,1,"B");
    p->Divide(p,total,1,1,"B");
    e->Divide(e,total,1,1,"B");
    g->Divide(g,total,1,1,"B");
    
    num->Draw("ebox");
    cin.get();    

    for(int i=0; i < 8*4*4; i++)
    {
        if(!num->IsBinUnderflow(i) && !num->IsBinOverflow(i))
        {
            cout.precision(4);
            cout << std::fixed;
            cout << num2->GetBinContent(i) << "+/-" << num->GetBinError(i) << "\t\t";
            cout << num->GetBinContent(i) << "+/-" << num->GetBinError(i) << "\t\t";
            cout << pi->GetBinContent(i) << "+/-" << pi->GetBinError(i) << "\t\t";
            cout << k->GetBinContent(i) << "+/-" << k->GetBinError(i) << "\t\t";
            cout << p->GetBinContent(i) << "+/-" << p->GetBinError(i) << "\t\t";
            cout << e->GetBinContent(i) << "+/-" << e->GetBinError(i) << endl;
            cout << g->GetBinContent(i) << "+/-" << g->GetBinError(i) << endl;
        }
        if(i==0)
        {
            cout << "g->mu eff" << '\t'<<'\t';
            cout << "g->mu" << '\t'<<'\t';
            cout << "g->pi" << "\t\t"<<'\t';
            cout << "g->k" << "\t\t"<<'\t';
            cout << "g->p" << "\t\t"<<'\t';
            cout << "g->e" << "\t\t"<<'\t';
            cout << "g->g" << endl;
        }
    }

    TFile q("PerfHists_Ghost.root","RECREATE");
    q.Add(num);
    q.Add(num2);
    q.Add(pi);
    q.Add(k);
    q.Add(p);
    q.Add(e);
    q.Add(g);
    num->Write();
    num2->Write();
    pi->Write();
    k->Write();
    p->Write();
    e->Write();
    g->Write();
    ///q.Close();

}

