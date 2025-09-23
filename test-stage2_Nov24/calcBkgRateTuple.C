/**********************************************************************
 Created on : 02/07/2025
 Purpose    : Calculate background rate from CMMSW ntuples
 Author     : Indranil Das, Research Associate, Imperial
 Email      : indranil.das@cern.ch | indra.ehep@gmail.com
**********************************************************************/

#include <iostream>
#include <fstream>
#include <map>

#include <TChain.h>
#include <TTree.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH2.h>
#include <TF1.h>
#include "TEfficiency.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraphErrors.h"

using namespace std;

int calcBkgRateTuple()
{
  
  const char *server = "root://eosuser.cern.ch/";
  // const char *infile = "/eos/cms/store/group/dpg_hgcal/comm_hgcal/TPG/stage2_emulator_ntuples_semiemulator_clusterProperties/minbias_PU200_2kFiles/ntuple_0.root" ;
  // std::unique_ptr<TFile> fin(TFile::Open(Form("%s/%s",server,infile)));
  // std::unique_ptr<TTree> tr((TTree*)fin->Get("l1tHGCalTriggerNtuplizer/HGCalTriggerNtuple"));
  //std::unique_ptr<TTree> tr((TTree*)fin->Get("l1tHGCalTriggerNtuplizer/HGCalTriggerNtuple"));
  
  TChain *trcmssw = new TChain("l1tHGCalTriggerNtuplizer/HGCalTriggerNtuple");
  // TChain *trcmssw = new TChain("TPG_Reco");
  // //const char *inputfname = "/tmp/remote_fs.txt";
  const char *inputfname = "/tmp/local_fs.txt";
  // //const char *inputfname = "/tmp/local_standalone_fs_16.txt";
  // const char *inputfname = "/tmp/local_standalone_fs_30.txt";
  std::string s;
  ifstream fin(inputfname);
  int nfiles = 0;
  while(std::getline(fin,s)){
    std::cout << "Filename : " << s << std::endl;
    //trcmssw->Add(Form("%s/%s",server,s.c_str()));
    trcmssw->Add(Form("%s",s.c_str()));
    nfiles++;
  }
  std::cout << "Nof added files : " << nfiles << std::endl;

  const char *infile16 = "/home/hep/idas/stage2_emulation_results/ResultBkg_iter42/minbias_PU200_2kFiles/TPGReco_tree_ntuples_16_merged.root";
  const char *infile30 = "/home/hep/idas/stage2_emulation_results/ResultBkg_iter42/minbias_PU200_2kFiles/TPGReco_tree_ntuples_30_merged.root";
  const char *infile45 = "/home/hep/idas/stage2_emulation_results/ResultBkg_iter42/minbias_PU200_2kFiles/TPGReco_tree_ntuples_45_merged.root";
  std::unique_ptr<TFile> fin16(TFile::Open(Form("%s",infile16)));
  std::unique_ptr<TFile> fin30(TFile::Open(Form("%s",infile30)));
  std::unique_ptr<TFile> fin45(TFile::Open(Form("%s",infile45)));

  std::unique_ptr<TTree> tr16((TTree*)fin16->Get("TPG_Reco"));
  std::unique_ptr<TTree> tr30((TTree*)fin30->Get("TPG_Reco"));
  std::unique_ptr<TTree> tr45((TTree*)fin45->Get("TPG_Reco"));
  
  trcmssw->SetBranchStatus("*",0);
  tr16->SetBranchStatus("*",0);
  tr30->SetBranchStatus("*",0);
  tr45->SetBranchStatus("*",0);  
  
  //////////////////// CMSSW ///////////////////////
  std::vector<float>  *clusE_a16_cmssw = 0 ;
  trcmssw->SetBranchStatus("cl3d_p16Tri_pt",1);
  trcmssw->SetBranchAddress("cl3d_p16Tri_pt" , &clusE_a16_cmssw);

  std::vector<float>  *clusE_a30_cmssw = 0 ;
  trcmssw->SetBranchStatus("cl3d_p03Tri_pt",1);
  trcmssw->SetBranchAddress("cl3d_p03Tri_pt" , &clusE_a30_cmssw);

  std::vector<float>  *clusE_a45_cmssw = 0 ;
  trcmssw->SetBranchStatus("cl3d_p045Tri_pt",1);
  trcmssw->SetBranchAddress("cl3d_p045Tri_pt" , &clusE_a45_cmssw);
  //////////////////// CMSSW ///////////////////////
  
  ///////////////////// 16 /////////////////////////
  std::vector<float>  *clusE_a16_uncorr = 0 ;
  tr16->SetBranchStatus("clus_pt",1);
  tr16->SetBranchAddress("clus_pt" , &clusE_a16_uncorr);
  
  std::vector<float>  *clusE_a16_corr1D = 0 ;
  tr16->SetBranchStatus("clus_pt_corr1D",1);
  tr16->SetBranchAddress("clus_pt_corr1D" , &clusE_a16_corr1D);
  
  std::vector<float>  *clusE_a16_corr2D = 0 ;
  tr16->SetBranchStatus("clus_pt_corr2D",1);
  tr16->SetBranchAddress("clus_pt_corr2D" , &clusE_a16_corr2D);
  
  std::vector<unsigned int>  *cluspass_a16 = 0 ;
  tr16->SetBranchStatus("clus_pass",1);
  tr16->SetBranchAddress("clus_pass" , &cluspass_a16);
  ///////////////////// 16 /////////////////////////
  
  ///////////////////// 30 /////////////////////////
  std::vector<float>  *clusE_a30_uncorr = 0 ;
  tr30->SetBranchStatus("clus_pt",1);
  tr30->SetBranchAddress("clus_pt" , &clusE_a30_uncorr);
  
  std::vector<float>  *clusE_a30_corr1D = 0 ;
  tr30->SetBranchStatus("clus_pt_corr1D",1);
  tr30->SetBranchAddress("clus_pt_corr1D" , &clusE_a30_corr1D);
  
  std::vector<float>  *clusE_a30_corr2D = 0 ;
  tr30->SetBranchStatus("clus_pt_corr2D",1);
  tr30->SetBranchAddress("clus_pt_corr2D" , &clusE_a30_corr2D);
  
  std::vector<unsigned int>  *cluspass_a30 = 0 ;
  tr30->SetBranchStatus("clus_pass",1);
  tr30->SetBranchAddress("clus_pass" , &cluspass_a30);  
  ///////////////////// 30 /////////////////////////
  
  ///////////////////// 45 /////////////////////////
  std::vector<float>  *clusE_a45_uncorr = 0 ;
  tr45->SetBranchStatus("clus_pt",1);
  tr45->SetBranchAddress("clus_pt" , &clusE_a45_uncorr);
  
  std::vector<float>  *clusE_a45_corr1D = 0 ;
  tr45->SetBranchStatus("clus_pt_corr1D",1);
  tr45->SetBranchAddress("clus_pt_corr1D" , &clusE_a45_corr1D);
  
  std::vector<float>  *clusE_a45_corr2D = 0 ;
  tr45->SetBranchStatus("clus_pt_corr2D",1);
  tr45->SetBranchAddress("clus_pt_corr2D" , &clusE_a45_corr2D);
  
  std::vector<unsigned int>  *cluspass_a45 = 0 ;
  tr45->SetBranchStatus("clus_pass",1);
  tr45->SetBranchAddress("clus_pass" , &cluspass_a45);  
  ///////////////////// 45 /////////////////////////


  // //////////////////// Find leading pt ///////////////////////
  // std::cout << "Total number of Events in CMSSW tree : " << trcmssw->GetEntries() << std::endl;
  Long64_t nofEvents = 318880 ;
  if(nofEvents==0) nofEvents = trcmssw->GetEntries();
  std::map<Long64_t,double> maxpt16_cmssw, maxpt30_cmssw, maxpt45_cmssw;
  std::map<Long64_t,double> maxptmap16, maxptmap30, maxptmap45;
  
  // for (Long64_t ievent = 0 ; ievent < nofEvents ; ievent++ ) {
  //   trcmssw->GetEntry(ievent) ;
    
  //   double maxpt = 0;
  //   for(int iclus = 0; iclus < clusE_a16_cmssw->size() ; iclus++){
  //     if(clusE_a16_cmssw->at(iclus)>maxpt and (cluspass_a16->at(iclus)==0 or cluspass_a16->at(iclus)==1)) maxpt = clusE_a16_cmssw->at(iclus);
  //   }
  //   maxpt16_cmssw[ievent] = maxpt;

  //   maxpt = 0;
  //   for(int iclus = 0; iclus < clusE_a30_cmssw->size() ; iclus++){
  //     if(clusE_a30_cmssw->at(iclus)>maxpt and (cluspass_a30->at(iclus)==0 or cluspass_a16->at(iclus)==1)) maxpt = clusE_a30_cmssw->at(iclus);
  //   }
  //   maxpt30_cmssw[ievent] = maxpt;

  //   maxpt = 0;
  //   for(int iclus = 0; iclus < clusE_a45_cmssw->size() ; iclus++){
  //     if(clusE_a45_cmssw->at(iclus)>maxpt and (cluspass_a45->at(iclus)==0 or cluspass_a16->at(iclus)==1)) maxpt = clusE_a45_cmssw->at(iclus);
  //   }
  //   maxpt45_cmssw[ievent] = maxpt;

    
  //   clusE_a16_cmssw->clear();
  //   clusE_a30_cmssw->clear();
  //   clusE_a45_cmssw->clear();
  // }
  
  //// Find leading pt ///////////////////////
  std::cout << "Total number of Events in tr16 tree : " << tr16->GetEntries() << std::endl;
  nofEvents = tr16->GetEntries() ;
  for (Long64_t ievent = 0 ; ievent < nofEvents ; ievent++ ) {
    tr16->GetEntry(ievent) ;
    tr30->GetEntry(ievent) ;
    tr45->GetEntry(ievent) ;
    
    double maxpt = 0;
    for(int iclus = 0; iclus < int(clusE_a16_uncorr->size()) ; iclus++){
      if(clusE_a16_uncorr->at(iclus)>maxpt ) maxpt = clusE_a16_uncorr->at(iclus);
    }
    maxptmap16[ievent] = maxpt;

    maxpt = 0;
    for(int iclus = 0; iclus < int(clusE_a30_uncorr->size()) ; iclus++){
      if(clusE_a30_uncorr->at(iclus)>maxpt and (cluspass_a30->at(iclus)==0 or cluspass_a16->at(iclus)==1)) maxpt = clusE_a30_uncorr->at(iclus);
    }
    maxptmap30[ievent] = maxpt;
    
    maxpt = 0;
    for(int iclus = 0; iclus < int(clusE_a45_uncorr->size()) ; iclus++){
      if(clusE_a45_uncorr->at(iclus)>maxpt and (cluspass_a45->at(iclus)==0 or cluspass_a16->at(iclus)==1)) maxpt = clusE_a45_uncorr->at(iclus);
    }
    maxptmap45[ievent] = maxpt;

    if(ievent%10000==0)
      std::cout<<"Event : "<< ievent
	       <<", nof Clusters16 : "<< clusE_a16_uncorr->size()
	       <<", nof Clusters30 : "<< clusE_a30_uncorr->size()
	       <<", nof Clusters45 : "<< clusE_a45_uncorr->size()
	       << ", maxpt: " << maxpt
	       << std::endl;
    
    clusE_a16_uncorr->clear();
    clusE_a30_uncorr->clear();
    clusE_a45_uncorr->clear();

    cluspass_a16->clear();
    cluspass_a30->clear();
    cluspass_a45->clear();
  }

  
  TGraphErrors *gr16 = new TGraphErrors(200);
  TGraphErrors *gr30 = new TGraphErrors(200);
  TGraphErrors *gr45 = new TGraphErrors(200);
  int ipoint = 0;
  for(int ithr = 20 ; ithr <= 200 ; ithr = ithr+6){
    int pass16 = 0, pass30 = 0, pass45 = 0;
    for (Long64_t ievent = 0 ; ievent < nofEvents ; ievent++ ) {
      if(maxptmap16[ievent]>double(ithr)) pass16++;
      if(maxptmap30[ievent]>double(ithr)) pass30++;
      if(maxptmap45[ievent]>double(ithr)) pass45++;
    }
    double rate16 = double(pass16) * 40.e3 * 2760. / 3564. / double(nofEvents) ;
    double rate30 = double(pass30) * 40.e3 * 2760. / 3564. / double(nofEvents) ;
    double rate45 = double(pass45) * 40.e3 * 2760. / 3564. / double(nofEvents) ;
    std::cout <<"Threshold :" << ithr << ", Rate16 : " << rate16 << ", Rate30 : " << rate30 << ", Rate45 : " << rate45 << std::endl;
    gr16->SetPoint(ipoint,double(ithr),rate16);
    gr30->SetPoint(ipoint,double(ithr),rate30);
    gr45->SetPoint(ipoint,double(ithr),rate45);
    
    double yval_err16 = rate16*fabs( sqrt(pass16)/double(pass16) - sqrt(nofEvents)/double(nofEvents) );
    double yval_err30 = rate30*fabs( sqrt(pass30)/double(pass30) - sqrt(nofEvents)/double(nofEvents) );
    double yval_err45 = rate45*fabs( sqrt(pass45)/double(pass45) - sqrt(nofEvents)/double(nofEvents) );
    gr16->SetPointError(ipoint, 0, yval_err16);
    gr30->SetPointError(ipoint, 0, yval_err30);
    gr45->SetPointError(ipoint, 0, yval_err45);
    ipoint++;
  }//threshold loop
  gr16->Set(ipoint);
  gr30->Set(ipoint);
  gr45->Set(ipoint);
  
  gr16->SetTitle("");
  gr30->SetTitle("");
  gr45->SetTitle("");
  
  gr16->SetFillColor(kBlack);
  gr30->SetFillColor(kBlue);
  gr45->SetFillColor(kRed);
  
  gr16->SetLineColor(kBlack);
  gr30->SetLineColor(kBlue);
  gr45->SetLineColor(kRed);
  
  gr16->SetLineWidth(3);
  gr30->SetLineWidth(3);
  gr45->SetLineWidth(3);
  
  gr16->SetMinimum(1);
  gr16->SetMaximum(5.e4);
  gr16->GetYaxis()->SetTitle("Single jet rate [kHz]");
  gr16->GetXaxis()->SetTitle("Threshold [GeV]");
  gr16->GetXaxis()->SetRangeUser(0,200);
  
  TLatex *texl = new TLatex(3.369606,58503.36,"CMS");
  texl->SetTextSize(0.035);
  TLatex *texp = new TLatex(26.9,59533.27,"Preliminary");
  texp->SetTextSize(0.025);
  texp->SetTextFont(52);
  auto leg0 = new TLegend(0.19,0.25,0.58,0.42);
  leg0->SetTextSize(0.028);
  leg0->SetHeader(Form("Minbias, PU 200"));
  leg0->AddEntry(gr16,"a = 0.016","lp");
  leg0->AddEntry(gr30,"a = 0.030","lp");
  leg0->AddEntry(gr45,"a = 0.045","lp");
  leg0->SetShadowColor(kWhite);
  
  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  c1->SetLogy();
  c1->SetGridx();
  c1->SetGridy();
  c1->SetTickx();
  c1->SetTicky();
  gr16->Draw("a3c");
  gr30->Draw("3c same");
  gr45->Draw("3c same");
  leg0->Draw();
  texl->Draw("same");
  texp->Draw("same");

  TFile *fout = new TFile("output.root","recreate");
  gr16->Write();
  gr30->Write();
  gr45->Write();
  c1->Write();
  fout->Close();
  
  maxptmap16.clear();
  maxptmap30.clear();
  maxptmap45.clear();

  // delete tr;
  // delete fout;
  
  return true;
}
