/**********************************************************************
 Created on : 27/07/2025
 Purpose    : Add clusters of mutiple passes
 Author     : Indranil Das, Research Associate, Imperial
 Email      : indranil.das@cern.ch | indra.ehep@gmail.com
**********************************************************************/
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH2.h>
#include <TF1.h>
#include "TEfficiency.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TTree.h"

int gindex = 1;

int addClusterPasses()
{
  int sl = 16;
  float sidelength = float(sl)/1000. ;
  float cir_rad = sidelength * TMath::Cos(TMath::Pi()/6.);
  int GetHexEdges(bool hextype, std::pair<float,float> xyCentre, float R, std::pair<float,float> []);//std::vector<std::pair<float,float>>& edglst);
  bool GetHistoBin(TH2F *h, float eta, float et, int& binX, int& binY);


  ////////////////// ============= cuts ===========///////////////////////
  float etaMin = 1.321;
  float etaMax = 3.152;
  
  float innerEtaMin = 1.7;
  float innerEtaMax = 2.9;
  
  double tcRoZCut = 0.15;
  double delRTh = 0.25;

  double pt_clusThresh = 150.;
  ////////////////// ============= cuts ===========///////////////////////

  ////////////////// ============= binning ===========///////////////////////
  const double effPtBin[29] = {0, 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 300., 400., 600.};
  const int nJetEtaBins = 6;
  Float_t jetEtaBin[nJetEtaBins+1] = {1.321, 1.7, 2.0, 2.3, 2.6, 2.9, 3.152} ;
  const int nJetPtBins = 43;
  Float_t jetPtBin[44] = {0, 15., 16., 18., 20., 22., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100.,  
                       110., 120., 130., 140., 150., 160., 170., 180., 190., 200.,
	               220., 240., 260., 280., 300.,
	               330., 360., 390., 420.,
                       440., 480., 520.};
  ////////////////// ============= binning ===========///////////////////////
  
  ////////////////// ============= histograms/functions ===========///////////////////////
  TH2F *hJetEtaPtBin = new TH2F("hJetEtaPtBin","hJetEtaPtBin", nJetEtaBins, jetEtaBin, nJetPtBins, jetPtBin);
  TEfficiency* effTrigGen = new TEfficiency("effTrigGen","Trigger efficiency;p_{T}(GeV);#epsilon",600,0,600.);
  TEfficiency* effTrigGenTDR = new TEfficiency("effTrigGenTDR","Trigger efficiency;p_{T}(GeV);#epsilon",28,effPtBin);
  TEfficiency* effTrigGenTDR_UnCorr = new TEfficiency("effTrigGenTDR_UnCorr","Trigger efficiency;p_{T}(GeV);#epsilon",28,effPtBin);
  TEfficiency* effTrigGenTDR_UnCorr_ClusSum = new TEfficiency("effTrigGenTDR_UnCorr_ClusSum","Trigger efficiency;p_{T}(GeV);#epsilon",28,effPtBin);
  TEfficiency* effTrigGenTDR_2DCorr_ClusSum = new TEfficiency("effTrigGenTDR_2DCorr_ClusSum","Trigger efficiency;p_{T}(GeV);#epsilon",28,effPtBin);
  TEfficiency* effTrigGenTDR_ClusSum = new TEfficiency("effTrigGenTDR_ClusSum","Trigger efficiency;p_{T}(GeV);#epsilon",28,effPtBin);

  TH1D *hNofClus = new TH1D("hNofClus","hNofClus", 400, -0.5, 399.5);
  TH1D *hNofClus3GeV = new TH1D("hNofClus3GeV","hNofClus3GeV", 400, -0.5, 399.5);  
  TH1D *hNassoClus = new TH1D("hNassoClus","hNassoClus", 20, -0.5, 19.5);
  TH1D *hclusDist = new TH1D("hclusDist","hclusDist", 400,  0.0, 0.6);  
  TH1D *hclosestClusE = new TH1D("hclosestClusE","hclosestClusE", 200, 0, 200);
  TH1D *hclosestClusE_UnCorr = new TH1D("hclosestClusE_UnCorr","hclosestClusE_UnCorr", 200, 0, 200);
  TH1D *hclosestClusE_passgt0 = new TH1D("hclosestClusE_passgt0","hclosestClusE_passgt0", 200, 0, 200);
  TH1D *hclosestClusE_corr2D_passgt0 = new TH1D("hclosestClusE_corr2D_passgt0","hclosestClusE_corr2D_passgt0", 200, 0, 200);
  TH1D *hclosestClusE_passgt0_corr2D = new TH1D("hclosestClusE_passgt0_corr2D","hclosestClusE_passgt0_corr2D", 200, 0, 200);
  
  
  float par0[nJetEtaBins] = {1.07956, 1.2916, 1.18022, 1.14248, 1.12945, 0.976157};
  float par1[nJetEtaBins] = {68.9515, 51.967, 52.3593, 47.0961, 39.7107, 32.1626};
  float par2[nJetEtaBins] = {18.5932, 5.89004, 8.34816, 10.0994, 9.27638, 16.5888};
  TF1* fClusEtaCorr[nJetEtaBins];
  for(int ieta=0;ieta<nJetEtaBins;ieta++) {
    fClusEtaCorr[ieta] = new TF1(Form("fClusEtaCorr_%d",ieta),"[0]+[1]/(x+[2])",1.0,500);
    fClusEtaCorr[ieta]->SetParameters(par0[ieta], par1[ieta], par2[ieta]); 
  }
  ////////////////// ============= histograms/functions ===========///////////////////////
  
  //const char *infile16 = "/home/indra/temp/multiple_passes/TPGS2Emu_tree_ntuples_16_0_vbf.root";
  //const char *infile16 = "/home/indra/temp/multiple_passes/TPGS2Emu_tree_ntuples_16_0_vbftc.root";
  //const char *infile16 = "/home/indra/temp/multiple_passes/TPGS2Emu_tree_ntuples_30_0_vbftc.root";
  //const char *infile16 = "/home/indra/temp/multiple_passes/TPGS2Emu_tree_ntuples_45_0_vbftc.root";
  //const char *infile16 = Form("/home/indra/temp/multiple_passes/TPGS2Emu_tree_ntuples_%d_0_pion-ideal.root",sl) ;
  //const char *infile16 = Form("/home/indra/temp/multiple_passes/TPGS2Emu_tree_ntuples_%d_0_vbf.root",sl) ;
  //const char *infile16 = "/home/indra/temp/multiple_passes/TPGReco_tree_ntuples_45_0_minbias.root";
  //const char *infile16 = "/home/indra/temp/multiple_passes/TPGReco_tree_ntuples_16_0.root";
  const char *infile16 = Form("/Data/root_files/stage2_emulation_results/Result_iter53/vbfHInv_0PU/TPGS2Emu_tree_ntuples_%d_merged.root",sl) ;
  const char *infile16_hist = Form("/Data/root_files/stage2_emulation_results/Result_iter53/vbfHInv_0PU/stage2SemiEmulator_ntuples_%d_merged.root",sl) ;
  std::unique_ptr<TFile> fin16(TFile::Open(Form("%s",infile16)));
  std::unique_ptr<TFile> fin16_hist(TFile::Open(Form("%s",infile16_hist)));
  std::unique_ptr<TTree> tr16((TTree*)fin16->Get("TPGS2Emu"));
  //std::unique_ptr<TTree> tr16((TTree*)fin16->Get("TPG_Reco"));
  
  tr16->SetBranchStatus("*",0);

  ////////////////////////////////////////////////////
  std::vector<float>  *genjet_pt = 0 ;
  tr16->SetBranchStatus("genjet_pt",1);
  tr16->SetBranchAddress("genjet_pt" , &genjet_pt);
  
  std::vector<float>  *genjet_eta = 0 ;
  tr16->SetBranchStatus("genjet_eta",1);
  tr16->SetBranchAddress("genjet_eta" , &genjet_eta);
  
  std::vector<float>  *genjet_phi = 0 ;
  tr16->SetBranchStatus("genjet_phi",1);
  tr16->SetBranchAddress("genjet_phi" , &genjet_phi);
  
  ////////////////////////////////////////////////////
  std::vector<float>  *tctot_ptsum = 0 ;
  tr16->SetBranchStatus("tctot_ptsum",1);
  tr16->SetBranchAddress("tctot_ptsum" , &tctot_ptsum);

  std::vector<float>  *tctot_xbyz = 0 ;
  tr16->SetBranchStatus("tctot_xbyz",1);
  tr16->SetBranchAddress("tctot_xbyz" , &tctot_xbyz);

  std::vector<float>  *tctot_ybyz = 0 ;
  tr16->SetBranchStatus("tctot_ybyz",1);
  tr16->SetBranchAddress("tctot_ybyz" , &tctot_ybyz);
  
  std::vector<float>  *tctot_eta = 0 ;
  tr16->SetBranchStatus("tctot_eta",1);
  tr16->SetBranchAddress("tctot_eta" , &tctot_eta);

  ////////////////////////////////////////////////////
  std::vector<float>  *tci_pt = 0 ;
  tr16->SetBranchStatus("tci_pt",1);
  tr16->SetBranchAddress("tci_pt" , &tci_pt);

  std::vector<float>  *tci_x = 0 ;
  tr16->SetBranchStatus("tci_x",1);
  tr16->SetBranchAddress("tci_x" , &tci_x);

  std::vector<float>  *tci_y = 0 ;
  tr16->SetBranchStatus("tci_y",1);
  tr16->SetBranchAddress("tci_y" , &tci_y);

  std::vector<float>  *tci_z = 0 ;
  tr16->SetBranchStatus("tci_z",1);
  tr16->SetBranchAddress("tci_z" , &tci_z);
  
  std::vector<float>  *tci_sect = 0 ;
  tr16->SetBranchStatus("tci_sect",1);
  tr16->SetBranchAddress("tci_sect" , &tci_sect);
  
  ////////////////////////////////////////////////////    
  std::vector<float>  *tcf_pt = 0 ;
  tr16->SetBranchStatus("tcf_pt",1);
  tr16->SetBranchAddress("tcf_pt" , &tcf_pt);

  std::vector<float>  *tcf_xbyz = 0 ;
  tr16->SetBranchStatus("tcf_xbyz",1);
  tr16->SetBranchAddress("tcf_xbyz" , &tcf_xbyz);

  std::vector<float>  *tcf_ybyz = 0 ;
  tr16->SetBranchStatus("tcf_ybyz",1);
  tr16->SetBranchAddress("tcf_ybyz" , &tcf_ybyz);
  
  std::vector<float>  *tcf_zcm = 0 ;
  tr16->SetBranchStatus("tcf_zcm",1);
  tr16->SetBranchAddress("tcf_zcm" , &tcf_zcm);
  
  std::vector<float>  *tcf_sect = 0 ;
  tr16->SetBranchStatus("tcf_sect",1);
  tr16->SetBranchAddress("tcf_sect" , &tcf_sect);  
  ////////////////////////////////////////////////////
  
  std::vector<float>  *clus_pt = 0 ;
  tr16->SetBranchStatus("clus_pt",1);
  tr16->SetBranchAddress("clus_pt" , &clus_pt);
  
  std::vector<float>  *clus_pt_corr1D = 0 ;
  tr16->SetBranchStatus("clus_pt_corr1D",1);
  tr16->SetBranchAddress("clus_pt_corr1D" , &clus_pt_corr1D);
  
  std::vector<float>  *clus_pt_corr2D = 0 ;
  tr16->SetBranchStatus("clus_pt_corr2D",1);
  tr16->SetBranchAddress("clus_pt_corr2D" , &clus_pt_corr2D);
  
  std::vector<float>  *clus_eta = 0 ;
  tr16->SetBranchStatus("clus_eta",1);
  tr16->SetBranchAddress("clus_eta" , &clus_eta);
  
  std::vector<float>  *clus_phi = 0 ;
  tr16->SetBranchStatus("clus_phi",1);
  tr16->SetBranchAddress("clus_phi" , &clus_phi);

  std::vector<float>  *clus_sect = 0 ;
  tr16->SetBranchStatus("clus_sect",1);
  tr16->SetBranchAddress("clus_sect" , &clus_sect);

  // std::vector<float>  *clus_zcm = 0 ;
  // tr16->SetBranchStatus("clus_zcm",1);
  // tr16->SetBranchAddress("clus_zcm" , &clus_zcm);
  
  std::vector<unsigned int>  *clus_pass = 0 ;
  tr16->SetBranchStatus("clus_pass",1);
  tr16->SetBranchAddress("clus_pass" , &clus_pass);  
  ////////////////////////////////////////////////////
  
  // TGraph *grGen = new TGraph(0);
  // TGraph *grClus0 = new TGraph(0);
  // TGraph *grClus1 = new TGraph(0);
  // TGraph *grClus2 = new TGraph(0);
  // TGraph *grClus3 = new TGraph(0);
  // TGraph *grClus4 = new TGraph(0);
  // TGraph *grClus5 = new TGraph(0);
  //Long64_t refEvent = 1;
  int refsect = 0;
  std::cout << std::fixed << std::setprecision(4) << std::setw(10) << std::endl;
  std::cout << "Nofevents: " << tr16->GetEntries() << std::endl;
  
  for (Long64_t ievent = 0 ; ievent < tr16->GetEntries() ; ievent++ ) {
    tr16->GetEntry(ievent) ;
    if(ievent%100000==0) std::cout << "ievent: " << ievent << ", nof genjets : " << genjet_pt->size() << std::endl;
    
    std::map<int,std::vector<int>> clus_asso;
    std::map<int,double> clus_totE, clus_totE_corr2D;
    int nofClus3GeV = 0;
    for(int iclus1 = 0; iclus1 < clus_pt->size() ; iclus1++){
      clus_asso[iclus1].resize(0);
      clus_totE[iclus1] = clus_pt->at(iclus1);
      clus_totE_corr2D[iclus1] = clus_pt_corr2D->at(iclus1);
      if(clus_pass->at(iclus1)!=0) continue;
      double clus_x1 = cos(clus_phi->at(iclus1))/sinh(clus_eta->at(iclus1));
      double clus_y1 = sin(clus_phi->at(iclus1))/sinh(clus_eta->at(iclus1));
      for(int iclus2 = 0; iclus2 < clus_pt->size() ; iclus2++){
	if(clus_pass->at(iclus2)<1) continue;
	double clus_x2 = cos(clus_phi->at(iclus2))/sinh(clus_eta->at(iclus2));
	double clus_y2 = sin(clus_phi->at(iclus2))/sinh(clus_eta->at(iclus2));      
	double dist = sqrt( (clus_x2-clus_x1)*(clus_x2-clus_x1) + (clus_y2-clus_y1)*(clus_y2-clus_y1) );
	hclusDist->Fill(dist);
	if(dist<2*sidelength and clus_pass->at(iclus2)==1){
	//if(dist< 0.04 and clus_pass->at(iclus2)>=1){
	//if(dist< 1.5*sidelength and clus_pass->at(iclus2)>=1 and clus_pass->at(iclus2)<=2){
	//if(dist< 2*sidelength and clus_pass->at(iclus2)>=1 and clus_pass->at(iclus2)<=2){
	//if(dist< 0.035 and clus_pass->at(iclus2)>=1){
	  clus_asso[iclus1].push_back(iclus2);
	  clus_totE[iclus1] = clus_totE[iclus1] + clus_pt->at(iclus2);
	  clus_totE_corr2D[iclus1] = clus_totE_corr2D[iclus1] + clus_pt_corr2D->at(iclus2);
	}
      }//clus2
      if(clus_pt->at(iclus1)>3.0) nofClus3GeV++;
      hNassoClus->Fill(clus_asso[iclus1].size());
    }//clus1
    hNofClus->Fill(clus_pt->size());
    hNofClus3GeV->Fill(nofClus3GeV);
    
    for(int igenj = 0; igenj < genjet_pt->size() ; igenj++){
      double genjet_x = cos(genjet_phi->at(igenj))/sinh(genjet_eta->at(igenj));
      double genjet_y = sin(genjet_phi->at(igenj))/sinh(genjet_eta->at(igenj));
      double minRoz = 1., closestClusE= 0., closestClusE_UnCorr= 0., closestClusE_passgt0 = 0., closestClusE_corr2D_passgt0 = 0., closestClusE_passgt0_corr2D = 0.;
      for(int iclus = 0; iclus < clus_pt->size() ; iclus++){
	if(clus_pt->at(iclus)<3.0) continue;
	if(genjet_eta->at(igenj)*clus_eta->at(iclus)<0) continue;
	double clus_x1 = cos(clus_phi->at(iclus))/sinh(clus_eta->at(iclus));
	double clus_y1 = sin(clus_phi->at(iclus))/sinh(clus_eta->at(iclus));
	double dClusGenRoz = sqrt( (clus_x1 - genjet_x)*(clus_x1 - genjet_x) + (clus_y1 - genjet_y)*(clus_y1 - genjet_y) );
	int binPtClus = -1, binEtaClus = -1;
	GetHistoBin(hJetEtaPtBin,clus_eta->at(iclus), clus_pt->at(iclus), binEtaClus, binPtClus);
	double ptcorr2D = fClusEtaCorr[binEtaClus-1]->Eval(clus_pt->at(iclus)) ;
	double ptcorr2D_tot = fClusEtaCorr[binEtaClus-1]->Eval(clus_totE[iclus]) ;
	if(dClusGenRoz<minRoz){
	  minRoz = dClusGenRoz;
	  closestClusE = clus_pt->at(iclus)*ptcorr2D;
	  closestClusE_UnCorr = clus_pt->at(iclus);
	  closestClusE_passgt0 = clus_totE[iclus];
	  closestClusE_corr2D_passgt0 = clus_totE_corr2D[iclus];
	  closestClusE_passgt0_corr2D = clus_totE[iclus]*ptcorr2D_tot;
	}	
      }
      
      bool hasClosestFound = (minRoz<1.)?true:false;
      bool hgcalInnerEta = (fabs(genjet_eta->at(igenj))>innerEtaMin and fabs(genjet_eta->at(igenj))<innerEtaMax) ? true : false;
      bool hgcalOuterEta = (fabs(genjet_eta->at(igenj))>etaMin and fabs(genjet_eta->at(igenj))<etaMax) ? true : false;
      bool hasClosestFoundPt = (hasClosestFound and closestClusE>pt_clusThresh)?true:false;
      bool hasClosestFoundPt_UnCorr = (hasClosestFound and closestClusE_UnCorr>pt_clusThresh)?true:false;
      bool hasClosestFoundPt_passgt0 = (hasClosestFound and closestClusE_passgt0>pt_clusThresh)?true:false;
      bool hasClosestFoundPt_corr2D_passgt0 = (hasClosestFound and closestClusE_corr2D_passgt0>pt_clusThresh)?true:false;
      bool hasClosestFoundPt_passgt0_corr2D = (hasClosestFound and closestClusE_passgt0_corr2D>pt_clusThresh)?true:false;
      
      if(hgcalOuterEta and hasClosestFound) {	
	if(hgcalInnerEta){
	  effTrigGen->Fill(hasClosestFoundPt,genjet_pt->at(igenj));
	  effTrigGenTDR->Fill(hasClosestFoundPt,genjet_pt->at(igenj));
	  effTrigGenTDR_UnCorr->Fill(hasClosestFoundPt_UnCorr,genjet_pt->at(igenj));
	  effTrigGenTDR_UnCorr_ClusSum->Fill(hasClosestFoundPt_passgt0,genjet_pt->at(igenj));
	  effTrigGenTDR_2DCorr_ClusSum->Fill(hasClosestFoundPt_corr2D_passgt0,genjet_pt->at(igenj));
	  effTrigGenTDR_ClusSum->Fill(hasClosestFoundPt_passgt0_corr2D,genjet_pt->at(igenj));
	  
	  hclosestClusE->Fill(closestClusE);
	  hclosestClusE_UnCorr->Fill(closestClusE_UnCorr);
	  hclosestClusE_passgt0->Fill(closestClusE_passgt0);
	  hclosestClusE_corr2D_passgt0->Fill(closestClusE_corr2D_passgt0);
	  hclosestClusE_passgt0_corr2D->Fill(closestClusE_passgt0_corr2D);
	}	
      }
      

    }//genj loop
    
    genjet_pt->clear();
    genjet_eta->clear();
    genjet_phi->clear();
    tctot_ptsum->clear();
    tctot_xbyz->clear();
    tctot_ybyz->clear();
    tctot_eta->clear();
    tci_pt->clear();
    tci_x->clear();
    tci_y->clear();
    tci_z->clear();
    tci_sect->clear();
    tcf_pt->clear();
    tcf_xbyz->clear();
    tcf_ybyz->clear();
    tcf_zcm->clear();
    tcf_sect->clear();
    clus_pt->clear();
    clus_pt_corr1D->clear();
    clus_pt_corr2D->clear();
    clus_eta->clear();
    clus_phi->clear();
    clus_sect->clear();
    //clus_zcm->clear();
    
  }//event loop

  TH1D *hSel_UnCorr = (TH1D *) effTrigGenTDR_UnCorr->GetCopyPassedHisto();
  TH1D *hTot_UnCorr = (TH1D *) effTrigGenTDR_UnCorr->GetCopyTotalHisto();
  hSel_UnCorr->SetName("hSel_UnCorr_new");
  hTot_UnCorr->SetName("hTot_UnCorr_new");


  TEfficiency *effTrigGenTDR_Uncorr_histfile = (TEfficiency *)fin16_hist->Get(effTrigGenTDR_UnCorr->GetName()) ;
  TH1D *hSel_UnCorr_hf = (TH1D *) effTrigGenTDR_Uncorr_histfile->GetCopyPassedHisto();
  TH1D *hTot_UnCorr_hf = (TH1D *) effTrigGenTDR_Uncorr_histfile->GetCopyTotalHisto();
  hSel_UnCorr_hf->SetName("hSel_UnCorr_old");
  hTot_UnCorr_hf->SetName("hTot_UnCorr_old");
  
  std::string outname = "output";
  TFile *fout = new TFile(Form("%s.root",outname.c_str()),"recreate");
  hSel_UnCorr->Write();
  hTot_UnCorr->Write();
  //////=== trigger efficiency histos ===////////
  effTrigGen->Write();
  effTrigGenTDR->Write();
  effTrigGenTDR_UnCorr->Write();
  effTrigGenTDR_UnCorr_ClusSum->Write();
  effTrigGenTDR_2DCorr_ClusSum->Write();
  effTrigGenTDR_ClusSum->Write();
  //////////////////////////////////////////////
  hSel_UnCorr_hf->Write();
  hTot_UnCorr_hf->Write();
  /////////////////////////////////////////////
  hNofClus->Write();
  hNofClus3GeV->Write();
  hNassoClus->Write();
  hclosestClusE->Write();
  hclosestClusE_UnCorr->Write();
  hclosestClusE_passgt0->Write();
  hclosestClusE_corr2D_passgt0->Write();
  hclosestClusE_passgt0_corr2D->Write();
  hclusDist->Write();
  fout->Close();
  delete fout;


  return true;
}

int GetHexEdges(bool hextype, std::pair<float,float> xyCentre, float R, std::pair<float,float> edges[6])//std::vector<std::pair<float,float>>& edglst)
{
  
  //hextype 1 for regular, 0 truncated type;
  float r = TMath::Cos(TMath::Pi()/6.)*R;
  float aby2 = TMath::Sin(TMath::Pi()/6.)*R;
  if(hextype){
    edges[0] = std::make_pair(xyCentre.first, xyCentre.second+R);
    edges[1] = std::make_pair(xyCentre.first+r, xyCentre.second+aby2);
    edges[2] = std::make_pair(xyCentre.first+r, xyCentre.second-aby2);
    edges[3] = std::make_pair(xyCentre.first, xyCentre.second-R);
    edges[4] = std::make_pair(xyCentre.first-r, xyCentre.second-aby2);
    edges[5] = std::make_pair(xyCentre.first-r, xyCentre.second+aby2);    
  }else{
    edges[0] = std::make_pair(xyCentre.first+aby2, xyCentre.second+r);
    edges[1] = std::make_pair(xyCentre.first+R, xyCentre.second);
    edges[2] = std::make_pair(xyCentre.first+aby2, xyCentre.second-r);
    edges[3] = std::make_pair(xyCentre.first-aby2, xyCentre.second-r);
    edges[4] = std::make_pair(xyCentre.first-R, xyCentre.second);
    edges[5] = std::make_pair(xyCentre.first-aby2, xyCentre.second+r);    
  }
  
  // edglst.clear();
  // edglst.push_back(edges[0]);
  
  return true;
}

double distanceXX(const double *ca, const double *cb){
  // std::cout << "Distance of " << ca[0] << ", " << ca[1]
  //	      << " from " << cb[0] << ", " << cb[1];
  
  double d(sqrt((ca[0] - cb[0]) * (ca[0] - cb[0]) +
		(ca[1] - cb[1]) * (ca[1] - cb[1])));
  
  // std::cout << " = " << d << std::flush << std::endl;
  return d;
}

bool GetHistoBin(TH2F *h, float eta, float et, int& binX, int& binY)
{
  eta = fabs(eta);
  double minX = h->GetXaxis()->GetBinCenter(1);
  double maxX = h->GetXaxis()->GetBinCenter(h->GetNbinsX());
  double minY = h->GetYaxis()->GetBinCenter(1);
  double maxY = h->GetYaxis()->GetBinCenter(h->GetNbinsY());

  binX = -1;
  if(eta <= minX)
    binX = 1;
  else if(eta >= maxX)
    binX = h->GetXaxis()->FindBin(maxX);
  else
    binX = h->GetXaxis()->FindBin(eta);

  binY = -1;
  if(et <= minY)
    binY = 1;
  else if(et >= maxY)
    binY = h->GetYaxis()->FindBin(maxY);
  else
    binY = h->GetYaxis()->FindBin(et);

  return true;
}
