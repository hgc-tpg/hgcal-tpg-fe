/**********************************************************************
 Created on : 11/05/2025
 Purpose    : Calculate the jet resolutions
 Author     : Indranil Das, Research Associate, Imperial
 Email      : indranil.das@cern.ch | indra.ehep@gmail.com
**********************************************************************/
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <TROOT.h>
#include <TChain.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH2.h>
#include <TF1.h>
#include "TParticlePDG.h"
#include "TDatabasePDG.h"
#include "TEfficiency.h"

#include "Stage2.hh"

#include "TPGTCFloats.hh"
#include "TPGLSBScales.hh"

typedef struct{
  std::string name;
  int index;
}JetPart;

typedef struct{
  uint32_t tcid, layer;
  float z,eta;
}TCOrder;

bool comp(TCOrder x, TCOrder y) {
  float diff = x.z - y.z ;
  bool condn = ((x.layer<y.layer) or ((x.layer==y.layer) and (fabs(diff)>1.e-5)) or ((x.layer==y.layer) and (fabs(diff)<1.e-5) and (fabs(x.eta)<fabs(y.eta))));
  return condn ;
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

int main(int argc, char** argv)
{
  
  TPGLSBScales::TPGStage2ClusterLSB lsbScales;
  lsbScales.print();  
  // return true;
  
  bool doPrint = 0;
  float pt_clusThresh = 30.; //GeV electron 30 GeV //pion = 100, VBF = 150 GeV
  float pt_clusThresh_etaeff = 50.; //GeV electron 50 GeV //pion = 150, VBF = 200 GeV
  
  float pt_TCThresh = 10.; //GeV
  std::cout << "========= Run as : ./stage2SemiEmulator.exe $input_file $index $nofevents ==========" << std::endl;
  
  if(argc < 2){
    std::cerr << argv[0] << ": no input file has been specified specified" << std::endl;
    return false;
  }
  std::string inputfile;
  std::istringstream issInfile(argv[1]);
  issInfile >> inputfile;
  std::string index = "0";
  if(argc >= 3){
    std::istringstream issnIndex(argv[2]);
    issnIndex >> index;
  }
  uint64_t nofEvents(0);
  if(argc >= 4){
    std::istringstream issnEvents(argv[3]);
    issnEvents >> nofEvents;
  }
  
  float sidelength = 0.016;
  if(argc >= 5){
    std::istringstream issSidelength(argv[4]);
    issSidelength >> sidelength;
  }
  
  std::string outputfile_extn = "";
  if(argc >= 6){
    std::istringstream issOutputfileExtn(argv[5]);
    issOutputfileExtn >> outputfile_extn;
  }

  std::string sampleType = "";
  if(argc >= 7){
    std::istringstream issSampleType(argv[6]);
    issSampleType >> sampleType;
  }
  
  pt_clusThresh = 30.; //GeV electron 30 GeV //pion = 100, VBF = 150 GeV
  pt_clusThresh_etaeff = 50.; //GeV electron 50 GeV //pion = 150, VBF = 200 GeV
  if(sampleType.find("vbfHInv")!=std::string::npos){

    pt_clusThresh = 150.; 
    pt_clusThresh_etaeff = 200.;
    int sdlg = int(sidelength*1000) ;
    //===============================================
    // @10 kHz
    //===============================================
    //gr16_corr2D_passall: a = 0.016:: Lower that rate 10 kHz: 140, and higher : 146, Xval: 143.352
    //gr30_corr2D_passall: a = 0.03:: Lower that rate 10 kHz: 152, and higher : 158, Xval: 152.457
    //gr45_corr2D_passall: a = 0.045:: Lower that rate 10 kHz: 158, and higher : 164, Xval: 163.346
    //===============================================
    if(sdlg==16)
      pt_clusThresh = 143.;
    else if(sdlg==30)
      pt_clusThresh = 152.;
    else if(sdlg==45)
      pt_clusThresh = 163.; 
    pt_clusThresh_etaeff = pt_clusThresh + 50.;
    //===============================================
    
    //===============================================
    // @100 kHz
    //===============================================
    //gr16_corr2D_passall: a = 0.016:: Lower that rate 100 kHz: 86, and higher : 92, Xval: 90.8449
    //gr30_corr2D_passall: a = 0.03:: Lower that rate 100 kHz: 92, and higher : 98, Xval: 95.3578
    //gr45_corr2D_passall: a = 0.045:: Lower that rate 100 kHz: 110, and higher : 116, Xval: 115.74
    //===============================================
    // if(sdlg==16)
    //   pt_clusThresh = 91.;
    // else if(sdlg==30)
    //   pt_clusThresh = 95.;
    // else if(sdlg==45)
    //   pt_clusThresh = 116.; 
    // pt_clusThresh_etaeff = pt_clusThresh + 50.;
    //===============================================
    
    //===============================================
    // @50 kHz
    //===============================================
    // gr16_corr2D_passall: a = 0.016:: Lower that rate 50 kHz: 98, and higher : 104, Xval: 100.644
    // gr30_corr2D_passall: a = 0.03:: Lower that rate 50 kHz: 104, and higher : 110, Xval: 107.071
    // gr45_corr2D_passall: a = 0.045:: Lower that rate 50 kHz: 122, and higher : 128, Xval: 126.683
    //===============================================
    // if(sdlg==16)
    //   pt_clusThresh = 101.;
    // else if(sdlg==30)
    //   pt_clusThresh = 107.;
    // else if(sdlg==45)
    //   pt_clusThresh = 127.; 
    // pt_clusThresh_etaeff = pt_clusThresh + 50.;
    //===============================================

  }else if(sampleType.find("singlePion")!=std::string::npos){
    // pt_clusThresh = 150.; 
    // pt_clusThresh_etaeff = 200.; 
    pt_clusThresh = 100.; 
    pt_clusThresh_etaeff = 150.; 
  }
  
  std::unique_ptr<TFile> fin( TFile::Open(Form(inputfile.c_str())) );
  //TFile *fin = TFile::Open(Form(inputfile.c_str()));
  std::cout<<"Input: filename : " << fin->GetName() << std::endl;
  std::cout << "Input: outfile index : " << index << std::endl;
  std::cout << "Input: nofEvents : " << nofEvents << std::endl;
  std::cout << "Input: sidelength : " << sidelength << std::endl;
  std::cout << "Input: MC Sample : " << sampleType << std::endl;
  std::cout << "pt_clusThresh for L1 : " << pt_clusThresh << std::endl;
  std::cout << "pt_clusThresh_etaeff for Trigger eta eff : " << pt_clusThresh_etaeff << std::endl;
  
  //TTree *tr  = (TTree*)fin->Get("l1tHGCalTriggerNtuplizer/HGCalTriggerNtuple");
  std::unique_ptr<TTree> tr((TTree*)fin->Get("l1tHGCalTriggerNtuplizer/HGCalTriggerNtuple"));  
  tr->SetBranchStatus("*",0);

  Int_t event ;
  tr->SetBranchStatus("event",1);
  tr->SetBranchAddress("event"	, &event);

  float vtx_z;
  tr->SetBranchStatus("vtx_z",1);
  tr->SetBranchAddress("vtx_z" , &vtx_z);
  
  Int_t gen_n = 0 ;
  tr->SetBranchStatus("gen_n",1);
  tr->SetBranchAddress("gen_n" , &gen_n);
  
  std::vector<float>  *gen_pt = 0 ;
  tr->SetBranchStatus("gen_pt",1);
  tr->SetBranchAddress("gen_pt" , &gen_pt);
  
  std::vector<float>  *gen_eta = 0 ;
  tr->SetBranchStatus("gen_eta",1);
  tr->SetBranchAddress("gen_eta" , &gen_eta);

  std::vector<float>  *gen_phi = 0 ;
  tr->SetBranchStatus("gen_phi",1);
  tr->SetBranchAddress("gen_phi" , &gen_phi);
  
  std::vector<float>  *gen_pdgid = 0 ;
  tr->SetBranchStatus("gen_pdgid",1);
  tr->SetBranchAddress("gen_pdgid" , &gen_pdgid);
  
  std::vector<float>  *gen_energy = 0 ;
  tr->SetBranchStatus("gen_energy",1);
  tr->SetBranchAddress("gen_energy" , &gen_energy);
  
  std::vector<float>  *gen_charge = 0 ;
  tr->SetBranchStatus("gen_charge",1);
  tr->SetBranchAddress("gen_charge" , &gen_charge);
  
  std::vector<float>  *gen_status = 0 ;
  tr->SetBranchStatus("gen_status",1);
  tr->SetBranchAddress("gen_status" , &gen_status);

  gInterpreter->GenerateDictionary("std::vector<std::vector<int> >", "vector");
  std::vector<std::vector<int>>  *gen_daughters = 0 ;
  tr->SetBranchStatus("gen_daughters",1);
  tr->SetBranchAddress("gen_daughters" , &gen_daughters);
  
  std::vector<float>  *genpart_pt = 0 ;
  tr->SetBranchStatus("genpart_pt",1);
  tr->SetBranchAddress("genpart_pt" , &genpart_pt);
  
  std::vector<float>  *genpart_eta = 0 ;
  tr->SetBranchStatus("genpart_eta",1);
  tr->SetBranchAddress("genpart_eta" , &genpart_eta);

  std::vector<float>  *genpart_phi = 0 ;
  tr->SetBranchStatus("genpart_phi",1);
  tr->SetBranchAddress("genpart_phi" , &genpart_phi);
  
  std::vector<float>  *genpart_pid = 0 ;
  tr->SetBranchStatus("genpart_pid",1);
  tr->SetBranchAddress("genpart_pid" , &genpart_pid);

  std::vector<float>  *genpart_mother = 0 ;
  tr->SetBranchStatus("genpart_mother",1);
  tr->SetBranchAddress("genpart_mother" , &genpart_mother);

  std::vector<float>  *genpart_energy = 0 ;
  tr->SetBranchStatus("genpart_energy",1);
  tr->SetBranchAddress("genpart_energy" , &genpart_energy);

  std::vector<float>  *genpart_gen = 0 ;
  tr->SetBranchStatus("genpart_gen",1);
  tr->SetBranchAddress("genpart_gen" , &genpart_gen);
  
  std::vector<float>  *genpart_fromBeamPipe = 0 ;
  tr->SetBranchStatus("genpart_fromBeamPipe",1);
  tr->SetBranchAddress("genpart_fromBeamPipe" , &genpart_fromBeamPipe);

  std::vector<float>  *genpart_exeta = 0 ;
  tr->SetBranchStatus("genpart_exeta",1);
  tr->SetBranchAddress("genpart_exeta" , &genpart_exeta);

  std::vector<float>  *genpart_exphi = 0 ;
  tr->SetBranchStatus("genpart_exphi",1);
  tr->SetBranchAddress("genpart_exphi" , &genpart_exphi);

  Int_t genjet_n = 0 ;
  tr->SetBranchStatus("genjet_n",1);
  tr->SetBranchAddress("genjet_n" , &genjet_n);

  std::vector<float>  *genjet_pt = 0 ;
  tr->SetBranchStatus("genjet_pt",1);
  tr->SetBranchAddress("genjet_pt" , &genjet_pt);

  std::vector<float>  *genjet_eta = 0 ;
  tr->SetBranchStatus("genjet_eta",1);
  tr->SetBranchAddress("genjet_eta" , &genjet_eta);
  
  std::vector<float>  *genjet_phi = 0 ;
  tr->SetBranchStatus("genjet_phi",1);
  tr->SetBranchAddress("genjet_phi" , &genjet_phi);

  std::vector<float>  *genjet_energy = 0 ;
  tr->SetBranchStatus("genjet_energy",1);
  tr->SetBranchAddress("genjet_energy" , &genjet_energy);

  std::vector<float>  *tc_pt = 0 ;
  tr->SetBranchStatus("tc_pt",1);
  tr->SetBranchAddress("tc_pt" , &tc_pt);

  std::vector<float>  *tc_mipPt = 0 ;
  tr->SetBranchStatus("tc_mipPt",1);
  tr->SetBranchAddress("tc_mipPt" , &tc_mipPt);

  std::vector<float>  *tc_eta = 0 ;
  tr->SetBranchStatus("tc_eta",1);
  tr->SetBranchAddress("tc_eta" , &tc_eta);

  std::vector<float>  *tc_phi = 0 ;
  tr->SetBranchStatus("tc_phi",1);
  tr->SetBranchAddress("tc_phi" , &tc_phi);

  std::vector<float>  *tc_layer = 0 ;
  tr->SetBranchStatus("tc_layer",1);
  tr->SetBranchAddress("tc_layer" , &tc_layer);
  
  std::vector<float>  *tc_x = 0 ;
  tr->SetBranchStatus("tc_x",1);
  tr->SetBranchAddress("tc_x" , &tc_x);
  
  std::vector<float>  *tc_y = 0 ;
  tr->SetBranchStatus("tc_y",1);
  tr->SetBranchAddress("tc_y" , &tc_y);
  
  std::vector<float>  *tc_z = 0 ;
  tr->SetBranchStatus("tc_z",1);
  tr->SetBranchAddress("tc_z" , &tc_z);
  
  std::vector<float>  *tc_energy = 0 ;
  tr->SetBranchStatus("tc_energy",1);
  tr->SetBranchAddress("tc_energy" , &tc_energy);
  
  std::vector<float>  *tc_data = 0 ;
  tr->SetBranchStatus("tc_data",1);
  tr->SetBranchAddress("tc_data" , &tc_data);

  std::vector<float>  *tc_subdet = 0 ;
  tr->SetBranchStatus("tc_subdet",1);
  tr->SetBranchAddress("tc_subdet" , &tc_subdet);
  
  std::vector<float>  *tc_wafertype = 0 ;
  tr->SetBranchStatus("tc_wafertype",1);
  tr->SetBranchAddress("tc_wafertype" , &tc_wafertype);

  std::vector<float>  *tc_uncompressedCharge = 0 ;
  tr->SetBranchStatus("tc_uncompressedCharge",1);
  tr->SetBranchAddress("tc_uncompressedCharge" , &tc_uncompressedCharge);

  Int_t cl3d_n = 0 ;
  tr->SetBranchStatus("cl3d_n",1);
  tr->SetBranchAddress("cl3d_n" , &cl3d_n);

  std::vector<float>  *cl3d_pt = 0 ;
  tr->SetBranchStatus("cl3d_pt",1);
  tr->SetBranchAddress("cl3d_pt" , &cl3d_pt);
  
  std::vector<float>  *cl3d_eta = 0 ;
  tr->SetBranchStatus("cl3d_eta",1);
  tr->SetBranchAddress("cl3d_eta" , &cl3d_eta);
  
  std::vector<float>  *cl3d_phi = 0 ;
  tr->SetBranchStatus("cl3d_phi",1);
  tr->SetBranchAddress("cl3d_phi" , &cl3d_phi);
  
  std::cout << "Total number of Events: " << tr->GetEntries() << std::endl;
  
  // nJetEtaBins = 52;
  // Float_t jetEtaBin[53] = {-2.500, -2.322, 
  // 			   -2.172, -2.043, -1.930, -1.830, -1.740, 
  // 			   -1.653, -1.566, -1.479, -1.392, -1.305, 
  // 			   -1.218, -1.131, -1.044, -0.957, -0.870, 
  // 			   -0.783, -0.696, -0.609, -0.522, -0.435, 
  // 			   -0.348, -0.261, -0.174, -0.087, 
  // 			   0.0, 0.087, 0.174, 0.261, 0.348, 
  // 			   0.435, 0.522, 0.609, 0.696, 0.783, 
  // 			   0.870, 0.957, 1.044, 1.131, 1.218, 
  // 			   1.305, 1.392, 1.479, 1.566, 1.653, 
  // 			   1.740, 1.830, 1.930, 2.043, 2.172, 
  // 			   2.322, 2.5};

  // nETBins = 42;
  // Float_t ETBin[43] = {15., 16., 18., 20., 22., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100.,  
  //                      110., 120., 130., 140., 150., 160., 170., 180., 190., 200.,
  // 	               220., 240., 260., 280., 300.,
  // 	               330., 360., 390., 420.,
  //                      440., 480., 520.}; 
  
  // const int nJetEtaBins = 31;
  // Float_t jetEtaBin[32] = { -3.212,
  //                           -3.034, -2.856, -2.678, -2.500, -2.322, 
  // 			    -2.172, -2.043, -1.930, -1.830, -1.740, 
  // 			    -1.653, -1.566, -1.479, -1.392, -1.305, 
  // 			    1.305, 1.392, 1.479, 1.566, 1.653, 
  // 			    1.740, 1.830, 1.930, 2.043, 2.172, 
  // 			    2.322, 2.5, 2.678, 2.856, 3.034,
  // 			    3.212};
  // const int nJetPtBins = 43;
  // Float_t jetPtBin[44] = {0, 15., 16., 18., 20., 22., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100.,  
  //                      110., 120., 130., 140., 150., 160., 170., 180., 190., 200.,
  // 	               220., 240., 260., 280., 300.,
  // 	               330., 360., 390., 420.,
  //                      440., 480., 520.};
  
  //Eta : outermost(Sci) = 1.321, innermost(Si) =  3.152
  //use delta-Eta = 0.0873 (or mutltiple) and make range for |\eta|

  // float deltaEta = 0.0873;
  // const int nJetEtaBins = 21;
  // Float_t jetEtaBin[22];
  // jetEtaBin[0] = 1.321;
  // for(int ieta=1;ieta<=nJetEtaBins;ieta++)
  //   jetEtaBin[ieta] = jetEtaBin[0] + ieta*deltaEta;
  
  const int nJetEtaBins = 6;
  Float_t jetEtaBin[7] = {1.321, 1.7, 2.0, 2.3, 2.6, 2.9, 3.152} ;
  
  const int nJetPtBins = 43;
  Float_t jetPtBin[44] = {0, 15., 16., 18., 20., 22., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100.,  
                       110., 120., 130., 140., 150., 160., 170., 180., 190., 200.,
	               220., 240., 260., 280., 300.,
	               330., 360., 390., 420.,
                       440., 480., 520.};
  
  int nPtResoBins = 400;
  int nPtRatioBins = 1000;
  float etaMinCore = 1.7, etaMaxCore = 2.9;
  //float etaMinCore = 1.6, etaMaxCore = 3.0;
  //////////////// 2D pt,eta histos ////////////////////////////
  TH1F ***hGenClusPtReso, ***hGenTCPtReso, ***hTCClusPtReso;  
  hGenClusPtReso = new TH1F**[nJetEtaBins];
  hGenTCPtReso = new TH1F**[nJetEtaBins];
  hTCClusPtReso = new TH1F**[nJetEtaBins];
  for(int ieta=0;ieta<nJetEtaBins;ieta++){
    hGenClusPtReso[ieta] = new TH1F*[nJetPtBins];
    hGenTCPtReso[ieta] = new TH1F*[nJetPtBins];
    hTCClusPtReso[ieta] = new TH1F*[nJetPtBins];
    for(int ipt=0;ipt<nJetPtBins;ipt++){
      hGenClusPtReso[ieta][ipt] = new TH1F(Form("hGenClusPtReso_%d_%d",ieta,ipt),
					   Form("hGenClusPtReso #eta:(%4.3f-%4.3f) p_{T}:(%2.0f-%2.0f) GeV",jetEtaBin[ieta],jetEtaBin[ieta+1],jetPtBin[ipt],jetPtBin[ipt+1]),
					   nPtResoBins,0.,40.);
      hGenTCPtReso[ieta][ipt] = new TH1F(Form("hGenTCPtReso_%d_%d",ieta,ipt),
					   Form("hGenTCPtReso #eta:(%4.3f-%4.3f) p_{T}:(%2.0f-%2.0f) GeV",jetEtaBin[ieta],jetEtaBin[ieta+1],jetPtBin[ipt],jetPtBin[ipt+1]),
					   nPtResoBins,0.,40.);
      hTCClusPtReso[ieta][ipt] = new TH1F(Form("hTCClusPtReso_%d_%d",ieta,ipt),
					  Form("hTCClusPtReso #eta:(%4.3f-%4.3f) p_{T}:(%2.0f-%2.0f) GeV",jetEtaBin[ieta],jetEtaBin[ieta+1],jetPtBin[ipt],jetPtBin[ipt+1]),
					  nPtResoBins,0.,40.);
    }// jet et
  }//jet eta
  //////////////////////////////////////////////////////////

  //////////////// 2D pt,eta TC ELoss ////////////////////////////
  TH1F ***hTCElossCE_E,***hTCElossCE_H;
  hTCElossCE_E = new TH1F**[nJetEtaBins];
  hTCElossCE_H = new TH1F**[nJetEtaBins];
  for(int ieta=0;ieta<nJetEtaBins;ieta++){
    hTCElossCE_E[ieta] = new TH1F*[nJetPtBins];
    hTCElossCE_H[ieta] = new TH1F*[nJetPtBins];
    for(int ipt=0;ipt<nJetPtBins;ipt++){
      hTCElossCE_E[ieta][ipt] = new TH1F(Form("hTCElossCE_E_%d_%d",ieta,ipt),
					   Form("hTCElossCE_E #eta:(%4.3f-%4.3f) p_{T}:(%2.0f-%2.0f) GeV",jetEtaBin[ieta],jetEtaBin[ieta+1],jetPtBin[ipt],jetPtBin[ipt+1]),
					   nPtResoBins,-400.,400.);
      hTCElossCE_H[ieta][ipt] = new TH1F(Form("hTCElossCE_H_%d_%d",ieta,ipt),
					  Form("hTCElossCE_H #eta:(%4.3f-%4.3f) p_{T}:(%2.0f-%2.0f) GeV",jetEtaBin[ieta],jetEtaBin[ieta+1],jetPtBin[ipt],jetPtBin[ipt+1]),
					  nPtResoBins,-400.,400.);
    }// jet et
  }//jet eta
  //////////////////////////////////////////////////////////
  
  // //////////////// 3D pt,eta TC ELoss ////////////////////////////
  // TH1F ****hTCElossCE_E_3D,****hTCElossCE_H_3D;
  // hTCElossCE_E_3D = new TH1F***[7];
  // hTCElossCE_H_3D = new TH1F***[7];
  // for(int idet=0;idet<7;idet++){
  //   hTCElossCE_E_3D[idet] = new TH1F**[nJetEtaBins];
  //   hTCElossCE_H_3D[idet] = new TH1F**[nJetEtaBins];
  //   for(int ieta=0;ieta<nJetEtaBins;ieta++){
  //     hTCElossCE_E_3D[idet][ieta] = new TH1F*[nJetPtBins];
  //     hTCElossCE_H_3D[idet][ieta] = new TH1F*[nJetPtBins];
  //     for(int ipt=0;ipt<nJetPtBins;ipt++){
  // 	hTCElossCE_E_3D[idet][ieta][ipt] = new TH1F(Form("hTCElossCE_E_3D_%d_%d_%d",idet,ieta,ipt),
  // 						 Form("hTCElossCE_E_3D det:%d, #eta:(%4.3f-%4.3f) p_{T}:(%2.0f-%2.0f) GeV",idet,jetEtaBin[ieta],jetEtaBin[ieta+1],jetPtBin[ipt],jetPtBin[ipt+1]),
  // 						 nPtResoBins,-400.,400.);
  // 	hTCElossCE_H_3D[idet][ieta][ipt] = new TH1F(Form("hTCElossCE_H_3D_%d_%d_%d",idet,ieta,ipt),
  // 						    Form("hTCElossCE_H_3D det:%d, #eta:(%4.3f-%4.3f) p_{T}:(%2.0f-%2.0f) GeV",idet,jetEtaBin[ieta],jetEtaBin[ieta+1],jetPtBin[ipt],jetPtBin[ipt+1]),
  // 						    nPtResoBins,-400.,400.);
  //     }// jet et
  //   }//jet eta
  // }//idet loop
  // //////////////////////////////////////////////////////////

  //////////////// 2D pt,eta TC ELoss ////////////////////////////
  TH1F **hTCElossCEEID, **hTCElossCEHID ;
  hTCElossCEEID = new TH1F*[nJetEtaBins];
  hTCElossCEHID = new TH1F*[nJetEtaBins];
  for(int ieta=0;ieta<nJetEtaBins;ieta++){
    hTCElossCEEID[ieta] = new TH1F(Form("hTCElossCEEID_%d",ieta),
				   Form("hTCElossCEEID #eta:(%4.3f-%4.3f)",jetEtaBin[ieta],jetEtaBin[ieta+1]),
				   nPtResoBins,-400.,400.);
    hTCElossCEHID[ieta] = new TH1F(Form("hTCElossCEHID_%d",ieta),
				   Form("hTCElossCEHID #eta:(%4.3f-%4.3f)",jetEtaBin[ieta],jetEtaBin[ieta+1]),
				   nPtResoBins,-400.,400.);
  }//jet eta
  //////////////////////////////////////////////////////////

  //////////////// 1D pt //////////////////////////////////
  TH1F **hGenClusPtReso1D, **hGenTCPtReso1D, **hTCClusPtReso1D;
  hGenClusPtReso1D = new TH1F*[nJetPtBins];
  hGenTCPtReso1D = new TH1F*[nJetPtBins];
  hTCClusPtReso1D = new TH1F*[nJetPtBins];
  for(int ipt=0;ipt<nJetPtBins;ipt++){
    // hGenClusPtReso1D[ipt] = new TH1F(Form("hGenClusPtReso1D_%d",ipt),
    // 				     Form("hGenClusPtReso1D p_{T}:(%2.0f-%2.0f) GeV",jetPtBin[ipt],jetPtBin[ipt+1]),
    // 				     nPtResoBins,-400.,400.);
    // hGenTCPtReso1D[ipt] = new TH1F(Form("hGenTCPtReso1D_%d",ipt),
    // 				     Form("hGenTCPtReso1D p_{T}:(%2.0f-%2.0f) GeV",jetPtBin[ipt],jetPtBin[ipt+1]),
    // 				     nPtResoBins,-400.,400.);
    // hTCClusPtReso1D[ipt] = new TH1F(Form("hTCClusPtReso1D_%d",ipt),
    // 				    Form("hTCClusPtReso1D p_{T}:(%2.0f-%2.0f) GeV",jetPtBin[ipt],jetPtBin[ipt+1]),
    // 				    nPtResoBins,-400.,400.);
    hGenClusPtReso1D[ipt] = new TH1F(Form("hGenClusPtReso1D_%d",ipt),
				     Form("hGenClusPtReso1D p_{T}:(%2.0f-%2.0f) GeV",jetPtBin[ipt],jetPtBin[ipt+1]),
				     nPtRatioBins,0.,10.);
    hGenTCPtReso1D[ipt] = new TH1F(Form("hGenTCPtReso1D_%d",ipt),
				     Form("hGenTCPtReso1D p_{T}:(%2.0f-%2.0f) GeV",jetPtBin[ipt],jetPtBin[ipt+1]),
				   nPtRatioBins,0.,10.);
    hTCClusPtReso1D[ipt] = new TH1F(Form("hTCClusPtReso1D_%d",ipt),
				    Form("hTCClusPtReso1D p_{T}:(%2.0f-%2.0f) GeV",jetPtBin[ipt],jetPtBin[ipt+1]),
				    nPtRatioBins,0.,10.);
  }// jet pt
  TH2D *h2GenVsClusPt = new TH2D("h2GenVsClusPt","h2GenVsClusPt", 100, 0., 500., 100, 0., 500.);
  TH2D *h2GenVsTCPt = new TH2D("h2GenVsTCPt","h2GenVsTCPt", 100, 0., 500., 100, 0., 500.);
  TH2D *h2TCVsClusPt = new TH2D("h2TCVsClusPt","h2TCVsClusPt", 100, 0., 500., 100, 0., 500.);
  ////////////////////////////////////////////////////////
  
  
  //////////////// 1D pt log x//////////////////////////////////
  double xm1 = -3;
  int nlog_tcptbins = 0;
  std::vector<double> logbins;
  int logmax = 1e6;
  logbins.push_back(0.0);
  for(int itcpt=1 ; itcpt < logmax ; itcpt++){
    double y = 0.1*itcpt;
    double x = log(y);
    if(x-xm1>=0.1){
      //std::cout << "nlog_tcptbins: " << nlog_tcptbins << ", itcpt: " << itcpt << ", y : " << y <<", x : " << x << std::endl;
      xm1 = x;
      logbins.push_back(y);
      nlog_tcptbins++;      
    }
  }
  
  TH1F **hGenjetByTCPtSumVsTCPt1DLog;
  hGenjetByTCPtSumVsTCPt1DLog = new TH1F*[int(logbins.size()-1)];
  for(int ipt=0;ipt<int(logbins.size()-1);ipt++){
    hGenjetByTCPtSumVsTCPt1DLog[ipt] = new TH1F(Form("hGenjetByTCPtSumVsTCPt1DLog_%d",ipt),
						Form("hGenjetByTCPtSumVsTCPt1DLog p_{T}:(%2.1f-%2.1f) MeV",logbins.at(ipt),logbins.at(ipt+1)),
						nPtResoBins,0.,40.);
  }
  //////////////// 1D pt //////////////////////////////////

  //////////////// 2D pt //////////////////////////////////
  const int nTCPtBins = 100;
  TH1F ***hGenjetByTCPtSumVsTCPt2DDetLog;
  hGenjetByTCPtSumVsTCPt2DDetLog = new TH1F**[7];
  for(int idet=0;idet<7;idet++){
    hGenjetByTCPtSumVsTCPt2DDetLog[idet] = new TH1F*[int(logbins.size()-1)];
    for(int ipt=0;ipt<int(logbins.size()-1);ipt++){
      hGenjetByTCPtSumVsTCPt2DDetLog[idet][ipt] = new TH1F(Form("hGenjetByTCPtSumVsTCPt2DDetLog_%d_%d",idet,ipt),
						 Form("hGenjetByTCPtSumVsTCPt2DDetLog idet: %d, p_{T}:(%2.0f-%2.0f) GeV",idet,logbins.at(ipt),logbins.at(ipt+1)),
						 nPtResoBins,0.,40.);
    }
  }
  TH1F ***hGenjetByTCPtSumVsTCPt2DLog;
  hGenjetByTCPtSumVsTCPt2DLog = new TH1F**[nJetEtaBins];
  for(int ieta=0;ieta<nJetEtaBins;ieta++){
    hGenjetByTCPtSumVsTCPt2DLog[ieta] = new TH1F*[int(logbins.size()-1)];
    for(int ipt=0;ipt<int(logbins.size()-1);ipt++){
      hGenjetByTCPtSumVsTCPt2DLog[ieta][ipt] = new TH1F(Form("hGenjetByTCPtSumVsTCPt2DLog_%d_%d",ieta,ipt),
							Form("hGenjetByTCPtSumVsTCPt2DLog #eta:(%4.3f-%4.3f), p_{T}:(%2.0f-%2.0f) GeV",jetEtaBin[ieta],jetEtaBin[ieta+1],logbins.at(ipt),logbins.at(ipt+1)),
							nPtResoBins,0.,40.);
    }
  }
  //////////////// 1D pt //////////////////////////////////
  
  // /////////////////// 1D pt CEH only correction /////////////////
  // TH1F **hHadGenClusPtReso1D, **hHadGenTCPtReso1D, **hHadTCClusPtReso1D;
  // hHadGenClusPtReso1D = new TH1F*[nJetPtBins];
  // hHadGenTCPtReso1D = new TH1F*[nJetPtBins];
  // hHadTCClusPtReso1D = new TH1F*[nJetPtBins];
  // for(int ipt=0;ipt<nJetPtBins;ipt++){
  //   hHadGenClusPtReso1D[ipt] = new TH1F(Form("hHadGenClusPtReso1D_%d",ipt),
  // 				     Form("hHadGenClusPtReso1D p_{T}:(%2.0f-%2.0f) GeV",jetPtBin[ipt],jetPtBin[ipt+1]),
  // 				     nPtResoBins,-400.,400.);
  //   hHadGenTCPtReso1D[ipt] = new TH1F(Form("hHadGenTCPtReso1D_%d",ipt),
  // 				     Form("hHadGenTCPtReso1D p_{T}:(%2.0f-%2.0f) GeV",jetPtBin[ipt],jetPtBin[ipt+1]),
  // 				     nPtResoBins,-400.,400.);
  //   hHadTCClusPtReso1D[ipt] = new TH1F(Form("hHadTCClusPtReso1D_%d",ipt),
  // 				    Form("hHadTCClusPtReso1D p_{T}:(%2.0f-%2.0f) GeV",jetPtBin[ipt],jetPtBin[ipt+1]),
  // 				    nPtResoBins,-400.,400.);
  // }// jet pt
  // ////////////////////////////////////////////////////////

  // /////////////////// 2D CEH only correction /////////////////
  // TH1F ***hHadGenClusPtReso2D;
  // hHadGenClusPtReso2D = new TH1F**[nJetEtaBins];
  // for(int ieta=0;ieta<nJetEtaBins;ieta++){
  //   hHadGenClusPtReso2D[ieta] = new TH1F*[nJetPtBins];
  //   for(int ipt=0;ipt<nJetPtBins;ipt++){
  //     hHadGenClusPtReso2D[ieta][ipt] = new TH1F(Form("hHadGenClusPtReso2D_%d_%d",ieta,ipt),
  // 					   Form("hHadGenClusPtReso2D #eta:(%4.3f-%4.3f) p_{T}:(%2.0f-%2.0f) GeV",jetEtaBin[ieta],jetEtaBin[ieta+1],jetPtBin[ipt],jetPtBin[ipt+1]),
  // 					   nPtResoBins,-400.,400.);
  //   }// jet et
  // }//jet eta
  // ////////////////////////////////////////////////////////
  
  // /////////////////// 3D CEH only correction /////////////////
  // TH1F ****hHadGenClusPtReso3D;
  // hHadGenClusPtReso3D = new TH1F***[7];
  // for(int idet=0;idet<7;idet++){
  //   hHadGenClusPtReso3D[idet] = new TH1F**[nJetEtaBins];
  //   for(int ieta=0;ieta<nJetEtaBins;ieta++){
  //     hHadGenClusPtReso3D[idet][ieta] = new TH1F*[nJetPtBins];
  //     for(int ipt=0;ipt<nJetPtBins;ipt++){
  // 	hHadGenClusPtReso3D[idet][ieta][ipt] = new TH1F(Form("hHadGenClusPtReso3D_%d_%d_%d",idet,ieta,ipt),
  // 						 Form("hHadGenClusPtReso3D det:%d, #eta:(%4.3f-%4.3f) p_{T}:(%2.0f-%2.0f) GeV",idet,jetEtaBin[ieta],jetEtaBin[ieta+1],jetPtBin[ipt],jetPtBin[ipt+1]),
  // 						 nPtResoBins,-400.,400.);
  //     }// jet et
  //   }//jet eta
  // }//idet loop
  // ////////////////////////////////////////////////////////
  
  //////////////// 3D pt,eta histos ////////////////////////////
  // TH1F ****hGenjetByTCPtSumVsTCPt3DLog;
  // hGenjetByTCPtSumVsTCPt3DLog = new TH1F***[7];
  // for(int idet=0;idet<7;idet++){
  //   hGenjetByTCPtSumVsTCPt3DLog[idet] = new TH1F**[nJetEtaBins];
  //   for(int ieta=0;ieta<nJetEtaBins;ieta++){
  //     hGenjetByTCPtSumVsTCPt3DLog[idet][ieta] = new TH1F*[int(logbins.size()-1)];
  //     for(int ipt=0;ipt<int(logbins.size()-1);ipt++){
  // 	hGenjetByTCPtSumVsTCPt3DLog[idet][ieta][ipt] = new TH1F(Form("hGenjetByTCPtSumVsTCPt3DLog_%d_%d_%d",idet,ieta,ipt),
  // 								Form("hGenjetByTCPtSumVsTCPt3DLog idet: %d, #eta:(%4.3f-%4.3f) , p_{T}:(%2.0f-%2.0f) GeV",idet,jetEtaBin[ieta],jetEtaBin[ieta+1],logbins.at(ipt),logbins.at(ipt+1)),
  // 								nPtResoBins,0.,40.);
  //     }
  //   }
  // }
  
  TH1F ****hGenClusPtReso3D; 
  hGenClusPtReso3D = new TH1F***[7];
  for(int idet=0;idet<7;idet++){
    hGenClusPtReso3D[idet] = new TH1F**[nJetEtaBins];
    for(int ieta=0;ieta<nJetEtaBins;ieta++){
      hGenClusPtReso3D[idet][ieta] = new TH1F*[nJetPtBins];
      for(int ipt=0;ipt<nJetPtBins;ipt++){
	hGenClusPtReso3D[idet][ieta][ipt] = new TH1F(Form("hGenClusPtReso3D_%d_%d_%d",idet,ieta,ipt),
						 Form("hGenClusPtReso3D det:%d, #eta:(%4.3f-%4.3f) p_{T}:(%2.0f-%2.0f) GeV",idet,jetEtaBin[ieta],jetEtaBin[ieta+1],jetPtBin[ipt],jetPtBin[ipt+1]),
						 nPtResoBins,0.,40.);
      }// jet et
    }//jet eta
  }//idet loop
  //////////////////////////////////////////////////////////
  
  ///////////Energy correction for Cl3d cluster //////////
  TH1F **hGenClusPtReso1Dcl3d, **hNewOldPtReso1Dcl3d;
  hGenClusPtReso1Dcl3d = new TH1F*[nJetPtBins];
  hNewOldPtReso1Dcl3d = new TH1F*[nJetPtBins];
  for(int ipt=0;ipt<nJetPtBins;ipt++){
    hGenClusPtReso1Dcl3d[ipt] = new TH1F(Form("hGenClusPtReso1Dcl3d_%d",ipt),
				     Form("hGenClusPtReso1Dcl3d p_{T}:(%2.0f-%2.0f) GeV",jetPtBin[ipt],jetPtBin[ipt+1]),
				     nPtRatioBins,0.,10.);
    hNewOldPtReso1Dcl3d[ipt] = new TH1F(Form("hNewOldPtReso1Dcl3d_%d",ipt),
				     Form("hNewOldPtReso1Dcl3d p_{T}:(%2.0f-%2.0f) GeV",jetPtBin[ipt],jetPtBin[ipt+1]),
				     nPtRatioBins,0.,10.);
  }// jet pt
  TH2D *h2GenVsClusPtcl3d = new TH2D("h2GenVsClusPtcl3d","h2GenVsClusPtcl3d", 100, 0., 500., 100, 0., 500.);
  TH2D *h2NewVsOldPtcl3d = new TH2D("h2NewVsOldPtcl3d","h2NewVsOldPtcl3d", 100, 0., 500., 100, 0., 500.);
  TProfile *pNewVsOldPtcl3d = new TProfile("pNewVsOldPtcl3d","pNewVsOldPtcl3d", 100, 0., 500., 0., 500.);
  ////////////////////////////////////////////////////////  
  TH2F *hJetEtaPtBin = new TH2F("hJetEtaPtBin","hJetEtaPtBin", nJetEtaBins, jetEtaBin, nJetPtBins, jetPtBin);
  TH2F *hJetEtaPtBinDet = new TH2F("hJetEtaPtBinDet","hJetEtaPtBinDEt", nJetEtaBins, jetEtaBin, nJetPtBins, jetPtBin);
  
  TH2F *hJetEtaTCPtBin = new TH2F("hJetEtaTCPtBin","hJetEtaTCPtBin", nJetEtaBins, jetEtaBin[0], jetEtaBin[nJetEtaBins-1], 100, 0.,100.);
  const Double_t *plogbins = logbins.data();
  TH2F *hJetEtaTCPtBinLog = new TH2F("hJetEtaTCPtBinLog","hJetEtaTCPtBinLog", nJetEtaBins, jetEtaBin[0], jetEtaBin[nJetEtaBins-1], int(logbins.size()-1), plogbins);
  
  ////////////// TC-vs-genpt /////////////////////////
  TH2D *h2TCvsGenPt = new TH2D("h2TCvsGenPt","h2TCvsGenPt", 250, 0., 500., 250, 0., 500.);
  TH2D *h2TCvsGenPt_CoreEta = new TH2D("h2TCvsGenPt_CoreEta","h2TCvsGenPt_CoreEta", 250, 0., 500., 250, 0., 500.);
  TH2D *h2TCvsGenPt_CoreEta_rozcut0p15 = new TH2D("h2TCvsGenPt_CoreEta_rozcut0p15","h2TCvsGenPt_CoreEta", 250, 0., 500., 250, 0., 500.);
  TH2D *h2TCvsGenPt_CoreEta_mpptcut2 = new TH2D("h2TCvsGenPt_CoreEta_mpptcut2","h2TCvsGenPt_CoreEta", 250, 0., 500., 250, 0., 500.);
  
  TH2D *h2TCvsGenPt_CEE = new TH2D("h2TCvsGenPt_CEE","h2TCvsGenPt_CEE", 250, 0., 500., 250, 0., 500.);
  TH2D *h2TCvsGenPt_CEE_CoreEta = new TH2D("h2TCvsGenPt_CEE_CoreEta","h2TCvsGenPt_CEE_CoreEta", 250, 0., 500., 250, 0., 500.);
  
  TProfile *pTCvsGenPt = new TProfile("pTCvsGenPt","pTCvsGenPt", 250, 0., 500., 0., 500.);
  TProfile *pTCvsGenPt_CoreEta = new TProfile("pTCvsGenPt_CoreEta","pTCvsGenPt_CoreEta", 250, 0., 500., 0., 500.);
  TProfile *pTCvsGenPt_CoreEta_rozcut0p15 = new TProfile("pTCvsGenPt_CoreEta_rozcut0p15","pTCvsGenPt_CoreEta_rozcut0p15", 250, 0., 500., 0., 500.);
  TProfile *pTCvsGenPt_CoreEta_mpptcut2 = new TProfile("pTCvsGenPt_CoreEta__mpptcut2","pTCvsGenPt_CoreEta__mpptcut2", 250, 0., 500., 0., 500.);
  TProfile *pTCvsGenPt_CEE = new TProfile("pTCvsGenPt_CEE","pTCvsGenPt_CEE", 250, 0., 500., 0., 500.);
  TProfile *pTCvsGenPt_CEE_CoreEta = new TProfile("pTCvsGenPt_CEE_CoreEta","pTCvsGenPt_CEE_CoreEta", 250, 0., 500., 0., 500.);
  
  ////////////// genpt/totTC-vs-(pt/eta) /////////////////////////
  TH2D *h2GenjetTCPtsumvsGenPt = new TH2D("h2GenjetTCPtsumvsGenPt","p_{T}^{genJet}/#Sigmap_{T}^{TC} vs p_{T}^{genJet}", 100, 0., 500., 200, 0., 10.);
  TH2D *h2GenjetTCPtsumvsGenEta = new TH2D("h2GenjetTCPtsumvsGenEta","p_{T}^{genJet}/#Sigmap_{T}^{TC} vs #eta^{genJet}",  90, 0., 3.93, 200, 0., 10.);
  TProfile *pGenjetTCPtsumvsGenPt = new TProfile("pGenjetTCPtsumvsGenPt","p_{T}^{genJet}/#Sigmap_{T}^{TC} vs p_{T}^{genJet}", 100, 0., 500., 0., 10.);
  TProfile *pGenjetTCPtsumvsGenEta = new TProfile("pGenjetTCPtsumvsGenEta","p_{T}^{genJet}/#Sigmap_{T}^{TC} vs #eta^{genJet}", 90, 0., 3.93, 0., 10.);

  TH2D *h2GenjetTCPtsumvsGenPt_CoreEta = new TH2D("h2GenjetTCPtsumvsGenPt_CoreEta","p_{T}^{genJet}/#Sigmap_{T}^{TC} vs p_{T}^{genJet}", 100, 0., 500., 200, 0., 10.);
  TProfile *pGenjetTCPtsumvsGenPt_CoreEta = new TProfile("pGenjetTCPtsumvsGenPt_CoreEta","p_{T}^{genJet}/#Sigmap_{T}^{TC} vs p_{T}^{genJet}", 100, 0., 500., 0., 10.);

  TH2D *h2GenjetTCPtsumvsGenPt_CEE = new TH2D("h2GenjetTCPtsumvsGenPt_CEE","p_{T}^{genJet}/#Sigmap_{T}^{TC} vs p_{T}^{genJet}", 100, 0., 500., 200, 0., 10.);
  TH2D *h2GenjetTCPtsumvsGenEta_CEE = new TH2D("h2GenjetTCPtsumvsGenEta_CEE","p_{T}^{genJet}/#Sigmap_{T}^{TC} vs #eta^{genJet}", 90, 0., 3.93, 200, 0., 10.);
  TH2D *h2GenjetTCPtsumvsGenPt_CEE_CoreEta = new TH2D("h2GenjetTCPtsumvsGenPt_CEE_CoreEta","p_{T}^{genJet}/#Sigmap_{T}^{TC} vs p_{T}^{genJet}", 100, 0., 500., 200, 0., 10.);
  TH2D *h2GenjetTCPtsumvsTCPtsum_CEE_CoreEta = new TH2D("h2GenjetTCPtsumvsTCPtsum_CEE_CoreEta","p_{T}^{genJet}/#Sigmap_{T}^{TC} vs #Sigmap_{T}^{TC}", 100, 0., 500., 200, 0., 10.);
  
  TProfile *pGenjetTCPtsumvsGenPt_CEE = new TProfile("pGenjetTCPtsumvsGenPt_CEE","p_{T}^{genJet}/#Sigmap_{T}^{TC} vs p_{T}^{genJet}", 100, 0., 500., 0., 10.);
  TProfile *pGenjetTCPtsumvsGenEta_CEE = new TProfile("pGenjetTCPtsumvsGenEta_CEE","p_{T}^{genJet}/#Sigmap_{T}^{TC} vs #eta^{genJet}", 90, 0., 3.93, 0., 10.);
  TProfile *pGenjetTCPtsumvsGenPt_CEE_CoreEta = new TProfile("pGenjetTCPtsumvsGenPt_CEE_CoreEta","p_{T}^{genJet}/#Sigmap_{T}^{TC} vs p_{T}^{genJet}", 100, 0., 500., 0., 10.);
  TProfile *pGenjetTCPtsumvsTCPtsum_CEE_CoreEta = new TProfile("pGenjetTCPtsumvsTCPtsum_CEE_CoreEta","p_{T}^{genJet}/#Sigmap_{T}^{TC} vs #Sigmap_{T}^{TC}", 100, 0., 500., 0., 10.);

  /////////////////// clus-vs-tc ////////////////////////////////
  TH2D *h2ClusVsTCPt = new TH2D("h2ClusVsTCPt","h2ClusVsTCPt", 250, 0., 500., 250, 0., 500.);
  TH2D *h2ClusVsTCPt_CoreEta = new TH2D("h2ClusVsTCPt_CoreEta","h2ClusVsTCPt_CoreEta", 250, 0., 500., 250, 0., 500.);
  TH2D *h2ClusVsTCPt_CoreEta_rozcut0p15 = new TH2D("h2ClusVsTCPt_CoreEta_rozcut0p15","h2ClusVsTCPt_CoreEta", 250, 0., 500., 250, 0., 500.);
  TH2D *h2TotClusVsTCPt = new TH2D("h2TotClusVsTCPt","h2TotClusVsTCPt", 250, 0., 500., 250, 0., 500.);
  TH2D *h2TotClusVsTCPt_CoreEta = new TH2D("h2TotClusVsTCPt_CoreEta","h2TotClusVsTCPt_CoreEta", 250, 0., 500., 250, 0., 500.);
  TH2D *h2TotClusVsTCPt_CoreEta_rozcut0p15 = new TH2D("h2TotClusVsTCPt_CoreEta_rozcut0p15","h2TotClusVsTCPt_CoreEta", 250, 0., 500., 250, 0., 500.);

  ///////////////////////////////////////////////////////////
  TH2D *h2TCPtvsGenPtByTCPtsum = new TH2D("h2TCPtvsGenPtByTCPtsum","h2TCPtvsGenPtByTCPtsum", 6000, 0., 60000., 200, -4., 4.); //TCpt in MeV
  TProfile *pTCPtvsGenPtByTCPtsum = new TProfile("pTCPtvsGenPtByTCPtsum","pTCPtvsGenPtByTCPtsum", 6000, 0., 60000., -4., 4.);
  TProfile *pTCPtvsGenPtByTCPtsumLogx = new TProfile("pTCPtvsGenPtByTCPtsumLogx","pTCPtvsGenPtByTCPtsumLogx", int(logbins.size()-1), plogbins, -4., 4.);
  TProfile *pTCPtvsGenPtByTCPtsumDet[7];
  for(int idet=0;idet<7;idet++){
    pTCPtvsGenPtByTCPtsumDet[idet] = new TProfile(Form("pTCPtvsGenPtByTCPtsumDet_%d",idet),Form("pTCPtvsGenPtByTCPtsumDet_%d",idet), 6000, 0., 60000., -4., 4.);
  }
  ////////////////////Unscaled Pt,eta, phi of TCs//////////////////////////////
  TH1D *hTCPt = new TH1D("hTCPt","p_{T} distribution of TCs",500,0.0,500.0);
  hTCPt->GetXaxis()->SetTitle("p_{T} (GeV)");
  hTCPt->GetYaxis()->SetTitle("Entries");
  TH1D *hTCEta = new TH1D("hTCEta","#eta distribution of TCs",1000,-4.15,4.15);
  hTCEta->GetXaxis()->SetTitle("#eta");
  hTCEta->GetYaxis()->SetTitle("Entries");
  TH1D *hTCPhiPos = new TH1D("hTCPhiPos","#phi distribution of TCs for #eta>0",400,-200.,200.);
  hTCPhiPos->GetXaxis()->SetTitle("#phi (degree)");
  hTCPhiPos->GetYaxis()->SetTitle("Entries");
  hTCPhiPos->SetLineColor(kBlue);
  TH1D *hTCPhiNeg = new TH1D("hTCPhiNeg","#phi distribution of TCs for #eta<0",400,-200.,200.);
  hTCPhiNeg->GetXaxis()->SetTitle("#phi (degree)");
  hTCPhiNeg->GetYaxis()->SetTitle("Entries");
  hTCPhiNeg->SetLineColor(kBlue);
  TH2D *h2TCEtaPhi = new TH2D("h2TCEtaPhi","#eta-#phi distribution of TCs", 400,-200.,200., 1000,-4.15,4.15);
  h2TCEtaPhi->GetXaxis()->SetTitle("#phi");
  h2TCEtaPhi->GetYaxis()->SetTitle("#eta");
  //////////////////////////////////////////////////////////
  TH1D *tcRoz = new TH1D("tcRoz",Form("Roz distribution of TCs w.r.t gen"),1000,0.0,1.0);
  TH1D *tcRozCore = new TH1D("tcRozCore",Form("Roz distribution of TCs w.r.t gen"),1000,0.0,1.0);
  TH1D *htcMipPt = new TH1D("htcMipPt","tc mipPt",2500,0.0,2500.0);
  TH1D *tcRozDet[7];
  for(int idet=0;idet<7;idet++){
    tcRozDet[idet] = new TH1D(Form("tcRozDet_%d",idet),Form("tcRozDet_%d",idet), 1000,0.0,1.0);
  }
  TH1D *tcRozPtwt = new TH1D("tcRozPtwt",Form("Roz distribution of TCs w.r.t gen"),1000,0.0,1.0);
  TH1D *clusRoz = new TH1D("clusRoz",Form("Roz distribution of Clusters w.r.t gen"),1000,0.0,1.0);
  TH1D *clusRozPtwt = new TH1D("clusRozPtwt",Form("Roz distribution of Clusters w.r.t gen"),1000,0.0,1.0);
  double tcRoZCut = 0.15;
  ///////////////////////////////////////////////////////////
  TH1D *hGenvsJetPtDiff = new TH1D("hGenvsJetPtDiff",Form("(p_{T}^{gen}-p_{T}^{genJet}) GeV"),200,-100.0,100.0);
  TProfile *pGenvsJetPtDiff = new TProfile("pGenvsJetPtDiff","p_{T}^{gen}/p_{T}^{genJet}", 100, 0., 500., 0., 10.);
  TH1D *hGenvsJetDR = new TH1D("hGenvsJetDR",Form("#DeltaR"),800,-0.1,0.7);
  TProfile *pPtDiffvsDR = new TProfile("pPtDiffvsDR",Form("(p_{T}^{gen}/p_{T}^{genJet}) vs #DeltaR"),800,-0.1,0.7,0.,10.);
  TH2D *h2PtDiffvsDR = new TH2D("h2PtDiffvsDR",Form("(p_{T}^{gen}/p_{T}^{genJet}) vs #DeltaR"),800,-0.1,0.7,100, 0.,10.);

  ////////////////////////////////////////////////////////////
  TH1D *hNJets = new TH1D("hNJets",Form("Nof VBF Jets"), 10,-0.5,9.5);
  TH1D *hNJetsF1 = new TH1D("hNJetsF1",Form("Nof VBF Jets (filtered for atmost 2)"), 10,-0.5,9.5);
  ////////////////////////////////////////////////////////////

  ///////////// Nof Clusters ////////////////////////////////
  TH1D *hNClus1GeV = new TH1D("hNClus1GeV","Nof cluster with p_{T}>1 GeV per side", 50,-0.5,49.5);
  TH1D *hNClus3GeV = new TH1D("hNClus3GeV","Nof cluster with p_{T}>3 GeV per side", 50,-0.5,49.5);
  TH1D *hNClus5GeV = new TH1D("hNClus5GeV","Nof cluster with p_{T}>5 GeV per side ", 50,-0.5,49.5);
  TH1D *hNClus10GeV = new TH1D("hNClus10GeV","Nof cluster with p_{T}>10 GeV per side", 50,-0.5,49.5);
  //////////////////////////////////////////////////////////
  
  ///////////// TC Correction ///////////////////////////////////
  TF1* f2 = new TF1("func","[0] - [1]/pow(x,1/3) - [2]/pow(x,4/3) + [3]/pow(x,2) - [4]*x ",3.,3.5e4);
  f2->SetParameters(1.58371, 0.406349, 0.374516, 0.636789, 2.9813e-06); //PU0 VBF
  
  // //f2->SetParameters(1.07448, 0.915572, 0.320687, 0.657203, -2.58768e-05); //PU200 VBF CB
  // TF1 *fhigh = new TF1("fhigh","[0] + [1]*x - [2]/x",0.4e3,12.e3);
  // fhigh->SetParameters(1.389, 7.78113e-05, 2,137);
  // TF1 *fhigh1 = new TF1("fhigh1","[0] * log (x-[1])",0.2e3,12.e3);
  // fhigh1->SetParameters(0.166198, -11.3154);
  
  // TF1* f2 = new TF1(Form("func"),"[0] - [1]/pow(x,1/3) - [2]/pow(x,4/3) + [3]/pow(x,2) - [4]*log (x-[5]) + [6]*x - [7]/x",0.1,7.5e4);
  // f2->SetParameters(0.816674, 1.17338, -136.175, 0.498167, -0.0698533, -1362.15, 9.2412e-06, 136.382); //PU200 VBF CB with userange 0 to 0.8  
  // TF1 *fhigh = new TF1("fhigh","[0] * log (x-[1]) + [2] + [3]*x",0.5e3,12.e3);
  // fhigh->SetParameters(0.0538316, 459.772, 0.935318, 7.96371e-05);
  
  ////////////////Cluster correction ///////////////////////////
  TF1* fClusCorr = new TF1("fClusCorr","[0]+[1]/(x+[2])",7,500);
  fClusCorr->SetParameters(1.22466, 22.8516, 0.401581); //PU0 VBF
  
  TF1* fcl3d = new TF1("fcl3d","[0]+[1]/(x+[2])",7,500);
  fcl3d->SetParameters(1.17576, 19.7599, 0.908577); //PU0 VBF
  
  //////////////////////////////////////////////////////////
  uint64_t totalEntries = tr->GetEntries();
  if(nofEvents==0 or nofEvents>totalEntries) nofEvents = totalEntries;
  
  std::cout << "Total Entries : " << totalEntries << std::endl;
  std::cout << "Loop for: " << nofEvents << " entries"  << std::endl;
  
  const auto default_precision{std::cout.precision()};
  
  // std::string board_config = "../hgcal-tpg-fe-data/local-input/stage2/firmware-data/vbf_Captures/CaptureStage2_250404_1118/Stage2Configuration.yaml" ;
  // TPGStage2Configuration::Stage2Board sb;
  // sb.readConfigYaml(board_config.c_str());
  // sb.print();
  
  TPGStage2Configuration::ClusPropLUT cplut;
  cplut.readMuEtaLUT("input/stage2/configuration/mean_eta_LUT.csv");
  cplut.readSigmaEtaLUT("input/stage2/configuration/sigma_eta_LUT.csv");
  
  TPGStage2Emulation::Stage2::_rOverZ = sidelength * sqrt(3.0) ;
  TPGStage2Emulation::Stage2 s2Clustering;
  s2Clustering.setClusPropLUT(&cplut);
  //s2Clustering.setConfiguration(&sb);
  s2Clustering.setkpower(0);
  std::cout << " Stage2::_rOverZ : " << s2Clustering.getROverZ() << std::endl;
  
  //TPGTriggerCellFloats tcf0,tcf1;
  TPGTCFloats tcf0,tcf1;
  for (Long64_t ievent = 0 ; ievent < nofEvents ; ievent++ ) {
    
    tr->GetEntry(ievent) ;
    
    if(doPrint) std::cout<<"Event : "<< ievent <<", nof gen "<< gen_n << ", nof_gen: " << genpart_pt->size()  <<", TCs : "<< tc_pt->size() << std::endl;
    if(ievent%1==0) std::cout<<"Event : "<< ievent <<", nof TCs : "<< tc_pt->size() << std::endl;
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::vector<int> taudlist,taugdlist;
    std::vector<JetPart> genlist;
    for(int igen=0; igen<gen_n; igen++ ){
      TParticlePDG *partPDG = TDatabasePDG::Instance()->GetParticle(gen_pdgid->at(igen));
      std::string gennom = (!partPDG)?"unknown":partPDG->GetName() ;
      std::string daughterlist = "(";
      for(int const& idx: gen_daughters->at(igen)) {
  	daughterlist += std::to_string(idx);
  	daughterlist += ",";
  	if(gennom.find("       pi+")!=std::string::npos) taudlist.push_back(idx);
  	if(taudlist.size()>0 and std::find(taudlist.begin(),taudlist.end(),igen)!=taudlist.end()) taugdlist.push_back(idx);
      }
      daughterlist += ")";
      
      JetPart p;
      p.name = (!partPDG)?"unknown":partPDG->GetName();
      p.index = igen;
      genlist.push_back(p);
      
      if(doPrint and igen<=10
	 // and
  	 // (gennom.find("h0")!=std::string::npos
  	 // or
  	 // gennom.find("tau")!=std::string::npos
  	 // or
  	 // (taudlist.size()>0 and std::find(taudlist.begin(),taudlist.end(),igen)!=taudlist.end())
  	 // or
  	 // (taugdlist.size()>0 and std::find(taugdlist.begin(),taugdlist.end(),igen)!=taugdlist.end()))
  	 ){
  	std::cout << "gen-ievent: " << std::setprecision(default_precision) << std::setw(4) <<ievent << ", igen: " << std::setw(4) << igen
  		  <<", pid: " << std::setw(5) << gen_pdgid->at(igen) << ", status: " << std::setw(5) << gen_status->at(igen)
  		  <<", Name: " << std::setw(10) << gennom.data() //((!partPDG)?"unknown":partPDG->GetName())
  		  <<" gen:(pt,eta,phi,energy) : (" << std::fixed << std::setprecision(2) << std::setw(10) << gen_pt->at(igen)
  		  << ", " << std::setw(10) << gen_eta->at(igen) << ", " << std::setw(10) << (TMath::RadToDeg()*gen_phi->at(igen)) << ", " << std::setw(10) << gen_energy->at(igen) << ") "
  		  <<" daughters: " << (gen_daughters->at(igen)).size()
  		  << ", " << daughterlist
	  //<< ", taudlist.size(): " << taudlist.size() 
  		  << std::defaultfloat
  		  << std::endl;
      }//print condition
    }//gen loop
    
    std::vector<JetPart> partlist;
    for(int ipart=0; ipart<genpart_pt->size(); ipart++ ){
      TParticlePDG *partPDG = TDatabasePDG::Instance()->GetParticle(genpart_pid->at(ipart));
      // if(
      // 	 (taudlist.size()>0 )//and std::find(taudlist.begin(),taudlist.end(),(genpart_gen->at(ipart)-1))!=taudlist.end())
      // 	 or
      // 	 (taugdlist.size()>0)//and std::find(taugdlist.begin(),taugdlist.end(),(genpart_gen->at(ipart)-1))!=taugdlist.end())
      // 	 ){
      
	JetPart p;
	p.name = (!partPDG)?"unknown":partPDG->GetName();
	p.index = ipart;
	partlist.push_back(p);
	// if(doPrint)
	//   std::cout << "genpart-ievent: " << std::setprecision(default_precision) << std::setw(4) <<ievent << ", ipart: " << std::setw(4) << ipart
	// 	    <<", pid: " << std::setw(5) << genpart_pid->at(ipart) << ", mother: " << std::setw(4) << genpart_mother->at(ipart)
	//     //<<", status: " << std::setw(5) << gen_status->at(ipart)
	// 	    <<", gen: " << std::setw(4) << (genpart_gen->at(ipart)-1) //<< ", fromBeamPipe: " << std::setw(4) << genpart_fromBeamPipe->at(ipart)
	// 	    <<", Name: " << std::setw(10) << ((!partPDG)?"unknown":partPDG->GetName())
	// 	    <<" part:(pt,eta,phi,energy) : (" << std::fixed << std::setprecision(2) << std::setw(8) << genpart_pt->at(ipart)
	// 	    << ", " << std::setw(8) << genpart_eta->at(ipart) << ", " //<< std::setw(5) << genpart_exeta->at(ipart)
	// 	    << ", " << std::setw(8) << (TMath::RadToDeg()*genpart_phi->at(ipart)) //<< ", " << std::setw(5) << (TMath::RadToDeg()*genpart_exphi->at(ipart))
	// 	    << ", " << std::setw(8) << genpart_energy->at(ipart) << ") "
	// 	    << std::defaultfloat
	// 	    << std::endl;
      // }
    }//genpart loop
    
    std::vector<JetPart> jetlist;
    double delRTh = 0.25;
    int nofHGCalJetsPos = 0, nofHGCalJetsNeg = 0;
    int nofJetsLEPhiPos = 0, nofJetsLEPhiNeg = 0;
    for(int ijet=0; ijet<genjet_n; ijet++ ){
      if(fabs(genjet_eta->at(ijet))<=1.321 or fabs(genjet_eta->at(ijet))>=3.152) continue;
      double phiDeg = TMath::RadToDeg()*genjet_phi->at(ijet);
      if(genjet_eta->at(ijet)>0.){
	nofHGCalJetsPos++;
	if((phiDeg>=-95. and phiDeg<=-85.) or (phiDeg>=25. and phiDeg<=35.) or (phiDeg>=145. and phiDeg<=155.)) nofJetsLEPhiPos++;
      }else{
	nofHGCalJetsNeg++;
	if((phiDeg>=-95. and phiDeg<=-85.) or (phiDeg>=25. and phiDeg<=35.) or (phiDeg>=145. and phiDeg<=155.)) nofJetsLEPhiNeg++;
      }
      JetPart jet;
      // jet.name = "";
      // jet.index = ijet ;
      // jetlist.push_back(jet);
      int minDeltaRindex = -1;
      float minDeltaR = 1.0; 
      float minDR[3] = {1., 1., 1.};      
      // for(int ipart=0; ipart<partlist.size(); ipart++){
      // 	int refpart = partlist.at(ipart).index;
      // double deltaR = TMath::Sqrt((genpart_eta->at(refpart)-genjet_eta->at(ijet))*(genpart_eta->at(refpart)-genjet_eta->at(ijet)) + (genpart_phi->at(refpart)-genjet_phi->at(ijet))*(genpart_phi->at(refpart)-genjet_phi->at(ijet)));
      for(int ipart=0; ipart<genlist.size(); ipart++){
	int refpart = genlist.at(ipart).index;
	double deltaR = TMath::Sqrt((gen_eta->at(refpart)-genjet_eta->at(ijet))*(gen_eta->at(refpart)-genjet_eta->at(ijet)) + (gen_phi->at(refpart)-genjet_phi->at(ijet))*(gen_phi->at(refpart)-genjet_phi->at(ijet)));
	bool jetcondn = (deltaR<delRTh and deltaR<minDeltaR);
	double deltaR_gen5 = 1., deltaR_gen6 = 1., deltaR_gen7 = 1.;
	double ptratio = 1.;
	if(sampleType.find("vbfHInv")!=std::string::npos) {
	  deltaR_gen5 = TMath::Sqrt( (gen_eta->at(5)-genjet_eta->at(ijet))*(gen_eta->at(5)-genjet_eta->at(ijet)) + (gen_phi->at(5)-genjet_phi->at(ijet))*(gen_phi->at(5)-genjet_phi->at(ijet)) );
	  deltaR_gen6 = TMath::Sqrt( (gen_eta->at(6)-genjet_eta->at(ijet))*(gen_eta->at(6)-genjet_eta->at(ijet)) + (gen_phi->at(6)-genjet_phi->at(ijet))*(gen_phi->at(6)-genjet_phi->at(ijet)) );
	  deltaR_gen7 = TMath::Sqrt( (gen_eta->at(7)-genjet_eta->at(ijet))*(gen_eta->at(7)-genjet_eta->at(ijet)) + (gen_phi->at(7)-genjet_phi->at(ijet))*(gen_phi->at(7)-genjet_phi->at(ijet)) );
	  ptratio = gen_pt->at(refpart)/genjet_pt->at(ijet);
	  //jetcondn = (jetcondn and (deltaR_gen5<delRTh or deltaR_gen6<delRTh or deltaR_gen7<delRTh));
	  jetcondn = (jetcondn and (deltaR_gen5<delRTh or deltaR_gen6<delRTh or deltaR_gen7<delRTh) and (ptratio>0.8 and ptratio<1.2));
	}
	if(jetcondn){
	  //if(deltaR<delRTh and deltaR<minDeltaR and (deltaR_gen5<delRTh or deltaR_gen6<delRTh or deltaR_gen7<delRTh) ){
	  //if(deltaR<0.4 and deltaR<minDeltaR){
	  minDeltaR = deltaR;
	  minDeltaRindex = ipart;
	  minDR[0] = deltaR_gen5;
	  minDR[1] = deltaR_gen6;
	  minDR[2] = deltaR_gen7;
	}
      }
      //for(int ipart=0; ipart<partlist.size(); ipart++){
      for(int ipart=0; ipart<genlist.size(); ipart++){
	//JetPart p = partlist.at(ipart);
	JetPart p = genlist.at(ipart);
	if(ipart==minDeltaRindex){
	  JetPart jet;
	  jet.name = p.name;
	  jet.index = ijet ;
	  jetlist.push_back(jet);
	  
	  if(sampleType.find("vbfHInv")!=std::string::npos) {
	    int mindr = TMath::LocMin(3,minDR);
	    int minId = mindr + 5; //index of non-Higgs particles of the hardest subprocess in VBF sample
	    hGenvsJetPtDiff->Fill(gen_pt->at(minId)-genjet_pt->at(ijet));	  
	    pGenvsJetPtDiff->Fill( gen_pt->at(minId),  gen_pt->at(minId)/genjet_pt->at(ijet) );
	    pPtDiffvsDR->Fill( minDR[mindr],  gen_pt->at(minId)/genjet_pt->at(ijet) );
	    h2PtDiffvsDR->Fill( minDR[mindr],  gen_pt->at(minId)/genjet_pt->at(ijet) );
	    hGenvsJetDR->Fill( minDR[mindr]  );
	  }else{
	    hGenvsJetPtDiff->Fill(gen_pt->at(ipart)-genjet_pt->at(ijet));	  
	    pGenvsJetPtDiff->Fill( gen_pt->at(ipart),  gen_pt->at(ipart)/genjet_pt->at(ijet) );
	    pPtDiffvsDR->Fill( minDeltaR,  gen_pt->at(ipart)/genjet_pt->at(ijet) );
	    h2PtDiffvsDR->Fill( minDeltaR,  gen_pt->at(ipart)/genjet_pt->at(ijet) );
	    hGenvsJetDR->Fill( minDeltaR  );
	  }
	  if(doPrint){
	    if(sampleType.find("vbfHInv")!=std::string::npos)
	      std::cout << "ijet-ievent: " << std::setprecision(default_precision) << std::setw(4) <<ievent << ", ijet: " << std::setw(4) << ijet
			<<", Name: " << std::setw(10) << p.name
			<< ", deltaR : " << std::setprecision(3) << std::setw(5) << minDeltaR
			<< ", deltaR5 : " << std::setprecision(3) << std::setw(5) << minDR[0]
			<< ", deltaR6 : " << std::setprecision(3) << std::setw(5) << minDR[1]
			<< ", deltaR7 : " << std::setprecision(3) << std::setw(5) << minDR[2]
			<<" jet:(pt,eta,phi,energy) : (" << std::fixed << std::setprecision(2) << std::setw(8) << genjet_pt->at(ijet)
			<< ", " << std::setw(8) << genjet_eta->at(ijet) << ", " //<< std::setw(5) << genjet_exeta->at(ijet)
			<< ", " << std::setw(8) << (TMath::RadToDeg()*genjet_phi->at(ijet)) //<< ", " << std::setw(5) << (TMath::RadToDeg()*genjet_exphi->at(ijet))
			<< ", " << std::setw(8) << genjet_energy->at(ijet) << ") "
			<< std::defaultfloat
			<< std::endl;
	    else
	      std::cout << "ijet-ievent: " << std::setprecision(default_precision) << std::setw(4) <<ievent << ", ijet: " << std::setw(4) << ijet
			<<", Name: " << std::setw(10) << p.name
			<< ", deltaR : " << std::setprecision(3) << std::setw(5) << minDeltaR
			<<" jet:(pt,eta,phi,energy) : (" << std::fixed << std::setprecision(2) << std::setw(8) << genjet_pt->at(ijet)
			<< ", " << std::setw(8) << genjet_eta->at(ijet) << ", " //<< std::setw(5) << genjet_exeta->at(ijet)
			<< ", " << std::setw(8) << (TMath::RadToDeg()*genjet_phi->at(ijet)) //<< ", " << std::setw(5) << (TMath::RadToDeg()*genjet_exphi->at(ijet))
			<< ", " << std::setw(8) << genjet_energy->at(ijet) << ") "
			<< std::defaultfloat
			<< std::endl;

	  }//print condn
	}//minDelta condition
      }//part loop
    }//jet loop
    /////////////////////////////////////////////////////////////////////////////////////////
    hNJets->Fill(float(jetlist.size()));

    ////////////////////////Conditions for two jets ////////////////////////////////////
    bool isOneperSide = true ;
    if(jetlist.size()==2){
      isOneperSide = ((genjet_eta->at(jetlist.at(0).index) * genjet_eta->at(jetlist.at(1).index))<0.)?true:false;
      // double *allJets = new double[jetlist.size()];
      // int *sortedIdx = new int[jetlist.size()];
      // for(int ipjet=0; ipjet < jetlist.size() ; ipjet++ ){
      // 	int ijet = jetlist.at(ipjet).index ;
      // 	allJets[ipjet] = genjet_pt->at(ijet);
      // }
      // TMath::Sort(int(jetlist.size()), allJets, sortedIdx);
      // int HptIdx = jetlist.at(sortedIdx[0]).index ;
      // int MptIdx = jetlist.at(sortedIdx[1]).index ;
      // double eta12 = genjet_eta->at(HptIdx) * genjet_eta->at(MptIdx) ;
      // std::vector<JetPart> tJetlist;
      // for(int ipjet=0; ipjet < jetlist.size() ; ipjet++ ){
      // 	if(eta12<0.){
      // 	  if(ipjet==sortedIdx[0] or ipjet==sortedIdx[1]){
      // 	    JetPart jet;
      // 	    jet.name = jetlist.at(ipjet).name ;
      // 	    jet.index = jetlist.at(ipjet).index ;
      // 	    tJetlist.push_back(jet);
      // 	  }
      // 	}else{
      // 	  if(ipjet==sortedIdx[0]){
      // 	    JetPart jet;
      // 	    jet.name = jetlist.at(ipjet).name ;
      // 	    jet.index = jetlist.at(ipjet).index ;
      // 	    tJetlist.push_back(jet);
      // 	  }	  
      // 	}
      // }
      // jetlist.clear();
      // for(int ipjet=0; ipjet < tJetlist.size() ; ipjet++ ){
      // 	int ijet = tJetlist.at(ipjet).index ;
      // 	JetPart jet;
      // 	jet.name = tJetlist.at(ipjet).name ;
      // 	jet.index = tJetlist.at(ipjet).index ;
      // 	jetlist.push_back(jet);
      // 	if(doPrint)
      // 	  std::cout << "filtered-ijet-ievent: " << std::setprecision(default_precision) << std::setw(4) <<ievent << ", ijet: " << std::setw(4) << ijet
      // 		    <<", Name: " << std::setw(10) << jet.name
      // 		    <<" jet:(pt,eta,phi,energy) : (" << std::fixed << std::setprecision(2) << std::setw(8) << genjet_pt->at(ijet)
      // 		    << ", " << std::setw(8) << genjet_eta->at(ijet) << ", " //<< std::setw(5) << genjet_exeta->at(ijet)
      // 		    << ", " << std::setw(8) << (TMath::RadToDeg()*genjet_phi->at(ijet)) //<< ", " << std::setw(5) << (TMath::RadToDeg()*genjet_exphi->at(ijet))
      // 		    << ", " << std::setw(8) << genjet_energy->at(ijet) << ") "
      // 		    << std::defaultfloat
      // 		    << std::endl;
      // }      
      // tJetlist.clear();      
      // delete allJets;
      // delete sortedIdx;
    }
    ////////////////////////////////////////////////////////////////////////////////////////
    
    // ///////////////////////// Conditions for three jets ////////////////////////////////////
    // if(jetlist.size()>2){
    //   double *allJets = new double[jetlist.size()];
    //   int *sortedIdx = new int[jetlist.size()];
    //   for(int ipjet=0; ipjet < jetlist.size() ; ipjet++ ){
    // 	int ijet = jetlist.at(ipjet).index ;
    // 	allJets[ipjet] = genjet_pt->at(ijet);
    //   }
    //   TMath::Sort(int(jetlist.size()), allJets, sortedIdx);
    //   int HptIdx = jetlist.at(sortedIdx[0]).index ;
    //   int MptIdx = jetlist.at(sortedIdx[1]).index ;
    //   int SptIdx = jetlist.at(sortedIdx[2]).index ;
    //   double eta12 = genjet_eta->at(HptIdx) * genjet_eta->at(MptIdx) ;
    //   double eta13 = genjet_eta->at(HptIdx) * genjet_eta->at(SptIdx) ;
    //   std::vector<JetPart> tJetlist;
    //   for(int ipjet=0; ipjet < jetlist.size() ; ipjet++ ){
    // 	if(eta12<0.){
    // 	  if(ipjet==sortedIdx[0] or ipjet==sortedIdx[1]){
    // 	    JetPart jet;
    // 	    jet.name = jetlist.at(ipjet).name ;
    // 	    jet.index = jetlist.at(ipjet).index ;
    // 	    tJetlist.push_back(jet);
    // 	  }
    // 	}else if(eta13<0.){
    // 	  if(ipjet==sortedIdx[0] or ipjet==sortedIdx[2]){
    // 	    JetPart jet;
    // 	    jet.name = jetlist.at(ipjet).name ;
    // 	    jet.index = jetlist.at(ipjet).index ;
    // 	    tJetlist.push_back(jet);
    // 	  }
    // 	}else{
    // 	  if(ipjet==sortedIdx[0]){
    // 	    JetPart jet;
    // 	    jet.name = jetlist.at(ipjet).name ;
    // 	    jet.index = jetlist.at(ipjet).index ;
    // 	    tJetlist.push_back(jet);
    // 	  }	  
    // 	}
    //   }
    //   jetlist.clear();
    //   for(int ipjet=0; ipjet < tJetlist.size() ; ipjet++ ){
    // 	int ijet = tJetlist.at(ipjet).index ;
    // 	JetPart jet;
    // 	jet.name = tJetlist.at(ipjet).name ;
    // 	jet.index = tJetlist.at(ipjet).index ;
    // 	jetlist.push_back(jet);
    // 	if(doPrint)
    // 	  std::cout << "filtered-ijet-ievent: " << std::setprecision(default_precision) << std::setw(4) <<ievent << ", ijet: " << std::setw(4) << ijet
    // 		    <<", Name: " << std::setw(10) << jet.name
    // 		    <<" jet:(pt,eta,phi,energy) : (" << std::fixed << std::setprecision(2) << std::setw(8) << genjet_pt->at(ijet)
    // 		    << ", " << std::setw(8) << genjet_eta->at(ijet) << ", " //<< std::setw(5) << genjet_exeta->at(ijet)
    // 		    << ", " << std::setw(8) << (TMath::RadToDeg()*genjet_phi->at(ijet)) //<< ", " << std::setw(5) << (TMath::RadToDeg()*genjet_exphi->at(ijet))
    // 		    << ", " << std::setw(8) << genjet_energy->at(ijet) << ") "
    // 		    << std::defaultfloat
    // 		    << std::endl;
    //   }      
    //   tJetlist.clear();      
    //   delete allJets;
    //   delete sortedIdx;
    // }
    // ////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////
    double xoz_posEta = -10., yoz_posEta = -10.;
    double xoz_negEta = -10., yoz_negEta = -10.;
    double genjetpt_posEta =  -1., genjetpt_negEta = -1.;
    int njets = 0, njetsp = 0, njetsm = 0;
    for(int ipjet=0; ipjet < jetlist.size() ; ipjet++ ){
      int ijet = jetlist.at(ipjet).index ;
      if(fabs(genjet_eta->at(ijet))<=1.321 or fabs(genjet_eta->at(ijet))>=3.152) continue;
      double roz = 1/sinh(genjet_eta->at(ijet)) ; //tan(2*atan(exp(-1.0*genjet_eta->at(ijet))));
      double xoz = roz*cos(genjet_phi->at(ijet));
      double yoz = roz*sin(genjet_phi->at(ijet));
      if(genjet_eta->at(ijet)>0. and nofHGCalJetsPos==1 ){ //and nofJetsLEPhiPos==0
	xoz_posEta = xoz; yoz_posEta = yoz;
	genjetpt_posEta = genjet_pt->at(ijet);
	njets++;
	njetsp++;
      }else if(genjet_eta->at(ijet)<0. and nofHGCalJetsNeg==1 ){ //and nofJetsLEPhiNeg==0
	xoz_negEta = xoz; yoz_negEta = yoz;
	genjetpt_negEta = genjet_pt->at(ijet);
	njets++;
	njetsm++;
      }
      if(doPrint)
	std::cout << "ijet1-ievent: " << std::setprecision(default_precision) << std::setw(4) <<ievent << ", ijet: " << std::setw(4) << ijet
		  <<", Name: " << std::setw(10) << jetlist.at(ipjet).name
		  <<" jet:(pt,eta,phi,energy) : (" << std::fixed << std::setprecision(2) << std::setw(8) << genjet_pt->at(ijet)
		  << ", " << std::setw(8) << genjet_eta->at(ijet) << ", " //<< std::setw(5) << genjet_exeta->at(ijet)
		  << ", " << std::setw(8) << (TMath::RadToDeg()*genjet_phi->at(ijet)) //<< ", " << std::setw(5) << (TMath::RadToDeg()*genjet_exphi->at(ijet))
		  << ", " << std::setw(8) << genjet_energy->at(ijet) << ") "
		  <<" jet:(roz,xoz,yoz) : (" << std::fixed << std::setprecision(2) << std::setw(5) << roz
		  <<", " << std::fixed << std::setprecision(2) << std::setw(5) << xoz
		  <<", " << std::fixed << std::setprecision(2) << std::setw(5) << yoz
		  << ") "
		  << std::defaultfloat
		  << std::endl;
      
    }
    ////////////////////////////////////////////////////////////////////////////////////////
    
    if(njets>2){
      std::cout << "Issue: Event: " << ievent << ", njetsp: " << njetsp << ", njetsm: " << njetsm << std::endl;
    }
    
    /////////////////////////// Conditions for 0 jets and others ////////////////////////////////////

    if(njets==0 or (jetlist.size()==2 and !isOneperSide) or njets>2){
      //if(jetlist.size()==0 or !isOneperSide or nofHGCalJetsPos>1 or nofHGCalJetsNeg>1){
      //if(jetlist.size()==0 or jetlist.size()>2 or !isOneperSide){
      //if(jetlist.size()==0){
      taudlist.clear();
      taugdlist.clear();
      genlist.clear();
      partlist.clear();
      //jetlist.clear();
      
      gen_pt->clear();
      gen_eta->clear();
      gen_phi->clear();
      gen_pdgid->clear();
      gen_energy->clear();
      gen_charge->clear();
      gen_status->clear();
      for(uint32_t i=0;i<gen_daughters->size();i++) (gen_daughters->at(i)).clear();
      gen_daughters->clear();
      genpart_pt->clear();
      genpart_eta->clear();
      genpart_phi->clear();
      genpart_pid->clear();
      genpart_mother->clear();
      genpart_energy->clear();
      genpart_gen->clear();
      genpart_fromBeamPipe->clear();
      genpart_exeta->clear();
      genpart_exphi->clear();
      genjet_pt->clear();
      genjet_eta->clear();
      genjet_phi->clear();
      genjet_energy->clear();
      tc_pt->clear();
      tc_mipPt->clear();
      tc_eta->clear();
      tc_phi->clear();
      tc_layer->clear();
      tc_x->clear();
      tc_y->clear();
      tc_z->clear();
      tc_energy->clear();
      tc_data->clear();
      tc_subdet->clear();
      tc_wafertype->clear();
      tc_uncompressedCharge->clear();
      cl3d_pt->clear();
      cl3d_phi->clear();
      cl3d_eta->clear();
      continue;
    }
    ////////////////////////////////////////////////////////////////////////////////////////
    
    hNJetsF1->Fill(float(njets));
    
    
    // std::vector<TCOrder> rearrgdtcs;
    // for(unsigned itc=0;itc<tc_pt->size();itc++){
    //   TCOrder tco;
    //   tco.tcid = itc;
    //   tco.layer = uint32_t(tc_layer->at(itc));
    //   tco.z = tc_z->at(itc);
    //   tco.eta = tc_eta->at(itc);
    //   rearrgdtcs.push_back(tco);
    // }
    // std::vector<TCOrder>::iterator itstart = rearrgdtcs.begin();
    // std::vector<TCOrder>::iterator itstop = rearrgdtcs.end();
    // std::sort(itstart,itstop,comp);
    
    /////////////////////////////////////////////////////////////////////////////////////////
    double totTCpt_posEta_iter1 = 0.,totTCpt_negEta_iter1 = 0.;
    double totTCpt_posEta_CEE_iter1 = 0.,totTCpt_negEta_CEE_iter1 = 0.;
    double totTCpt_posEta_core_iter1 = 0.,totTCpt_negEta_core_iter1 = 0.;
    double totTCpt_posEta_core_rozcut0p15_iter1 = 0.,totTCpt_negEta_core_rozcut0p15_iter1 = 0.;
    double totTCpt_posEta_core_mpptcut2_iter1 = 0.,totTCpt_negEta_core_mpptcut2_iter1 = 0.;
    double totTCpt_posEta_CEE_core_iter1 = 0.,totTCpt_negEta_CEE_core_iter1 = 0.;
    for(unsigned itc=0;itc<tc_pt->size();itc++){
      
      double minRoz = 1., clstPt = 0, clstEta = 0;
      for(int ipjet=0; ipjet < jetlist.size() ; ipjet++ ){
	int ijet = jetlist.at(ipjet).index ;
	double gjroz = 1/sinh(genjet_eta->at(ijet)) ; 
	double gjxoz = gjroz*cos(genjet_phi->at(ijet));
	double gjyoz = gjroz*sin(genjet_phi->at(ijet));
	double dClusGenRoz = sqrt( (tc_x->at(itc)/tc_z->at(itc) - gjxoz)*(tc_x->at(itc)/tc_z->at(itc) - gjxoz) + (tc_y->at(itc)/tc_z->at(itc) - gjyoz)*(tc_y->at(itc)/tc_z->at(itc) - gjyoz) );
	if(dClusGenRoz<minRoz and tc_eta->at(itc)*genjet_eta->at(ijet)>0.){
	  minRoz = dClusGenRoz;
	  clstPt = genjet_pt->at(ijet);
	  clstEta = genjet_eta->at(ijet);
	}
      }
      double tcptMeV = tc_pt->at(itc) * 1.e3;
      double tcpt = (tcptMeV<3.)?3:tcptMeV; //VBF PU0      
      double ptcorr = (tc_layer->at(itc)<=26)?1.07:1.14;//1.13 for Hadronic part //f2->Eval(tcpt);
      //double ptcorr = f2->Eval(tcpt);
      
      // double tcpt = (tcptMeV<0.1)?0.1:tcptMeV; //VBF PU200
      // double ptcorr_high = fhigh->Eval(tcpt);
      // // double ptcorr_high1 = fhigh1->Eval(tcpt);
      // if(tcpt>400.) ptcorr = ptcorr* ptcorr_high;
      // // if(tcpt>200.) ptcorr = ptcorr* ptcorr_high1;
      
      tcRoz->Fill(minRoz);
      htcMipPt->Fill(tc_mipPt->at(itc));
      tcRozPtwt->Fill(minRoz, tc_pt->at(itc));

      //if(fabs(clstEta)>etaMinCore and fabs(clstEta)<etaMaxCore){ 
	if(tc_eta->at(itc)<0 and genjetpt_negEta>-1.){
	  totTCpt_negEta_iter1 += ptcorr*tc_pt->at(itc);
	  if(tc_layer->at(itc)<=26) totTCpt_negEta_CEE_iter1 += ptcorr*tc_pt->at(itc);    
	}else if(tc_eta->at(itc)>0. and genjetpt_posEta>-1.){
	  totTCpt_posEta_iter1 += ptcorr*tc_pt->at(itc);
	  if(tc_layer->at(itc)<=26) totTCpt_posEta_CEE_iter1 += ptcorr*tc_pt->at(itc);
	}
	//}
      
      if(fabs(clstEta)>etaMinCore and fabs(clstEta)<etaMaxCore){ //only add for the plateau region	
	if(tc_eta->at(itc)<0 and genjetpt_negEta>-1.){
	  totTCpt_negEta_core_iter1 += ptcorr*tc_pt->at(itc);
	  if(tc_layer->at(itc)<=26) totTCpt_negEta_CEE_core_iter1 += ptcorr*tc_pt->at(itc);
	}else if(tc_eta->at(itc)>0. and genjetpt_posEta>-1.){
	  totTCpt_posEta_core_iter1 += ptcorr*tc_pt->at(itc);
	  if(tc_layer->at(itc)<=26) totTCpt_posEta_CEE_core_iter1 += ptcorr*tc_pt->at(itc);
	}
      }
      
      if(fabs(clstEta)>etaMinCore and fabs(clstEta)<etaMaxCore and minRoz<tcRoZCut){ //only add for the plateau region
	tcRozCore->Fill(minRoz);
	if(tc_eta->at(itc)<0 and genjetpt_negEta>-1.){
	  totTCpt_negEta_core_rozcut0p15_iter1 += ptcorr*tc_pt->at(itc);
	}else if(tc_eta->at(itc)>0. and genjetpt_posEta>-1.){
	  totTCpt_posEta_core_rozcut0p15_iter1 += ptcorr*tc_pt->at(itc);
	}
      }
      
      if(fabs(clstEta)>etaMinCore and fabs(clstEta)<etaMaxCore and tc_mipPt->at(itc)>2.){ //only add for the plateau region	
	if(tc_eta->at(itc)<0 and genjetpt_negEta>-1.){
	  totTCpt_negEta_core_mpptcut2_iter1 += ptcorr*tc_pt->at(itc);	  
	}else if(tc_eta->at(itc)>0. and genjetpt_posEta>-1.){
	  totTCpt_posEta_core_mpptcut2_iter1 += ptcorr*tc_pt->at(itc);
	}
      }
      
      
    }    
    /////////////////////////////////////////////////////////////////////////////////////////
    double totTCpt_posEta = 0.,totTCpt_negEta = 0.;
    double totTCpt_posEta_CEE = 0.,totTCpt_negEta_CEE = 0.;
    double totTCpt_posEta_layer1 = 0.,totTCpt_negEta_layer1 = 0.;
    uint32_t nofPosEtaTCs = 0, nofNegEtaTCs = 0;
    double tcXoZposEta = 0.,tcYoZposEta = 0., tcXoZnegEta = 0.,tcYoZnegEta = 0.;
    double maxroz_posEta = -1.0, maxroz_negEta = -1.0;
    double tcPhi_posEta = 0.,tcPhi_negEta = 0.;
    double tcEta_posEta = 0.,tcEta_negEta = 0.;
    double maxTCPt_posEta = -1.0, maxTCPt_negEta = -1.0;
    int maxTCPt_posEta_layer = -1, maxTCPt_negEta_layer = -1;
    int maxTCPt_posEta_detType = -1, maxTCPt_negEta_detType = -1;
    double maxTCPt_posEta_genPt = -1, maxTCPt_negEta_genPt = -1;
    double maxTCPt_posEta_genEta = -1, maxTCPt_negEta_genEta = -1;
    int detType = -1;
    int subdet = -1, wafertype = -1;
    //std::vector<TPGTriggerCellWord> vTcw[6];
    std::vector<TPGTCBits> vTcw[6];
    // for(unsigned itco=0;itco<tc_pt->size();itco++){
    //   unsigned itc = rearrgdtcs.at(itco).tcid ; 
    for(unsigned itc=0;itc<tc_pt->size();itc++){
      double z(fabs(tc_z->at(itc)));
      float phi_deg = TMath::RadToDeg() * tc_phi->at(itc) ;
      uint16_t sec0(7),sec1(7);
      
      //The following segmentation is applied for trigger cells
      if(phi_deg>=0. and phi_deg<=60.) {sec0 = 2; sec1 = 0;}
      else if(phi_deg>60. and phi_deg<=120.) {sec0 = sec1 = 0;}
      else if(phi_deg>120. and phi_deg<=180.) {sec0 = 0; sec1 = 1;}
      else if(phi_deg>=-180. and phi_deg<=-120.) {sec0 = sec1 = 1;}
      else if(phi_deg>-120. and phi_deg<=-60.) {sec0 = 1; sec1 = 2;}
      else {sec0 = sec1 = 2;}
      
      if(tc_z->at(itc)<0.){
	switch(sec0){
	case 0:
	  sec0 += 3;
	  break;
	case 1:
	  sec0 += 4;
	  break;
	case 2:
	  sec0 += 2;
	  break;
	default:
	  ;
	}
	switch(sec1){
	case 0:
	  sec1 += 3;
	  break;
	case 1:
	  sec1 += 4;
	  break;
	case 2:
	  sec1 += 2;
	  break;
	default:
	  ;
	}	
      }
      detType = -1;
      subdet = tc_subdet->at(itc);
      wafertype = tc_wafertype->at(itc);
      if(subdet==1 and wafertype==0) detType = 1;
      else if(subdet==1 and wafertype==1) detType = 2;
      else if(subdet==1 and wafertype==2) detType = 3;
      else if(subdet==2 and wafertype==0) detType = 4;
      else if(subdet==2 and wafertype==1) detType = 5;
      else if(subdet==2 and wafertype==2) detType = 6;
      else if(subdet==10) detType = 7;
      
      double minRoz = 1., clstPt = 0, clstEta = 0;
      for(int ipjet=0; ipjet < jetlist.size() ; ipjet++ ){
	int ijet = jetlist.at(ipjet).index ;
	double gjroz = 1/sinh(genjet_eta->at(ijet)) ; 
	double gjxoz = gjroz*cos(genjet_phi->at(ijet));
	double gjyoz = gjroz*sin(genjet_phi->at(ijet));
	double dClusGenRoz = sqrt( (tc_x->at(itc)/tc_z->at(itc) - gjxoz)*(tc_x->at(itc)/tc_z->at(itc) - gjxoz) + (tc_y->at(itc)/tc_z->at(itc) - gjyoz)*(tc_y->at(itc)/tc_z->at(itc) - gjyoz) );
	if(dClusGenRoz<minRoz and tc_eta->at(itc)*genjet_eta->at(ijet)>0.){
	  minRoz = dClusGenRoz;
	  clstPt = genjet_pt->at(ijet);
	  clstEta = genjet_eta->at(ijet);
	}
      }
      tcRozDet[detType-1]->Fill(minRoz);
      
      double tcptMeV = tc_pt->at(itc) * 1.e3;
      double tcpt = (tcptMeV<3.)?3:tcptMeV;
      double ptcorr = (tc_layer->at(itc)<=26)?1.07:1.14;//1.07:1.13f2->Eval(tcpt);
      //double ptcorr = f2->Eval(tcpt);
      tcpt = ptcorr*tc_pt->at(itc);
      tcptMeV = tcpt*1.e3;
      
      // if(doPrint){
      // 	std::cout << "itc-ievent: " << std::fixed << std::setprecision(2) << std::setw(4) << ievent
      // 		  << ", itc: " << std::setw(5) << itc
      // 		  << ", layer: " << std::setw(2) << int(tc_layer->at(itc))
      // 		  <<", (x,y,z): (" << std::setw(7) << tc_x->at(itc) << ", " << std::setw(7) << tc_y->at(itc) << ", " << std::setw(7) << tc_z->at(itc) << ")"
      // 		  <<", (pt,eta,phi,energy): (" << std::fixed << std::setprecision(5) << std::setw(8) << tc_pt->at(itc)
      // 		  << ", " << std::setw(8) << tc_eta->at(itc) 
      // 		  << ", " << std::setw(8) << (TMath::RadToDeg()*tc_phi->at(itc)) 
      // 		  << ", " << std::setw(8) << tc_energy->at(itc)*1.e3 << ") "
      // 		  <<", (reso, ptcorr,new_pt) : ( " << std::fixed << std::setprecision(5) << std::setw(8) << reso
      // 		  << ", " << std::setw(8) << ptcorr
      // 		  << ", " << std::setw(8) << tc_pt->at(itc)*ptcorr << ") "
      // 		  << std::defaultfloat
      // 		  << std::endl;
      // }
      
      for (uint32_t isect = 0 ; isect < 3 ; isect++ ){
	uint32_t addisect = 0;
	if(tc_z->at(itc)>0.0) addisect = 3;
	tcf0.setZero();
	tcf0.setROverZPhiF(tc_x->at(itc)/z,tc_y->at(itc)/z,isect+addisect);
	if(tcf0.getXOverZF()>=0.0){
	  float scale = ptcorr;
	  //if(tc_layer->at(itc)==27 or tc_layer->at(itc)==29) scale = 1.5;
	  tcf0.setEnergyGeV(scale * tc_pt->at(itc));
	  tcf0.setLayer(tc_layer->at(itc));
	  //if(doPrint) tcf0.print();
	  //vTcw[isect+addisect].push_back(tcf0);
	  if((tc_eta->at(itc)<0. and genjetpt_negEta>-1.) or (tc_eta->at(itc)>0. and genjetpt_posEta>-1.))
	    vTcw[isect+addisect].push_back(tcf0);
	}
      }
      
      // bool hgcalEta = (fabs(tc_eta->at(itc))>1.321 and fabs(tc_eta->at(itc))<3.152) ? true : false;
      //bool hgcalInnerEta = (fabs(tc_eta->at(itc))>etaMinCore and fabs(tc_eta->at(itc))<etaMaxCore) ? true : false;
      bool hgcalInnerEta = (fabs(clstEta)>etaMinCore and fabs(clstEta)<etaMaxCore) ? true : false;
      int binEta = 1, binPt = 1, binPtLog = 1;
      GetHistoBin(hJetEtaTCPtBin, tc_eta->at(itc), tcpt, binEta, binPt);
      GetHistoBin(hJetEtaTCPtBinLog, tc_eta->at(itc), tcptMeV, binEta, binPtLog);
      
      hTCPt->Fill(tc_pt->at(itc));
      hTCEta->Fill(tc_eta->at(itc));
      h2TCEtaPhi->Fill(TMath::RadToDeg()*tc_phi->at(itc), tc_eta->at(itc));
      if(tc_eta->at(itc)<0. and genjetpt_negEta>-1.){
	hTCPhiNeg->Fill(TMath::RadToDeg()*tc_phi->at(itc));
	//if(hgcalInnerEta){
	// totTCpt_negEta += tc_pt->at(itc);
	// if(tc_layer->at(itc)<=26) totTCpt_negEta_CEE += tc_pt->at(itc);
	//}
	//if(minRoz<tcRoZCut){
	totTCpt_negEta += tcpt;
	if(tc_layer->at(itc)<=26) totTCpt_negEta_CEE += tcpt;
	//}
	tcXoZnegEta += tcpt*tc_x->at(itc)/tc_z->at(itc) ;
	tcYoZnegEta += tcpt*tc_y->at(itc)/tc_z->at(itc) ;
	tcPhi_negEta += tcpt*tc_phi->at(itc) ;
	tcEta_negEta += tcpt*tc_eta->at(itc) ;
	nofNegEtaTCs++;
	double droz = sqrt( (tc_x->at(itc)/tc_z->at(itc)-xoz_negEta)*(tc_x->at(itc)/tc_z->at(itc)-xoz_negEta) + (tc_y->at(itc)/tc_z->at(itc)-yoz_negEta)*(tc_y->at(itc)/tc_z->at(itc)-yoz_negEta) );
	if(droz>maxroz_negEta) maxroz_negEta = droz;
	if(tcpt>maxTCPt_negEta){
	  //if(tc_uncompressedCharge->at(itc)>maxTCPt_negEta){
	  //maxTCPt_negEta = tc_uncompressedCharge->at(itc);//tcpt;
	  maxTCPt_negEta = tcpt;
	  maxTCPt_negEta_layer = tc_layer->at(itc);
	  maxTCPt_negEta_detType = detType;
	  maxTCPt_negEta_genPt = clstPt;
	  maxTCPt_negEta_genEta = clstEta;
	}
	//if(hgcalInnerEta and minRoz<tcRoZCut and clstPt>(pt_clusThresh*0.5)){
	if(hgcalInnerEta and minRoz<tcRoZCut){
	  h2TCPtvsGenPtByTCPtsum->Fill(tcptMeV, (clstPt/totTCpt_negEta_iter1));
	  pTCPtvsGenPtByTCPtsum->Fill(tcptMeV, (clstPt/totTCpt_negEta_iter1));
	  pTCPtvsGenPtByTCPtsumLogx->Fill(tcptMeV, (clstPt/totTCpt_negEta_iter1));
	  pTCPtvsGenPtByTCPtsumDet[detType-1]->Fill(tcptMeV, (clstPt/totTCpt_negEta_iter1));
	  hGenjetByTCPtSumVsTCPt1DLog[binPtLog-1]->Fill( clstPt/totTCpt_negEta_iter1 );
     	  hGenjetByTCPtSumVsTCPt2DDetLog[detType-1][binPtLog-1]->Fill( clstPt/totTCpt_negEta_iter1 );
	  hGenjetByTCPtSumVsTCPt2DLog[binEta-1][binPtLog-1]->Fill( clstPt/totTCpt_negEta_iter1 );
	  //hGenjetByTCPtSumVsTCPt3DLog[detType-1][binEta-1][binPtLog-1]->Fill( clstPt/totTCpt_negEta_iter1 );
	}
	h2GenjetTCPtsumvsGenPt->Fill(clstPt, clstPt/totTCpt_negEta_iter1); 
	h2GenjetTCPtsumvsGenEta->Fill(fabs(clstEta), clstPt/totTCpt_negEta_iter1);
	pGenjetTCPtsumvsGenPt->Fill(clstPt, clstPt/totTCpt_negEta_iter1); 
	pGenjetTCPtsumvsGenEta->Fill(fabs(clstEta), clstPt/totTCpt_negEta_iter1);
	if(hgcalInnerEta){
	  h2GenjetTCPtsumvsGenPt_CoreEta->Fill(clstPt, clstPt/totTCpt_negEta_CEE_core_iter1);
	  pGenjetTCPtsumvsGenPt_CoreEta->Fill(clstPt, clstPt/totTCpt_negEta_CEE_core_iter1);
	}
	if(tc_layer->at(itc)<=26){
	  h2GenjetTCPtsumvsGenPt_CEE->Fill(clstPt, clstPt/totTCpt_negEta_CEE_iter1);
	  h2GenjetTCPtsumvsGenEta_CEE->Fill(fabs(clstEta), clstPt/totTCpt_negEta_CEE_iter1);
	  pGenjetTCPtsumvsGenPt_CEE->Fill(clstPt, clstPt/totTCpt_negEta_CEE_iter1);
	  pGenjetTCPtsumvsGenEta_CEE->Fill(fabs(clstEta), clstPt/totTCpt_negEta_CEE_iter1);
	  if(hgcalInnerEta){
	    h2GenjetTCPtsumvsGenPt_CEE_CoreEta->Fill(clstPt, clstPt/totTCpt_negEta_CEE_core_iter1);
	    pGenjetTCPtsumvsGenPt_CEE_CoreEta->Fill(clstPt, clstPt/totTCpt_negEta_CEE_core_iter1);
	    h2GenjetTCPtsumvsTCPtsum_CEE_CoreEta->Fill(totTCpt_negEta_CEE_core_iter1, clstPt/totTCpt_negEta_CEE_core_iter1);
	    pGenjetTCPtsumvsTCPtsum_CEE_CoreEta->Fill(totTCpt_negEta_CEE_core_iter1, clstPt/totTCpt_negEta_CEE_core_iter1);	    
	  }
	}
      }else if(tc_eta->at(itc)>0. and genjetpt_posEta>-1.){
	//if(hgcalInnerEta){
	// totTCpt_posEta += tc_pt->at(itc);
	// if(tc_layer->at(itc)<=26) totTCpt_posEta_CEE += tc_pt->at(itc);
	//}
	//if(minRoz<tcRoZCut){
	hTCPhiPos->Fill(TMath::RadToDeg()*tc_phi->at(itc));
	totTCpt_posEta += tcpt;
	if(tc_layer->at(itc)<=26) totTCpt_posEta_CEE += tcpt;
	//}
	tcXoZposEta += tcpt*tc_x->at(itc)/tc_z->at(itc) ;
	tcYoZposEta += tcpt*tc_y->at(itc)/tc_z->at(itc) ;
	tcPhi_posEta += tcpt*tc_phi->at(itc) ;
	tcEta_posEta += tcpt*tc_eta->at(itc) ;
	nofPosEtaTCs++;
	double droz = sqrt( (tc_x->at(itc)/tc_z->at(itc)-xoz_posEta)*(tc_x->at(itc)/tc_z->at(itc)-xoz_posEta) + (tc_y->at(itc)/tc_z->at(itc)-yoz_posEta)*(tc_y->at(itc)/tc_z->at(itc)-yoz_posEta) );
	if(droz>maxroz_posEta) maxroz_posEta = droz;
	if(tcpt>maxTCPt_posEta){
	  //if(tc_uncompressedCharge->at(itc)>maxTCPt_posEta){
	  //maxTCPt_posEta = tc_uncompressedCharge->at(itc);//tcpt;
	  maxTCPt_posEta = tcpt;
	  maxTCPt_posEta_layer = tc_layer->at(itc);
	  maxTCPt_posEta_detType = detType;
	  maxTCPt_posEta_genPt = clstPt;
	  maxTCPt_posEta_genEta = clstEta;
	}
	//if(hgcalInnerEta and minRoz<tcRoZCut and clstPt>(pt_clusThresh*0.5)){
	if(hgcalInnerEta and minRoz<tcRoZCut){
	  h2TCPtvsGenPtByTCPtsum->Fill(tcptMeV, (clstPt/totTCpt_posEta_iter1));
	  pTCPtvsGenPtByTCPtsum->Fill(tcptMeV, (clstPt/totTCpt_posEta_iter1));
	  pTCPtvsGenPtByTCPtsumLogx->Fill(tcptMeV, (clstPt/totTCpt_posEta_iter1));
	  pTCPtvsGenPtByTCPtsumDet[detType-1]->Fill(tcptMeV, (clstPt/totTCpt_posEta_iter1));
	  hGenjetByTCPtSumVsTCPt1DLog[binPtLog-1]->Fill( clstPt/totTCpt_posEta_iter1 );
     	  hGenjetByTCPtSumVsTCPt2DDetLog[detType-1][binPtLog-1]->Fill( clstPt/totTCpt_posEta_iter1 );
	  hGenjetByTCPtSumVsTCPt2DLog[binEta-1][binPtLog-1]->Fill( clstPt/totTCpt_posEta_iter1 );
	  //hGenjetByTCPtSumVsTCPt3DLog[detType-1][binEta-1][binPtLog-1]->Fill( clstPt/totTCpt_posEta_iter1 );
	}
	h2GenjetTCPtsumvsGenPt->Fill(clstPt, clstPt/totTCpt_posEta_iter1); 
	h2GenjetTCPtsumvsGenEta->Fill(fabs(clstEta), clstPt/totTCpt_posEta_iter1);
	pGenjetTCPtsumvsGenPt->Fill(clstPt, clstPt/totTCpt_posEta_iter1); 
	pGenjetTCPtsumvsGenEta->Fill(fabs(clstEta), clstPt/totTCpt_posEta_iter1);
	if(hgcalInnerEta){
	  h2GenjetTCPtsumvsGenPt_CoreEta->Fill(clstPt, clstPt/totTCpt_posEta_CEE_core_iter1);
	  pGenjetTCPtsumvsGenPt_CoreEta->Fill(clstPt, clstPt/totTCpt_posEta_CEE_core_iter1);
	}
	if(tc_layer->at(itc)<=26){
	  h2GenjetTCPtsumvsGenPt_CEE->Fill(clstPt, clstPt/totTCpt_posEta_CEE_iter1);
	  h2GenjetTCPtsumvsGenEta_CEE->Fill(fabs(clstEta), clstPt/totTCpt_posEta_CEE_iter1);
	  pGenjetTCPtsumvsGenPt_CEE->Fill(clstPt, clstPt/totTCpt_posEta_CEE_iter1);
	  pGenjetTCPtsumvsGenEta_CEE->Fill(fabs(clstEta), clstPt/totTCpt_posEta_CEE_iter1);
	  if(hgcalInnerEta){
	    h2GenjetTCPtsumvsGenPt_CEE_CoreEta->Fill(clstPt, clstPt/totTCpt_posEta_CEE_core_iter1);
	    pGenjetTCPtsumvsGenPt_CEE_CoreEta->Fill(clstPt, clstPt/totTCpt_posEta_CEE_core_iter1);
	    h2GenjetTCPtsumvsTCPtsum_CEE_CoreEta->Fill(totTCpt_posEta_CEE_core_iter1, clstPt/totTCpt_posEta_CEE_core_iter1);
	    pGenjetTCPtsumvsTCPtsum_CEE_CoreEta->Fill(totTCpt_posEta_CEE_core_iter1, clstPt/totTCpt_posEta_CEE_core_iter1);
	  }
	}
      }
      
    }//end of TC loop
    
    //rearrgdtcs.clear();
    if(doPrint){
      int nofTCTots = 0;
      if(totTCpt_negEta_iter1>0.) nofTCTots++;
      if(totTCpt_posEta_iter1>0.) nofTCTots++;
      std::cout << "itc-ievent: " << std::fixed << std::setprecision(2) << std::setw(4) << ievent
		<< ", nofTCTots: " << std::setw(5) << nofTCTots
		<< ", totTCpt_negEta_iter1: " << std::setw(5) << totTCpt_negEta_iter1
		<< ", totTCpt_posEta_iter1: " << std::setw(5) << totTCpt_posEta_iter1
		<< std::defaultfloat
		<< std::endl;
      
    }
    
    for(int ipjet=0; ipjet < jetlist.size() ; ipjet++ ){
      
      int ijet = jetlist.at(ipjet).index ;

      //if(fabs(genjet_eta->at(ijet))>etaMinCore and fabs(genjet_eta->at(ijet))<etaMaxCore){
      //if(fabs(genjet_eta->at(ijet))>1.6 and fabs(genjet_eta->at(ijet))<3.0){
	if(genjet_eta->at(ijet)<0 and genjetpt_negEta>-1. and totTCpt_negEta_iter1>0.){	
	  h2TCvsGenPt->Fill(genjetpt_negEta, totTCpt_negEta_iter1);
	  pTCvsGenPt->Fill(genjetpt_negEta, totTCpt_negEta_iter1);
	  h2TCvsGenPt_CEE->Fill(genjetpt_negEta, totTCpt_negEta_CEE_iter1);
	  pTCvsGenPt_CEE->Fill(genjetpt_negEta, totTCpt_negEta_CEE_iter1);
	}else if(genjet_eta->at(ijet)>0 and genjetpt_posEta>-1. and totTCpt_posEta_iter1>0.){
	  h2TCvsGenPt->Fill(genjetpt_posEta, totTCpt_posEta_iter1);
	  pTCvsGenPt->Fill(genjetpt_posEta, totTCpt_posEta_iter1);
	  h2TCvsGenPt_CEE->Fill(genjetpt_posEta, totTCpt_posEta_CEE_iter1);
	  pTCvsGenPt_CEE->Fill(genjetpt_posEta, totTCpt_posEta_CEE_iter1);
	}
	//}
      
      if(fabs(genjet_eta->at(ijet))>etaMinCore and fabs(genjet_eta->at(ijet))<etaMaxCore){
	if(genjet_eta->at(ijet)<0 and genjetpt_negEta>-1. and totTCpt_negEta_core_iter1>0.){
	  h2TCvsGenPt_CoreEta->Fill(genjetpt_negEta, totTCpt_negEta_core_iter1);
	  pTCvsGenPt_CoreEta->Fill(genjetpt_negEta, totTCpt_negEta_core_iter1);
	  h2TCvsGenPt_CEE_CoreEta->Fill(genjetpt_negEta, totTCpt_negEta_CEE_core_iter1);
	  pTCvsGenPt_CEE_CoreEta->Fill(genjetpt_negEta, totTCpt_negEta_CEE_core_iter1);
	}else if(genjet_eta->at(ijet)>0 and genjetpt_posEta>-1. and totTCpt_posEta_core_iter1>0.){
	  h2TCvsGenPt_CoreEta->Fill(genjetpt_posEta, totTCpt_posEta_core_iter1);
	  pTCvsGenPt_CoreEta->Fill(genjetpt_posEta, totTCpt_posEta_core_iter1);
	  h2TCvsGenPt_CEE_CoreEta->Fill(genjetpt_posEta, totTCpt_posEta_CEE_core_iter1);
	  pTCvsGenPt_CEE_CoreEta->Fill(genjetpt_posEta, totTCpt_posEta_CEE_core_iter1);
	}
      }

      if(fabs(genjet_eta->at(ijet))>etaMinCore and fabs(genjet_eta->at(ijet))<etaMaxCore){
	if(genjet_eta->at(ijet)<0 and genjetpt_negEta>-1. and totTCpt_negEta_core_rozcut0p15_iter1>0.){
	  h2TCvsGenPt_CoreEta_rozcut0p15->Fill(genjetpt_negEta, totTCpt_negEta_core_rozcut0p15_iter1);
	  pTCvsGenPt_CoreEta_rozcut0p15->Fill(genjetpt_negEta, totTCpt_negEta_core_rozcut0p15_iter1);
	}else if(genjet_eta->at(ijet)>0 and genjetpt_posEta>-1. and totTCpt_posEta_core_rozcut0p15_iter1>0.){
	  h2TCvsGenPt_CoreEta_rozcut0p15->Fill(genjetpt_posEta, totTCpt_posEta_core_rozcut0p15_iter1);
	  pTCvsGenPt_CoreEta_rozcut0p15->Fill(genjetpt_posEta, totTCpt_posEta_core_rozcut0p15_iter1);
	}
      }
      
      if(fabs(genjet_eta->at(ijet))>etaMinCore and fabs(genjet_eta->at(ijet))<etaMaxCore){
	if(genjet_eta->at(ijet)<0 and genjetpt_negEta>-1. and totTCpt_negEta_core_mpptcut2_iter1>0.){
	  h2TCvsGenPt_CoreEta_mpptcut2->Fill(genjetpt_negEta, totTCpt_negEta_core_mpptcut2_iter1);
	  pTCvsGenPt_CoreEta_mpptcut2->Fill(genjetpt_negEta, totTCpt_negEta_core_mpptcut2_iter1);
	}else if(genjet_eta->at(ijet)>0 and genjetpt_posEta>-1. and totTCpt_posEta_core_mpptcut2_iter1>0.){
	  h2TCvsGenPt_CoreEta_mpptcut2->Fill(genjetpt_posEta, totTCpt_posEta_core_mpptcut2_iter1);
	  pTCvsGenPt_CoreEta_mpptcut2->Fill(genjetpt_posEta, totTCpt_posEta_core_mpptcut2_iter1);
	}
      }
      
      
    }
    
    
    ///////////////////=========== Emulation ============= ///////////////////////
    //std::vector<TPGClusterData> vCld[6];
    std::vector<TPGCluster> vCld[6];    
    for (uint32_t isect = 0 ; isect < 6 ; isect++ ) {
      s2Clustering.run(vTcw[isect],vCld[isect]);
      if(doPrint) std::cout << isect << ", Size of Tcs: " << vTcw[isect].size() << ", Size of Clusters: " << vCld[isect].size() << std::endl;
    }
    ///////////////////=========== Emulation ============= ///////////////////////

    int nofClus1GeVp = 0, nofClus1GeVm = 0;
    int nofClus3GeVp = 0, nofClus3GeVm = 0;
    int nofClus5GeVp = 0, nofClus5GeVm = 0;
    int nofClus10GeVp = 0, nofClus10GeVm = 0;
    for (uint32_t isect = 0 ; isect < 6 ; isect++ ){
      //for(TPGCluster const& clf : vCld[isect]){
      if(isect==0){
	nofClus1GeVm = 0;
	nofClus3GeVm = 0;
	nofClus5GeVm = 0;
	nofClus10GeVm = 0;
      }
      if(isect==3){
	nofClus1GeVp = 0;
	nofClus3GeVp = 0;
	nofClus5GeVp = 0;
	nofClus10GeVp = 0;
      }
      int iclus = 0;
      for(TPGCluster const& clf : vCld[isect]){
	double tcPtSum = (clf.getGlobalEtaRad(isect) < 0) ? totTCpt_negEta : totTCpt_posEta ;
	if(doPrint and clf.getEnergyGeV()>3.0){
	  std::cout << "iclus-ievent: " << std::fixed << std::setprecision(2) << std::setw(4) << ievent
		    << ", sector: " << std::setw(5) << isect
		    << ", iclus: " << std::setw(5) << iclus++
		    <<", (x,y,z): (" << std::setw(7) << (clf.getGlobalXOverZF(isect) * clf.getZCm()) << ", " << std::setw(7) << (clf.getGlobalYOverZF(isect) * clf.getZCm()) << ", " << std::setw(7) << clf.getGlobalZCm(isect) << ")"
		    <<", (xoz,yoz,roz): (" << std::setw(5) << clf.getGlobalXOverZF(isect) << ", " << std::setw(5) << clf.getGlobalYOverZF(isect) << ", " << std::setw(5) << clf.getGlobalRhoOverZF(isect) << ")"
		    <<", (pt,eta,phi): (" << std::fixed << std::setprecision(2) << std::setw(8) << clf.getEnergyGeV()
		    << ", " << std::setw(8) << clf.getGlobalEtaRad(isect)
		    << ", " << std::setw(8) << (TMath::RadToDeg()* clf.getGlobalPhiRad(isect))
	            //<< ", " << std::setw(8) << tc_energy->at(itc)*1.e6
		    << ") "
		    // <<", (tcPtSum,tcPtSum/lsb,energy,energy*lsb): (" << std::setw(8) << tcPtSum
		    // << ", " << std::setw(8) << uint32_t(tcPtSum/lsbScales.LSB_E_TC()+0.5)
		    // << ", " << std::setw(8) << clf.getEnergy()
		    // << ", " << std::setw(8) << (clf.getEnergy()*lsbScales.LSB_E_TC())
		    // << ") "
		    << std::defaultfloat
		    << std::endl;
	}

	if(isect<=2){
	  if(clf.getEnergyGeV()>1.0) nofClus1GeVm++;
	  if(clf.getEnergyGeV()>3.0) nofClus3GeVm++;	
	  if(clf.getEnergyGeV()>5.0) nofClus5GeVm++;
	  if(clf.getEnergyGeV()>10.0) nofClus10GeVm++;
	}else{
	  if(clf.getEnergyGeV()>1.0) nofClus1GeVp++;
	  if(clf.getEnergyGeV()>3.0) nofClus3GeVp++;	
	  if(clf.getEnergyGeV()>5.0) nofClus5GeVp++;
	  if(clf.getEnergyGeV()>10.0) nofClus10GeVp++;
	}
      }//cluster loop
    }//sector loop
    if(nofClus1GeVm>0) hNClus1GeV->Fill(nofClus1GeVm);
    if(nofClus1GeVp>0) hNClus1GeV->Fill(nofClus1GeVp);
    if(nofClus3GeVm>0) hNClus3GeV->Fill(nofClus3GeVm);
    if(nofClus3GeVp>0) hNClus3GeV->Fill(nofClus3GeVp);
    if(nofClus5GeVm>0) hNClus5GeV->Fill(nofClus5GeVm);
    if(nofClus5GeVp>0) hNClus5GeV->Fill(nofClus5GeVp);
    if(nofClus10GeVm>0) hNClus10GeV->Fill(nofClus10GeVm);
    if(nofClus10GeVp>0) hNClus10GeV->Fill(nofClus10GeVp);
    
    for(int ipjet=0; ipjet < jetlist.size() ; ipjet++ ){
      int ijet = jetlist.at(ipjet).index ;
      //hEt->Fill(genjet_pt->at(ijet));
      double genjetpt = (genjet_eta->at(ijet) < 0) ? genjetpt_negEta : genjetpt_posEta;
      if(genjetpt==-1.) continue;      
      int minsect = (genjet_eta->at(ijet) < 0) ? 0 : 3;
      int maxsect = (genjet_eta->at(ijet) < 0) ? 3 : 6;
      double tcPtSum = (genjet_eta->at(ijet) < 0) ? totTCpt_negEta : totTCpt_posEta ;
      double tcPtSum_core = (genjet_eta->at(ijet) < 0) ? totTCpt_negEta_core_iter1 : totTCpt_posEta_core_iter1 ;
      double tcPtSum_core_tcroz = (genjet_eta->at(ijet) < 0) ? totTCpt_negEta_core_rozcut0p15_iter1 : totTCpt_posEta_core_rozcut0p15_iter1 ;
      double tcPtSum_CEE = (genjet_eta->at(ijet) < 0) ? totTCpt_negEta_CEE : totTCpt_posEta_CEE ;
      int tcDetType = (genjet_eta->at(ijet) < 0) ? maxTCPt_negEta_detType : maxTCPt_posEta_detType ;
      double tcPt = (genjet_eta->at(ijet) < 0) ? maxTCPt_negEta : maxTCPt_posEta ;
      double avgXoZ = (genjet_eta->at(ijet) < 0) ? tcXoZnegEta/tcPtSum : tcXoZposEta/tcPtSum ;
      double avgYoZ = (genjet_eta->at(ijet) < 0) ? tcYoZnegEta/tcPtSum : tcYoZposEta/tcPtSum ;
      double avgPhi = (genjet_eta->at(ijet) < 0) ? tcPhi_negEta/tcPtSum : tcPhi_posEta/tcPtSum ;
      double avgEta = (genjet_eta->at(ijet) < 0) ? tcEta_negEta/tcPtSum : tcEta_posEta/tcPtSum ;
      double gjroz = 1/sinh(genjet_eta->at(ijet)) ; //tan(2*atan(exp(-1.0*genjet_eta->at(ijet))));
      double gjxoz = gjroz*cos(genjet_phi->at(ijet));
      double gjyoz = gjroz*sin(genjet_phi->at(ijet));
      double avgRoZ = sqrt((avgXoZ-gjxoz)*(avgXoZ-gjxoz) + (avgYoZ-gjyoz)*(avgYoZ-gjyoz));
      bool hasFound = false, hasFoundPt = false;
      int binPt, binEta;
      int binPtClus, binEtaClus;
      int binPt_cl3d, binEta_cl3d;
      double cluspt = genjet_pt->at(ijet);
      int nofClus = 0;
      //hTrigGen->Fill(genjet_pt->at(ijet));
      double minRoz = 1., closestClusE= 100., closestClusEDiff= 100., totClusE = 0., closestGenjetByClus = 0.0, closestClusEta = 0.0;
      for (uint32_t isect = minsect ; isect < maxsect ; isect++ ){
	for(TPGCluster const& clf : vCld[isect]){	  
	  if(clf.getEnergyGeV()<3.0) continue;
	  double dClusGenRoz = sqrt( (clf.getGlobalXOverZF(isect) - gjxoz)*(clf.getGlobalXOverZF(isect) - gjxoz) + (clf.getGlobalYOverZF(isect) - gjyoz)*(clf.getGlobalYOverZF(isect) - gjyoz) );
	  if(dClusGenRoz<minRoz){
	    minRoz = dClusGenRoz;
	    closestClusE = clf.getEnergyGeV();
	    closestClusEDiff = clf.getEnergyGeV() - genjet_pt->at(ijet);
	    closestGenjetByClus = genjet_pt->at(ijet)/clf.getEnergyGeV() ;
	    closestClusEta = clf.getGlobalEtaRad(isect); 
	  }//minRoZ condition
	  totClusE += clf.getEnergyGeV() ; 
	}//cluster loop
      }//sector loop
      //closestClusE *= fClusCorr->Eval(closestClusE);
      
      double minRoz_cl3d = 1., closestClusE_cl3d= 0., closestClusEta_cl3d= 0.;
      for(int icl3=0;icl3<cl3d_n;icl3++){	  
	if(cl3d_pt->at(icl3)<3.0) continue;
	double clusXoZ = cos(cl3d_phi->at(icl3))/sinh(cl3d_eta->at(icl3)) ; 
	double clusYoZ = sin(cl3d_phi->at(icl3))/sinh(cl3d_eta->at(icl3)); 
	double dClusGenRoz = sqrt( (clusXoZ - gjxoz)*(clusXoZ - gjxoz) + (clusYoZ - gjyoz)*(clusYoZ - gjyoz) );
	if(dClusGenRoz<minRoz_cl3d){
	  minRoz_cl3d = dClusGenRoz;
	  closestClusE_cl3d = cl3d_pt->at(icl3);
	  closestClusEta_cl3d = cl3d_eta->at(icl3);
	}//minRoZ condition
      }//cluster loop
      //closestClusE_cl3d *= fcl3d->Eval(closestClusE_cl3d);
      
      bool hasClosestFound = (minRoz<1.)?true:false;
      bool hasClosestFound_cl3d = (minRoz_cl3d<1.)?true:false;
      // bool hasClosestFoundPt = (hasClosestFound and closestClusE>pt_clusThresh)?true:false;
      // bool hasClosestFoundPtTot = (hasClosestFound and totClusE>pt_clusThresh)?true:false;
      //Eta : outermost(Sci) = 1.321, innermost(Si) =  3.152
      bool hgcalEta = (fabs(genjet_eta->at(ijet))>1.321 and fabs(genjet_eta->at(ijet))<3.152) ? true : false;
      bool hgcalInnerEta = (fabs(genjet_eta->at(ijet))>etaMinCore and fabs(genjet_eta->at(ijet))<etaMaxCore) ? true : false;
      // bool hgcalOuterEta = (fabs(genjet_eta->at(ijet))>etaMinCore and fabs(genjet_eta->at(ijet))<etaMaxCore5) ? true : false;
      
      clusRoz->Fill(minRoz);
      clusRozPtwt->Fill(minRoz, genjet_pt->at(ijet));
      hJetEtaPtBin->Fill(genjet_eta->at(ijet), genjet_pt->at(ijet));
      GetHistoBin(hJetEtaPtBin,genjet_eta->at(ijet), genjet_pt->at(ijet), binEta, binPt);
      GetHistoBin(hJetEtaPtBin,closestClusEta, closestClusE, binEtaClus, binPtClus);
      GetHistoBin(hJetEtaPtBin,closestClusEta_cl3d, closestClusE_cl3d, binEta_cl3d, binPt_cl3d);
      
      if(hgcalEta and hasClosestFound) {
	hJetEtaPtBinDet->Fill( fabs(genjet_eta->at(ijet)), genjet_pt->at(ijet));
	
	hGenClusPtReso[binEtaClus-1][binPtClus-1]->Fill( genjet_pt->at(ijet)/closestClusE );
	hGenTCPtReso[binEtaClus-1][binPtClus-1]->Fill( genjet_pt->at(ijet)/tcPtSum );
	hTCClusPtReso[binEtaClus-1][binPtClus-1]->Fill( tcPtSum/closestClusE );
	
	hGenClusPtReso3D[tcDetType-1][binEtaClus-1][binPtClus-1]->Fill( genjet_pt->at(ijet)/closestClusE );
	
	if(tcDetType<=3){
	  hTCElossCEEID[binEtaClus-1]->Fill( tcPtSum_CEE );
	  hTCElossCE_E[binEtaClus-1][binPtClus-1]->Fill( tcPtSum_CEE );
	  // hTCElossCE_E_3D[tcDetType-1][binEtaClus-1][binPtClus-1]->Fill( tcPtSum_CEE );
	}else{
	  hTCElossCEHID[binEtaClus-1]->Fill( (tcPtSum-tcPtSum_CEE) );
	  hTCElossCE_H[binEtaClus-1][binPtClus-1]->Fill( (tcPtSum-tcPtSum_CEE) );
	  // hTCElossCE_H_3D[tcDetType-1][binEtaClus-1][binPtClus-1]->Fill( (tcPtSum-tcPtSum_CEE) );

	  //hHadGenClusPtReso2D[binEtaClus-1][binPtClus-1]->Fill( ((genjet_pt->at(ijet)-tcPtSum_CEE) - closestClusE) );
	  // hHadGenClusPtReso3D[tcDetType-1][binEtaClus-1][binPtClus-1]->Fill( ((genjet_pt->at(ijet)-tcPtSum_CEE) - closestClusE) );
	}

	h2ClusVsTCPt->Fill(tcPtSum, closestClusE);
	h2TotClusVsTCPt->Fill(tcPtSum, totClusE);
	if(totClusE>tcPtSum){
	  std::cerr<<"Issue: totClusE>tcPtSum in file : " << fin->GetName() << ", Event: " << ievent << ", ijet: " << ijet << ", genjetpt: " << genjetpt << ", genjeteta: " << genjet_eta->at(ijet) << ", phi: " << genjet_phi->at(ijet) << ", tcPtSum: " << tcPtSum << ", totClusE:" << totClusE  << std::endl;
	}
	if(hgcalInnerEta){
	  
	  hGenClusPtReso1D[binPtClus-1]->Fill( genjet_pt->at(ijet)/closestClusE );
	  hGenTCPtReso1D[binPtClus-1]->Fill( genjet_pt->at(ijet)/tcPtSum );
	  hTCClusPtReso1D[binPtClus-1]->Fill( tcPtSum/closestClusE );	  
	  
	  h2GenVsClusPt->Fill(closestClusE, genjet_pt->at(ijet));
	  h2GenVsTCPt->Fill(tcPtSum, genjet_pt->at(ijet));
	  h2TCVsClusPt->Fill( closestClusE, tcPtSum);
	  
	  h2ClusVsTCPt_CoreEta->Fill(tcPtSum_core, closestClusE);
	  h2ClusVsTCPt_CoreEta_rozcut0p15->Fill(tcPtSum_core_tcroz, closestClusE);

	  h2TotClusVsTCPt_CoreEta->Fill(tcPtSum_core, totClusE);
	  h2TotClusVsTCPt_CoreEta_rozcut0p15->Fill(tcPtSum_core_tcroz, totClusE);

	  // GetHistoBin(hJetEtaTCPtBin,genjet_eta->at(ijet), tcPt, binEta, binPt);
	  // hGenTCPtResovsTCPt1D[binPt-1]->Fill( (genjet_pt->at(ijet) - tcPtSum) );
	  // hGenTCPtResovsTCPt2D[tcDetType-1][binPt-1]->Fill( (genjet_pt->at(ijet) - tcPtSum) );
	  
	  if(hasClosestFound_cl3d){
	    hNewOldPtReso1Dcl3d[binPtClus-1]->Fill( closestClusE/closestClusE_cl3d );
	    h2NewVsOldPtcl3d->Fill( closestClusE, closestClusE_cl3d);
	    pNewVsOldPtcl3d->Fill( closestClusE, closestClusE_cl3d);
	  }

	}
      }
      if(hgcalInnerEta and hasClosestFound_cl3d) {
	hGenClusPtReso1Dcl3d[binPt_cl3d-1]->Fill( genjet_pt->at(ijet)/closestClusE_cl3d );
	h2GenVsClusPtcl3d->Fill(closestClusE_cl3d, genjet_pt->at(ijet));
      }
      
    }//jet loop 
    
    
    
    for(uint16_t i=0;i<6;i++) {
      vCld[i].clear();
      vTcw[i].clear();
    }
    
    taudlist.clear();
    taugdlist.clear();
    genlist.clear();
    partlist.clear();
    jetlist.clear();

    gen_pt->clear();
    gen_eta->clear();
    gen_phi->clear();
    gen_pdgid->clear();
    gen_energy->clear();
    gen_charge->clear();
    gen_status->clear();
    for(uint32_t i=0;i<gen_daughters->size();i++) (gen_daughters->at(i)).clear();
    gen_daughters->clear();
    genpart_pt->clear();
    genpart_eta->clear();
    genpart_phi->clear();
    genpart_pid->clear();
    genpart_mother->clear();
    genpart_energy->clear();
    genpart_gen->clear();
    genpart_fromBeamPipe->clear();
    genpart_exeta->clear();
    genpart_exphi->clear();
    genjet_pt->clear();
    genjet_eta->clear();
    genjet_phi->clear();
    genjet_energy->clear();
    tc_pt->clear();
    tc_mipPt->clear();
    tc_eta->clear();
    tc_phi->clear();
    tc_layer->clear();
    tc_x->clear();
    tc_y->clear();
    tc_z->clear();
    tc_energy->clear();
    tc_data->clear();
    tc_subdet->clear();
    tc_wafertype->clear();
    tc_uncompressedCharge->clear();
    cl3d_pt->clear();
    cl3d_phi->clear();
    cl3d_eta->clear();
  }//event loop
  
  std::string outname = "CalcResolution_" + outputfile_extn + "_" + index;

  TH1D *tcRoz_cumul = (TH1D *)tcRoz->GetCumulative(kTRUE,"_cumul");
  tcRoz_cumul->Scale(1./tcRoz_cumul->GetBinContent(tcRoz_cumul->GetMaximumBin()));

  TH1D *tcRozPtwt_cumul = (TH1D *)tcRozPtwt->GetCumulative(kTRUE,"_cumul");
  tcRozPtwt_cumul->Scale(1./tcRozPtwt_cumul->GetBinContent(tcRozPtwt_cumul->GetMaximumBin()));

  TH1D *clusRoz_cumul = (TH1D *)clusRoz->GetCumulative(kTRUE,"_cumul");
  clusRoz_cumul->Scale(1./clusRoz_cumul->GetBinContent(clusRoz_cumul->GetMaximumBin()));

  TFile *fout = new TFile(Form("%s.root",outname.c_str()),"recreate");
  hGenvsJetPtDiff->Write();
  pGenvsJetPtDiff->Write();
  hGenvsJetDR->Write();
  pPtDiffvsDR->Write();
  h2PtDiffvsDR->Write();
  hJetEtaPtBin->Write();
  hJetEtaPtBinDet->Write();
  h2TCPtvsGenPtByTCPtsum->Write();
  pTCPtvsGenPtByTCPtsum->Write();
  pTCPtvsGenPtByTCPtsumLogx->Write();
  tcRoz->Write();
  htcMipPt->Write();
  tcRozCore->Write();
  tcRoz_cumul->Write();
  tcRozPtwt->Write();
  tcRozPtwt_cumul->Write();
  clusRoz->Write();
  clusRoz_cumul->Write();
  clusRozPtwt->Write();
  hNJets->Write();
  hNJetsF1->Write();
  hNClus1GeV->Write();
  hNClus3GeV->Write();
  hNClus5GeV->Write();
  hNClus10GeV->Write();
  //=============================================
  TDirectory *dir = fout->mkdir("TC_Scale");
  dir->cd();
  //=
  h2TCvsGenPt->Write();
  h2TCvsGenPt_CEE->Write();
  h2TCvsGenPt_CEE_CoreEta->Write();
  pTCvsGenPt->Write();
  pTCvsGenPt_CEE->Write();
  pTCvsGenPt_CEE_CoreEta->Write();  
  //=
  h2TCvsGenPt_CoreEta->Write();
  pTCvsGenPt_CoreEta->Write();
  h2TCvsGenPt_CoreEta_rozcut0p15->Write();
  pTCvsGenPt_CoreEta_rozcut0p15->Write();
  h2TCvsGenPt_CoreEta_mpptcut2->Write();
  pTCvsGenPt_CoreEta_mpptcut2->Write();
  //=
  h2GenjetTCPtsumvsGenPt->Write();
  h2GenjetTCPtsumvsGenEta->Write();
  pGenjetTCPtsumvsGenPt->Write();
  pGenjetTCPtsumvsGenEta->Write();
  //=
  h2GenjetTCPtsumvsGenPt_CEE->Write();
  h2GenjetTCPtsumvsGenEta_CEE->Write();
  h2GenjetTCPtsumvsGenPt_CEE_CoreEta->Write();
  pGenjetTCPtsumvsGenPt_CEE->Write();
  pGenjetTCPtsumvsGenEta_CEE->Write();
  pGenjetTCPtsumvsGenPt_CEE_CoreEta->Write();
  //=
  h2GenjetTCPtsumvsTCPtsum_CEE_CoreEta->Write();
  pGenjetTCPtsumvsTCPtsum_CEE_CoreEta->Write();
  //=
  h2GenjetTCPtsumvsGenPt_CoreEta->Write();
  pGenjetTCPtsumvsGenPt_CoreEta->Write();
  //=
  fout->cd();
  //=============================================
  dir = fout->mkdir("TC_histos");
  dir->cd();
  hTCPt->Write();
  hTCEta->Write();
  hTCPhiPos->Write();
  hTCPhiNeg->Write();
  h2TCEtaPhi->Write();
  fout->cd();
  //=============================================
  dir = fout->mkdir("Clus");
  dir->cd();
  h2GenVsClusPt->Write();
  h2GenVsTCPt->Write();
  h2TCVsClusPt->Write();
  //=
  h2GenVsClusPtcl3d->Write();  
  h2NewVsOldPtcl3d->Write();  
  pNewVsOldPtcl3d->Write();
  //=
  h2ClusVsTCPt->Write();  
  h2ClusVsTCPt_CoreEta->Write();  
  h2ClusVsTCPt_CoreEta_rozcut0p15->Write();
  h2TotClusVsTCPt->Write();  
  h2TotClusVsTCPt_CoreEta->Write();  
  h2TotClusVsTCPt_CoreEta_rozcut0p15->Write();  
  //=
  fout->cd();
  //=============================================
  dir = fout->mkdir("Reso_1D");
  dir->cd();
  for(int ipt=0;ipt<nJetPtBins;ipt++) hGenClusPtReso1D[ipt]->Write();
  for(int ipt=0;ipt<nJetPtBins;ipt++) hGenTCPtReso1D[ipt]->Write();
  for(int ipt=0;ipt<nJetPtBins;ipt++) hTCClusPtReso1D[ipt]->Write();
  // for(int ipt=0;ipt<nJetPtBins;ipt++) hHadGenClusPtReso1D[ipt]->Write();
  // for(int ipt=0;ipt<nJetPtBins;ipt++) hHadGenTCPtReso1D[ipt]->Write();
  // for(int ipt=0;ipt<nJetPtBins;ipt++) hHadTCClusPtReso1D[ipt]->Write();
  for(int ipt=0;ipt<nJetPtBins;ipt++) hGenClusPtReso1Dcl3d[ipt]->Write();
  for(int ipt=0;ipt<nJetPtBins;ipt++) hNewOldPtReso1Dcl3d[ipt]->Write();
  for(int ieta=0;ieta<nJetEtaBins;ieta++) hTCElossCEEID[ieta]->Write();
  for(int ieta=0;ieta<nJetEtaBins;ieta++) hTCElossCEHID[ieta]->Write();
  for(int ipt=0;ipt<int(logbins.size()-1);ipt++) hGenjetByTCPtSumVsTCPt1DLog[ipt]->Write();
  for(int idet=0;idet<7;idet++) pTCPtvsGenPtByTCPtsumDet[idet]->Write();
  for(int idet=0;idet<7;idet++) tcRozDet[idet]->Write();
  fout->cd();
  //=============================================
  dir = fout->mkdir("Reso_2D");
  dir->cd();
  for(int ieta=0;ieta<nJetEtaBins;ieta++){
    for(int ipt=0;ipt<nJetPtBins;ipt++){
      hGenClusPtReso[ieta][ipt]->Write();
    }// jet et
  }//jet eta
  
  for(int ieta=0;ieta<nJetEtaBins;ieta++){
    for(int ipt=0;ipt<nJetPtBins;ipt++){
      hTCClusPtReso[ieta][ipt]->Write();
    }// jet et
  }//jet eta
  
  for(int ieta=0;ieta<nJetEtaBins;ieta++){
    for(int ipt=0;ipt<nJetPtBins;ipt++){
      hGenTCPtReso[ieta][ipt]->Write();
    }// jet et
  }//jet eta
  
  for(int ieta=0;ieta<nJetEtaBins;ieta++){
    for(int ipt=0;ipt<nJetPtBins;ipt++){
      hTCElossCE_E[ieta][ipt]->Write();
    }// jet et
  }//jet eta
  
  for(int ieta=0;ieta<nJetEtaBins;ieta++){
    for(int ipt=0;ipt<nJetPtBins;ipt++){
      hTCElossCE_H[ieta][ipt]->Write();
    }// jet et
  }//jet eta

  for(int idet=0;idet<7;idet++)
    for(int ipt=0;ipt<int(logbins.size()-1);ipt++)
      hGenjetByTCPtSumVsTCPt2DDetLog[idet][ipt]->Write();

  for(int ieta=0;ieta<nJetEtaBins;ieta++)
    for(int ipt=0;ipt<int(logbins.size()-1);ipt++)
      hGenjetByTCPtSumVsTCPt2DLog[ieta][ipt]->Write();
  
  fout->cd();
  //=============================================
  dir = fout->mkdir("Reso_3D");
  dir->cd();
  for(int idet=0;idet<7;idet++){
    for(int ieta=0;ieta<nJetEtaBins;ieta++){
      for(int ipt=0;ipt<nJetPtBins;ipt++){
	hGenClusPtReso3D[idet][ieta][ipt]->Write();
      }// jet et
    }//jet eta
  }//idet


  fout->cd();
  // //=============================================

  // //=============================================
  // dir = fout->mkdir("Reso_3D_TC");
  // dir->cd();
  // for(int idet=0;idet<7;idet++){
  //   for(int ieta=0;ieta<nJetEtaBins;ieta++){
  //     for(int ipt=0;ipt<int(logbins.size()-1);ipt++){
  // 	hGenjetByTCPtSumVsTCPt3DLog[idet][ieta][ipt]->Write();
  //     }// jet et
  //   }//jet eta
  // }//idet
  // fout->cd();
  // //=============================================
  // dir = fout->mkdir("CEE_3D");
  // dir->cd();
  // for(int idet=0;idet<7;idet++){
  //   for(int ieta=0;ieta<nJetEtaBins;ieta++){
  //     for(int ipt=0;ipt<nJetPtBins;ipt++){
  // 	hTCElossCE_E_3D[idet][ieta][ipt]->Write();
  //     }// jet et
  //   }//jet eta
  // }//idet
  
  // fout->cd();
  // dir = fout->mkdir("CEH_3D");
  // dir->cd();
  // for(int idet=0;idet<7;idet++){
  //   for(int ieta=0;ieta<nJetEtaBins;ieta++){
  //     for(int ipt=0;ipt<nJetPtBins;ipt++){
  // 	hTCElossCE_H_3D[idet][ieta][ipt]->Write();
  //     }// jet et
  //   }//jet eta
  // }//idet  
  // fout->cd();
  // //=============================================
  
  fout->Close();
  delete fout;

  logbins.clear();
  // fin->Close();
  // delete fin;

  return true;
}
