/**********************************************************************
 Created on : 24/04/2025
 Purpose    : Stage2 semi-emulator test
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
#include <TGraph.h>

#include "TParticlePDG.h"
#include "TDatabasePDG.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

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
  }else if(sampleType.find("singlePion")!=std::string::npos){
    pt_clusThresh = 100.; 
    pt_clusThresh_etaeff = 200.; 
    // pt_clusThresh = 150.; 
    // pt_clusThresh_etaeff = 200.; 
  }else if(sampleType.find("minbias")!=std::string::npos){
    pt_clusThresh = 150.; 
    pt_clusThresh_etaeff = 200.; 
  }
  
  ////////////////// ============= cuts ===========///////////////////////
  float etaMin = 1.321;
  float etaMax = 3.152;
  
  float innerEtaMin = 1.7;
  float innerEtaMax = 2.9;
  
  double tcRoZCut = 0.15;
  double delRTh = 0.25;
  ////////////////// ============= cuts ===========///////////////////////
  
  
  std::unique_ptr<TFile> fin( TFile::Open(Form(inputfile.c_str())) );
  //TFile *fin = TFile::Open(Form(inputfile.c_str()));
  std::cout<<"Input: filename : " << fin->GetName() << std::endl;
  std::cout << "Input: outfile index : " << index << std::endl;
  std::cout << "Input: nofEvents : " << nofEvents << std::endl;
  std::cout << "Input: sidelength : " << sidelength << std::endl;
  std::cout << "Input: MC Sample : " << sampleType << std::endl;
  std::cout << "pt_clusThresh for L1 : " << pt_clusThresh << std::endl;
  std::cout << "pt_clusThresh_etaeff for Trigger eta eff : " << pt_clusThresh_etaeff << std::endl;
  int  resoIndex = int(1000*sidelength);
  std::cout << "resoIndex : " << resoIndex << std::endl;
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

  // Int_t cl3d_n = 0 ;
  // tr->SetBranchStatus("cl3d_n",1);
  // tr->SetBranchAddress("cl3d_n" , &cl3d_n);

  // std::vector<float>  *cl3d_pt = 0 ;
  // tr->SetBranchStatus("cl3d_pt",1);
  // tr->SetBranchAddress("cl3d_pt" , &cl3d_pt);
  
  // std::vector<float>  *cl3d_eta = 0 ;
  // tr->SetBranchStatus("cl3d_eta",1);
  // tr->SetBranchAddress("cl3d_eta" , &cl3d_eta);
  
  // std::vector<float>  *cl3d_phi = 0 ;
  // tr->SetBranchStatus("cl3d_phi",1);
  // tr->SetBranchAddress("cl3d_phi" , &cl3d_phi);

  std::cout << "Total number of Events: " << tr->GetEntries() << std::endl;
  // std::cout << "Sizeof: " << sizeof(TPGTCBits) << std::endl;
  // std::cout << "Sizeof: " << sizeof(uint32_t) << std::endl;

  // TPGTCBits tcbit;
  // tcbit.setEnergy(262142);
  // tcbit.setROverZ(4094);
  // tcbit.setPhi(4094);
  // tcbit.setLayer(62);
  // tcbit.setU(32);
  // tcbit.setUbar(62);
  // tcbit.print();

  ///////////////========== Pt,eta hist bins =========///////////////
  
  const double effPtBin[29] = {0, 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 300., 400., 600.};
  
  const int nJetEtaBins = 6;
  Float_t jetEtaBin[nJetEtaBins+1] = {1.321, 1.7, 2.0, 2.3, 2.6, 2.9, 3.152} ;
  const int nJetPtBins = 43;
  Float_t jetPtBin[44] = {0, 15., 16., 18., 20., 22., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100.,  
                       110., 120., 130., 140., 150., 160., 170., 180., 190., 200.,
	               220., 240., 260., 280., 300.,
	               330., 360., 390., 420.,
                       440., 480., 520.};
  //////////////////////////////////////////////////////////////////
  
  // TF1* f2 = new TF1(Form("func"),"0.01+[0]/(x^[1]-[2])",10,1000);
  // f2->SetParameters(-0.00428022,-0.297381,0.6113);

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
  //fClusCorr->SetParameters(1.22466, 22.8516, 0.401581); //PU0 VBF
  fClusCorr->SetParameters(1.23, 16.95, 15.84); //PU0 pion size =16 pt 100
  //fClusCorr->SetParameters(1.0, 53.02, 13.91); //PU0 VBF size =16
  //fClusCorr->SetParameters(1.03, 24.92, 2.54); //PU0 VBF size =30
  //fClusCorr->SetParameters(1.06, 18.96, 4.28); //PU0 VBF size =45
  
  TF1* fcl3d = new TF1("fcl3d","[0]+[1]/(x+[2])",7,500);
  fcl3d->SetParameters(1.17576, 19.7599, 0.908577); //PU0 VBF

  //VBF correction 16
  float par0[nJetEtaBins] = {1.07956, 1.2916, 1.18022, 1.14248, 1.12945, 0.976157};
  float par1[nJetEtaBins] = {68.9515, 51.967, 52.3593, 47.0961, 39.7107, 32.1626};
  float par2[nJetEtaBins] = {18.5932, 5.89004, 8.34816, 10.0994, 9.27638, 16.5888};
  
  // //VBF correction 30
  // float par0[6] = {1.03452, 1.17763, 1.04363, 1.11672, 1.11673, 1.00139};
  // float par1[6] = {41.6436, 34.3258, 47.9796, 26.7711, 20.9924, 27.4734};
  // float par2[6] = {16.7377, 1.34016, 16.0725, 4.91101, 3.39069, 24.9721};
  
  // //VBF correction 45
  // float par0[6] = {0.9205, 1.09595, 1.01855, 1.09636, 1.06454, 0.91201};
  // float par1[6] = {53.1531, 28.7612, 35.9239, 18.509, 18.412, 24.4667};
  // float par2[6] = {30.0757, 4.80461, 13.5817, 1.28464, 6.60814, 20.6122};
  
  // //Poin correction (pt150) [size =16]
  // float par0[nJetEtaBins] = {1.07567, 1.17335, 1.18182, 1.19163, 1.17032, 1.1243};
  // float par1[nJetEtaBins] = {26.7441, 16.8533, 14.6125, 13.4602, 12.0038, 17.2314};
  // float par2[nJetEtaBins] = {11.3612, 17.8989, 15.6896, 14.1095, 13.8239, 7.32437};

  // //Poin correction (pt100) [size =16]
  // float par0[6] = {1.07567, 1.17335, 1.18182, 1.19163, 1.17032, 1.1243};
  // float par1[6] = {26.7441, 16.8533, 14.6125, 13.4602, 12.0038, 17.2314};
  // float par2[6] = {11.3612, 17.8989, 15.6896, 14.1095, 13.8239, 7.32437};

  TF1* fClusEtaCorr[nJetEtaBins];
  for(int ieta=0;ieta<nJetEtaBins;ieta++) {
    fClusEtaCorr[ieta] = new TF1(Form("fClusEtaCorr_%d",ieta),"[0]+[1]/(x+[2])",1.0,500);
    fClusEtaCorr[ieta]->SetParameters(par0[ieta], par1[ieta], par2[ieta]); 
  }
  //////////////////////////////////////////////////////////////
  ///////////==== pt,eta bin finder =================//////////
  TH2F *hJetEtaPtBin = new TH2F("hJetEtaPtBin","hJetEtaPtBin", nJetEtaBins, jetEtaBin, nJetPtBins, jetPtBin);
  //////////////////////////////////////////////////////////////
  
  ////////////////////////////////////========== Book histogram ===========//////////////////////////////////////////
  // TH1D *hEt = new TH1D("hEt","hEt",1000,0,1000.);
  // TH1D *hPhi = new TH1D("hPhi","hPhi",2*TMath::TwoPi(),-1.0*TMath::TwoPi(),TMath::TwoPi());
  // TH1D *hPhiDeg = new TH1D("hPhiDeg","hPhiDeg",2*TMath::RadToDeg()*TMath::TwoPi(),-1.0*TMath::RadToDeg()*TMath::TwoPi(),TMath::RadToDeg()*TMath::TwoPi());

  ////////////////////////////////////////////////////////////
  TH1D *hNJets = new TH1D("hNJets",Form("Nof VBF Jets"), 10,-0.5,9.5);
  TH1D *hNJetsF1 = new TH1D("hNJetsF1",Form("Nof VBF Jets (filtered for atmost 2)"), 10,-0.5,9.5);
  ////////////////////////////////////////////////////////////
  
  TH2D *h2GenVsClusPt = new TH2D("h2GenVsClusPt","h2GenVsClusPt", 100, 0., 500., 100, 0., 500.);
  TH2D *h2GenVsTCPt = new TH2D("h2GenVsTCPt","h2GenVsTCPt", 100, 0., 500., 100, 0., 500.);
  TH2D *h2TCVsClusPt = new TH2D("h2TCVsClusPt","h2TCVsClusPt", 100, 0., 500., 100, 0., 500.);
  TH2D *h2GenVsClusPtcl3d = new TH2D("h2GenVsClusPtcl3d","h2GenVsClusPtcl3d", 100, 0., 500., 100, 0., 500.);
  TH2D *h2NewVsOldPtcl3d = new TH2D("h2NewVsOldPtcl3d","h2NewVsOldPtcl3d", 100, 0., 500., 100, 0., 500.);
  
  TEfficiency* effTrigGen = new TEfficiency("effTrigGen","Trigger efficiency;p_{T}(GeV);#epsilon",600,0,600.);
  TEfficiency* effTrigGenEta = new TEfficiency("effTrigGenEta","Trigger efficiency;#eta;#epsilon",400,0,4.);
  TEfficiency* effTrigGenPhi = new TEfficiency("effTrigGenPhi","Trigger efficiency;#phi;#epsilon",50,-200.,200.);
  TEfficiency* effTrigGenTDR = new TEfficiency("effTrigGenTDR","Trigger efficiency;p_{T}(GeV);#epsilon",28,effPtBin);
  TEfficiency* effTrigGenTDR1D = new TEfficiency("effTrigGenTDR1D","Trigger efficiency;p_{T}(GeV);#epsilon",28,effPtBin);
  TEfficiency* effTrigGenTDR_UnCorr = new TEfficiency("effTrigGenTDR_UnCorr","Trigger efficiency;p_{T}(GeV);#epsilon",28,effPtBin);
  TEfficiency* effTrigGenTDR_cl3d = new TEfficiency("effTrigGenTDR_cl3d","Trigger efficiency;p_{T}(GeV);#epsilon",28,effPtBin);
  TEfficiency* effTrigGenTDR_cl3d_UnCorr = new TEfficiency("effTrigGenTDR_cl3d_UnCorr","Trigger efficiency;p_{T}(GeV);#epsilon",28,effPtBin);

  TEfficiency* rateSingles1D = new TEfficiency("rateSingles1D","Trigger rate;p_{T}(GeV);rate",600,0,600.);
  TEfficiency* rateSingles2D = new TEfficiency("rateSingles2D","Trigger rate;p_{T}(GeV);rate",600,0,600.);
  TEfficiency* rateSingles_Uncorr = new TEfficiency("rateSingles_Uncorr","Trigger rate;p_{T}(GeV);rate",600,0,600.);
  const int nRateBins = 200;
  TH1D *hBkgRateSingle1D = new TH1D("hBkgRateSingle1D","hBkgRateSingle1D", nRateBins, 0., float(nRateBins));
  TH1D *hBkgRateSingle2D = new TH1D("hBkgRateSingle2D","hBkgRateSingle2D", nRateBins, 0., float(nRateBins));
  TH1D *hBkgRateSingle_Uncorr = new TH1D("hBkgRateSingle_Uncorr","hBkgRateSingle_Uncorr", nRateBins, 0., float(nRateBins));
  
  TH1D *hBkgRateSingle1D_MaxPt = new TH1D("hBkgRateSingle1D_MaxPt","hBkgRateSingle1D_MaxPt", nRateBins, 0., float(nRateBins));
  TH1D *hBkgRateSingle2D_MaxPt = new TH1D("hBkgRateSingle2D_MaxPt","hBkgRateSingle2D_MaxPt", nRateBins, 0., float(nRateBins));  
  TH1D *hBkgRateSingle_Uncorr_MaxPt = new TH1D("hBkgRateSingle_Uncorr_MaxPt","hBkgRateSingle_Uncorr_MaxPt", nRateBins, 0., float(nRateBins));
  
  
  TH1D *hSelEvents = new TH1D("hSelEvents","hSelEvents", 10,-0.5,9.5);
  TH1D *hTotClus = new TH1D("hTotClus","hTotClus", 600,0,600.);
  TH1D *hClusPtNoCut = new TH1D("hClusPtNoCut","hClusPtNoCut", 600,0,600.);
  TH1D *hClusPt1D = new TH1D("hClusPt1D","hClusPt1D", 600,0,600.);
  TH1D *hClusPt2D = new TH1D("hClusPt2D","hClusPt2D", 600,0,600.);
  TH1D *hClusPt_Uncorr = new TH1D("hClusPt_Uncorr","hClusPt_Uncorr", 600,0,600.);
  ////////////////////////////////////========== Book histogram ===========//////////////////////////////////////////
  
  uint64_t totalEntries = tr->GetEntries();
  if(nofEvents==0 or nofEvents>totalEntries) nofEvents = totalEntries;
  
  std::cout << "Total Entries : " << totalEntries << std::endl;
  std::cout << "Loop for: " << nofEvents << " entries"  << std::endl;

  
  const auto default_precision{std::cout.precision()};
    
  TPGStage2Configuration::ClusPropLUT cplut;
  cplut.readMuEtaLUT("input/stage2/configuration/mean_eta_LUT.csv");
  cplut.readSigmaEtaLUT("input/stage2/configuration/sigma_eta_LUT.csv");

  TPGStage2Emulation::Stage2::_rOverZ = sidelength * sqrt(3.0) ;
  TPGStage2Emulation::Stage2 s2Clustering;
  s2Clustering.setClusPropLUT(&cplut);
  //s2Clustering.setConfiguration(&sb);
  s2Clustering.setkpower(0);
  std::cout << " Stage2::_rOverZ : " << s2Clustering.getROverZ() << std::endl;
  
  std::string outtreefname = "TPGReco_tree_" + outputfile_extn + "_" + index;
  TDirectory *savedir = gDirectory;
  //std::vector<TPGCluster> tpgClus;
  std::vector<float> clus_pt,clus_eta,clus_phi;
  std::vector<float> clus_pt_corr1D, clus_pt_corr2D;
  std::vector<unsigned int> clus_pass, clus_sect;
  std::vector<float> clus_local_eta,clus_local_phi;
  std::vector<float> tci_pt,tci_layer,tci_x,tci_y,tci_z,tci_sect;
  std::vector<float> tcf_pt,tcf_layer,tcf_xbyz,tcf_ybyz,tcf_zcm,tcf_sect;
  
  TFile *foutTree = new TFile(Form("%s.root",outtreefname.c_str()),"recreate");
  TTree *outputTree = new TTree("TPG_Reco","TPG_Reco");
  outputTree->SetAutoSave();
  //outputTree->Branch("tpgClus", &tpgClus);
  outputTree->Branch("clus_pass", &clus_pass);
  outputTree->Branch("clus_pt", &clus_pt);
  outputTree->Branch("clus_pt_corr1D", &clus_pt_corr1D);
  outputTree->Branch("clus_pt_corr2D", &clus_pt_corr2D);
  outputTree->Branch("clus_eta", &clus_eta);
  outputTree->Branch("clus_phi", &clus_phi);
  outputTree->Branch("clus_sect", &clus_sect);
  outputTree->Branch("clus_local_eta", &clus_local_eta);
  outputTree->Branch("clus_local_phi", &clus_local_phi);
  // outputTree->Branch("tci_pt", &tci_pt);
  // outputTree->Branch("tci_layer", &tci_layer);
  // outputTree->Branch("tci_x", &tci_x);
  // outputTree->Branch("tci_y", &tci_y);
  // outputTree->Branch("tci_z", &tci_z);
  // outputTree->Branch("tci_sect", &tci_sect);
  // outputTree->Branch("tcf_pt", &tcf_pt);
  // outputTree->Branch("tcf_layer", &tcf_layer);
  // outputTree->Branch("tcf_xbyz", &tcf_xbyz);
  // outputTree->Branch("tcf_ybyz", &tcf_ybyz);
  // outputTree->Branch("tcf_zcm", &tcf_zcm);
  // outputTree->Branch("tcf_sect", &tcf_sect);
  savedir->cd();
  
  //TPGTriggerCellFloats tcf0,tcf1;
  TPGTCFloats tcf0,tcf1;
  for (Long64_t ievent = 0 ; ievent < nofEvents ; ievent++ ) {
    
    tr->GetEntry(ievent) ;
    
    if(doPrint) std::cout<<"Event : "<< ievent <<", nof gen "<< gen_n << ", nof_gen: " << genpart_pt->size()  <<", TCs : "<< tc_pt->size() << std::endl;
    if(ievent%1==0) std::cout<<"Event : "<< ievent <<", nof TCs : "<< tc_pt->size() << std::endl;
    hSelEvents->Fill(0);
    /////////////////////////////============ Pythia Gen ===================////////////////////////////////
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

      if(doPrint
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
    /////////////////////////////============ Pythia Gen ===================////////////////////////////////
    
    
    ////////////////////////////================ GEANT genpart ===================////////////////////////////
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
    ////////////////////////////================ GEANT genpart ===================////////////////////////////
    
    
    /////////////////////////////////=========== Genjet =============//////////////////////////////////////////    
    std::vector<JetPart> jetlist;
    int nofHGCalJetsPos = 0, nofHGCalJetsNeg = 0;
    for(int ijet=0; ijet<genjet_n; ijet++ ){
      if(fabs(genjet_eta->at(ijet))<=etaMin or fabs(genjet_eta->at(ijet))>=etaMax) continue;
      double phiDeg = TMath::RadToDeg()*genjet_phi->at(ijet);
      if(genjet_eta->at(ijet)>0.){
	nofHGCalJetsPos++;
      }else{
	nofHGCalJetsNeg++;
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
	  
	  // if(sampleType.find("vbfHInv")!=std::string::npos) {
	  //   int mindr = TMath::LocMin(3,minDR);
	  //   int minId = mindr + 5; //index of non-Higgs particles of the hardest subprocess in VBF sample
	  //   hGenvsJetPtDiff->Fill(gen_pt->at(minId)-genjet_pt->at(ijet));	  
	  //   pGenvsJetPtDiff->Fill( gen_pt->at(minId),  gen_pt->at(minId)/genjet_pt->at(ijet) );
	  //   pPtDiffvsDR->Fill( minDR[mindr],  gen_pt->at(minId)/genjet_pt->at(ijet) );
	  //   h2PtDiffvsDR->Fill( minDR[mindr],  gen_pt->at(minId)/genjet_pt->at(ijet) );
	  //   hGenvsJetDR->Fill( minDR[mindr]  );
	  // }else{
	  //   hGenvsJetPtDiff->Fill(gen_pt->at(ipart)-genjet_pt->at(ijet));	  
	  //   pGenvsJetPtDiff->Fill( gen_pt->at(ipart),  gen_pt->at(ipart)/genjet_pt->at(ijet) );
	  //   pPtDiffvsDR->Fill( minDeltaR,  gen_pt->at(ipart)/genjet_pt->at(ijet) );
	  //   h2PtDiffvsDR->Fill( minDeltaR,  gen_pt->at(ipart)/genjet_pt->at(ijet) );
	  //   hGenvsJetDR->Fill( minDeltaR  );
	  // }
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
    /////////////////////////////////=========== Genjet =============//////////////////////////////////////////
    hNJets->Fill(float(jetlist.size()));
    hSelEvents->Fill(1);
    // ////////////////////////Conditions for two jets ////////////////////////////////////
    // bool isOneperSide = true ;
    // if(jetlist.size()==2){
    //   isOneperSide = ((genjet_eta->at(jetlist.at(0).index) * genjet_eta->at(jetlist.at(1).index))<0.)?true:false;
    //   // double *allJets = new double[jetlist.size()];
    //   // int *sortedIdx = new int[jetlist.size()];
    //   // for(int ipjet=0; ipjet < jetlist.size() ; ipjet++ ){
    //   // 	int ijet = jetlist.at(ipjet).index ;
    //   // 	allJets[ipjet] = genjet_pt->at(ijet);
    //   // }
    //   // TMath::Sort(int(jetlist.size()), allJets, sortedIdx);
    //   // int HptIdx = jetlist.at(sortedIdx[0]).index ;
    //   // int MptIdx = jetlist.at(sortedIdx[1]).index ;
    //   // double eta12 = genjet_eta->at(HptIdx) * genjet_eta->at(MptIdx) ;
    //   // std::vector<JetPart> tJetlist;
    //   // for(int ipjet=0; ipjet < jetlist.size() ; ipjet++ ){
    //   // 	if(eta12<0.){
    //   // 	  if(ipjet==sortedIdx[0] or ipjet==sortedIdx[1]){
    //   // 	    JetPart jet;
    //   // 	    jet.name = jetlist.at(ipjet).name ;
    //   // 	    jet.index = jetlist.at(ipjet).index ;
    //   // 	    tJetlist.push_back(jet);
    //   // 	  }
    //   // 	}else{
    //   // 	  if(ipjet==sortedIdx[0]){
    //   // 	    JetPart jet;
    //   // 	    jet.name = jetlist.at(ipjet).name ;
    //   // 	    jet.index = jetlist.at(ipjet).index ;
    //   // 	    tJetlist.push_back(jet);
    //   // 	  }	  
    //   // 	}
    //   // }
    //   // jetlist.clear();
    //   // for(int ipjet=0; ipjet < tJetlist.size() ; ipjet++ ){
    //   // 	int ijet = tJetlist.at(ipjet).index ;
    //   // 	JetPart jet;
    //   // 	jet.name = tJetlist.at(ipjet).name ;
    //   // 	jet.index = tJetlist.at(ipjet).index ;
    //   // 	jetlist.push_back(jet);
    //   // 	if(doPrint)
    //   // 	  std::cout << "filtered-ijet-ievent: " << std::setprecision(default_precision) << std::setw(4) <<ievent << ", ijet: " << std::setw(4) << ijet
    //   // 		    <<", Name: " << std::setw(10) << jet.name
    //   // 		    <<" jet:(pt,eta,phi,energy) : (" << std::fixed << std::setprecision(2) << std::setw(8) << genjet_pt->at(ijet)
    //   // 		    << ", " << std::setw(8) << genjet_eta->at(ijet) << ", " //<< std::setw(5) << genjet_exeta->at(ijet)
    //   // 		    << ", " << std::setw(8) << (TMath::RadToDeg()*genjet_phi->at(ijet)) //<< ", " << std::setw(5) << (TMath::RadToDeg()*genjet_exphi->at(ijet))
    //   // 		    << ", " << std::setw(8) << genjet_energy->at(ijet) << ") "
    //   // 		    << std::defaultfloat
    //   // 		    << std::endl;
    //   // }      
    //   // tJetlist.clear();      
    //   // delete allJets;
    //   // delete sortedIdx;
    // }
    // ////////////////////////////////////////////////////////////////////////////////////////
    
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
    
    // /////////////////////////// Conditions for 0 jets and others ////////////////////////////////////

    // if(njets==0 or (jetlist.size()==2 and !isOneperSide) or njets>2){
    //   //if(jetlist.size()==0 or !isOneperSide or nofHGCalJetsPos>1 or nofHGCalJetsNeg>1){
    //   //if(jetlist.size()==0 or jetlist.size()>2 or !isOneperSide){
    //   //if(jetlist.size()==0){
    //   taudlist.clear();
    //   taugdlist.clear();
    //   genlist.clear();
    //   partlist.clear();
    //   //jetlist.clear();
      
    //   gen_pt->clear();
    //   gen_eta->clear();
    //   gen_phi->clear();
    //   gen_pdgid->clear();
    //   gen_energy->clear();
    //   gen_charge->clear();
    //   gen_status->clear();
    //   for(uint32_t i=0;i<gen_daughters->size();i++) (gen_daughters->at(i)).clear();
    //   gen_daughters->clear();
    //   genpart_pt->clear();
    //   genpart_eta->clear();
    //   genpart_phi->clear();
    //   genpart_pid->clear();
    //   genpart_mother->clear();
    //   genpart_energy->clear();
    //   genpart_gen->clear();
    //   genpart_fromBeamPipe->clear();
    //   genpart_exeta->clear();
    //   genpart_exphi->clear();
    //   genjet_pt->clear();
    //   genjet_eta->clear();
    //   genjet_phi->clear();
    //   genjet_energy->clear();
    //   tc_pt->clear();
    //   //tc_mipPt->clear();
    //   tc_eta->clear();
    //   tc_phi->clear();
    //   tc_layer->clear();
    //   tc_x->clear();
    //   tc_y->clear();
    //   tc_z->clear();
    //   tc_energy->clear();
    //   tc_data->clear();
    //   tc_subdet->clear();
    //   tc_wafertype->clear();
    //   tc_uncompressedCharge->clear();
    //   cl3d_pt->clear();
    //   cl3d_phi->clear();
    //   cl3d_eta->clear();
    //   continue;
    // }
    // ////////////////////////////////////////////////////////////////////////////////////////
    
    hNJetsF1->Fill(float(njets));
    hSelEvents->Fill(2);
    ///////////////////========= Trigger Cell ==========/////////////////////////
    double totTCpt_posEta = 0.,totTCpt_negEta = 0.;
    double totTCpt_posEta_layer1 = 0.,totTCpt_negEta_layer1 = 0.;
    uint32_t nofPosEtaTCs = 0, nofNegEtaTCs = 0;
    double tcXoZposEta = 0.,tcYoZposEta = 0., tcXoZnegEta = 0.,tcYoZnegEta = 0.;
    double maxroz_posEta = -1.0, maxroz_negEta = -1.0;
    double tcPhi_posEta = 0.,tcPhi_negEta = 0.;
    double tcEta_posEta = 0.,tcEta_negEta = 0.;
    int detType = -1;
    int subdet = -1, wafertype = -1;

    tci_pt.clear();
    tci_layer.clear();
    tci_x.clear();
    tci_y.clear();
    tci_z.clear();
    tci_sect.clear();
    std::vector<TPGTCBits> vTcw[6];
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
      
      double tcptMeV = tc_pt->at(itc) * 1.e3;
      double tcpt = (tcptMeV<3.)?3:tcptMeV;
      //double ptcorr = f2->Eval(tcpt);
      double ptcorr = (tc_layer->at(itc)<=26)?1.07:1.14;//1.07:1.13f2->Eval(tcpt);
      tcpt = ptcorr*tc_pt->at(itc);
      tcptMeV = tcpt*1.e3;
            
      if(doPrint){
	std::cout << "itc-ievent: " << std::fixed << std::setprecision(2) << std::setw(4) << ievent
		  << ", itc: " << std::setw(5) << itc
		  << ", layer: " << std::setw(2) << int(tc_layer->at(itc))
		  <<", (x,y,z): (" << std::setw(7) << tc_x->at(itc) << ", " << std::setw(7) << tc_y->at(itc) << ", " << std::setw(7) << tc_z->at(itc) << ")"
		  <<", (pt,eta,phi,energy): (" << std::fixed << std::setprecision(2) << std::setw(8) << tc_pt->at(itc)*1.e3
		  << ", " << std::setw(8) << tc_eta->at(itc) 
		  << ", " << std::setw(8) << (TMath::RadToDeg()*tc_phi->at(itc)) 
		  << ", " << std::setw(8) << tc_energy->at(itc)*1.e3 << ") "
		  << std::defaultfloat
		  << std::endl;
      }
      
      for (uint32_t isect = 0 ; isect < 3 ; isect++ ){
	uint32_t addisect = 0;
	if(tc_z->at(itc)>0.0) addisect = 3;
	tcf0.setZero();
	tcf0.setROverZPhiF(tc_x->at(itc)/z,tc_y->at(itc)/z,isect+addisect);
	if(tcf0.getXOverZF()>=0.0){
	  float scale = ptcorr;
	  //float scale = 1.0;
	  //if(tc_layer->at(itc)==27 or tc_layer->at(itc)==29) scale = 1.5;
	  tcf0.setEnergyGeV(scale * tc_pt->at(itc));
	  tcf0.setLayer(tc_layer->at(itc));
	  //if(doPrint) tcf0.print();
	  //if((tc_eta->at(itc)<0. and genjetpt_negEta>-1.) or (tc_eta->at(itc)>0. and genjetpt_posEta>-1.))
	  vTcw[isect+addisect].push_back(tcf0);
	  tci_pt.push_back(scale * tc_pt->at(itc));
	  tci_layer.push_back(tc_layer->at(itc));
	  tci_x.push_back(tc_x->at(itc));
	  tci_y.push_back(tc_y->at(itc));
	  tci_z.push_back(z);
	  tci_sect.push_back(isect+addisect);
	}
      }//sector loop
      
      if(tc_eta->at(itc)<0. and genjetpt_negEta>-1.){
	//if(minRoz<tcRoZCut)
	totTCpt_negEta += tcpt;
	tcXoZnegEta += tcpt*tc_x->at(itc)/tc_z->at(itc) ;
	tcYoZnegEta += tcpt*tc_y->at(itc)/tc_z->at(itc) ;
	tcPhi_negEta += tcpt*tc_phi->at(itc) ;
	tcEta_negEta += tcpt*tc_eta->at(itc) ;
	nofNegEtaTCs++;
	// tcDxyozsize->Fill( (tc_x->at(itc)/tc_z->at(itc)-xoz_negEta), (tc_y->at(itc)/tc_z->at(itc)-yoz_negEta) );
	double droz = sqrt( (tc_x->at(itc)/tc_z->at(itc)-xoz_negEta)*(tc_x->at(itc)/tc_z->at(itc)-xoz_negEta) + (tc_y->at(itc)/tc_z->at(itc)-yoz_negEta)*(tc_y->at(itc)/tc_z->at(itc)-yoz_negEta) );
	if(droz>maxroz_negEta) maxroz_negEta = droz;
	// tcRoz->Fill(droz,tcpt);
      }else if(tc_eta->at(itc)>0. and genjetpt_posEta>-1.){
	//if(minRoz<tcRoZCut)
	totTCpt_posEta += tcpt;
	tcXoZposEta += tcpt*tc_x->at(itc)/tc_z->at(itc) ;
	tcYoZposEta += tcpt*tc_y->at(itc)/tc_z->at(itc) ;
	tcPhi_posEta += tcpt*tc_phi->at(itc) ;
	tcEta_posEta += tcpt*tc_eta->at(itc) ;
	nofPosEtaTCs++;
	// tcDxyozsize->Fill( (tc_x->at(itc)/tc_z->at(itc)-xoz_posEta), (tc_y->at(itc)/tc_z->at(itc)-yoz_posEta) );
	double droz = sqrt( (tc_x->at(itc)/tc_z->at(itc)-xoz_posEta)*(tc_x->at(itc)/tc_z->at(itc)-xoz_posEta) + (tc_y->at(itc)/tc_z->at(itc)-yoz_posEta)*(tc_y->at(itc)/tc_z->at(itc)-yoz_posEta) );
	if(droz>maxroz_posEta) maxroz_posEta = droz;
	// tcRoz->Fill(droz,tcpt);
      }
      
    }//end of TC loop
    ///////////////////========= Trigger Cell ==========/////////////////////////
    
    hSelEvents->Fill(3);
    ///////////////////=========== Emulation ============= ///////////////////////
    //std::vector<TPGClusterData> vCld[6];
    std::vector<TPGCluster> vCld[6];    
    for (uint32_t isect = 0 ; isect < 6 ; isect++ ) {
      s2Clustering.run(vTcw[isect],vCld[isect]);
      if(doPrint) std::cout << isect << ", Size of Tcs: " << vTcw[isect].size() << ", Size of Clusters: " << vCld[isect].size() << std::endl;
    }
    ///////////////////=========== Emulation ============= ///////////////////////

    hSelEvents->Fill(4);
    ///////////////////=========== Fill Cluster only details ============= ///////////////////////
    int nofClus1GeV = 0;
    int nofClus3GeV = 0;
    int nofClus5GeV = 0;
    int nofClus10GeV = 0;

    //tpgClus.clear();
    clus_pass.clear();
    clus_pt.clear();
    clus_pt_corr1D.clear();
    clus_pt_corr2D.clear();
    clus_eta.clear();
    clus_phi.clear();
    clus_sect.clear();
    clus_local_eta.clear();
    clus_local_phi.clear();
    
    tcf_pt.clear();
    tcf_layer.clear();
    tcf_xbyz.clear();
    tcf_ybyz.clear();
    tcf_zcm.clear();
    tcf_sect.clear();

    std::vector<int> pt_uncorr, pt1d, pt2d;
    for (uint32_t isect = 0 ; isect < 6 ; isect++ ){
      //for(TPGCluster const& clf : vCld[isect]){
      if(isect==0 or isect==3){
	nofClus1GeV = 0;
	nofClus3GeV = 0;
	nofClus5GeV = 0;
	nofClus10GeV = 0;
      }
      for(TPGTCFloats const& tcf : vTcw[isect]){
	tcf_pt.push_back(tcf.getEnergyGeV());
	tcf_layer.push_back(float(tcf.getLayer().to_uint()));
	tcf_xbyz.push_back(tcf.getXOverZF());
	tcf_ybyz.push_back(tcf.getYOverZF());
	tcf_zcm.push_back(tcf.getZCm());
	tcf_sect.push_back(isect);
      }
      int iclus = 0;
      for(TPGCluster const& clf : vCld[isect]){
	// clusxyOverZ[isect]->Fill(clf.getLocalXOverZF(),clf.getLocalYOverZF());
	// clusxy[isect]->Fill(clf.getGlobalXOverZF(isect),clf.getGlobalYOverZF(isect));
	// hClusLocalPhi[isect]->Fill(clf.getLocalPhiRad());
	// hClusGlobalPhi[isect]->Fill(clf.getGlobalPhiRad(isect));
	// hClusLocalPhiBits[isect]->Fill(clf.getLocalPhi());
	// hClusLocalEtaBits[isect]->Fill(clf.getEta());
	// hClusGlobalEta[isect]->Fill(clf.getGlobalEtaRad(isect));
	//tpgClus.push_back(clf);
	double tcPtSum = (clf.getGlobalEtaRad(isect) < 0) ? totTCpt_negEta : totTCpt_posEta ;
	int binPtClus = -1, binEtaClus = -1;
	GetHistoBin(hJetEtaPtBin,clf.getGlobalEtaRad(isect), clf.getEnergyGeV(), binEtaClus, binPtClus);
	double ClusE= -1., ClusE1D = -1., ClusE_UnCorr= -1.;
	double ptcorr1D = fClusCorr->Eval(clf.getEnergyGeV()) ;
	double ptcorr2D = fClusEtaCorr[binEtaClus-1]->Eval(clf.getEnergyGeV()) ;
	ptcorr1D = (ptcorr1D>8)?8.:ptcorr1D; //to avoid very low energy over-correction
	ptcorr2D = (ptcorr2D>8)?8.:ptcorr2D; //to avoid very low energy over-correction
	ClusE_UnCorr = clf.getEnergyGeV();
	ClusE1D = clf.getEnergyGeV() * ptcorr1D;
	ClusE = clf.getEnergyGeV() * ptcorr2D;
	
	bool isPresent = false;
	for(uint32_t ipt=0;ipt<pt_uncorr.size();ipt++) if(pt_uncorr.at(ipt)==TMath::Nint(ClusE_UnCorr)) isPresent = true ;	
	if(!isPresent) pt_uncorr.push_back(TMath::Nint(ClusE_UnCorr));
	
	isPresent = false;
	for(uint32_t ipt=0;ipt<pt1d.size();ipt++) if(pt1d.at(ipt)==TMath::Nint(ClusE1D)) isPresent = true ;
	if(!isPresent) pt1d.push_back(TMath::Nint(ClusE1D));
	
	isPresent = false;
	for(uint32_t ipt=0;ipt<pt2d.size();ipt++) if(pt2d.at(ipt)==TMath::Nint(ClusE)) isPresent = true ;	
	if(!isPresent) pt2d.push_back(TMath::Nint(ClusE));
	//doPrint = true;
	if(doPrint){
	  std::cout << "iclus-ievent: " << std::fixed << std::setprecision(5) << std::setw(7) << ievent
		    << ", sector: " << std::setw(5) << isect
		    << ", iclus: " << std::setw(5) << iclus++
		    <<", (x,y,z): (" << std::setw(7) << (clf.getGlobalXOverZF(isect) * clf.getZCm()) << ", " << std::setw(7) << (clf.getGlobalYOverZF(isect) * clf.getZCm()) << ", " << std::setw(7) << clf.getGlobalZCm(isect) << ")"
		    <<", (xoz,yoz,roz): (" << std::setw(5) << clf.getGlobalXOverZF(isect) << ", " << std::setw(5) << clf.getGlobalYOverZF(isect) << ", " << std::setw(5) << clf.getGlobalRhoOverZF(isect) << ")"
		    <<", (xoz,yoz): (" << std::setw(5) << clf.getLocalXOverZF() << ", " << std::setw(5) << clf.getLocalYOverZF() << ")"
		    <<", (pt,eta,phi): (" << std::fixed << std::setprecision(5) << std::setw(8) << clf.getEnergyGeV()
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
	doPrint = false;
	clus_pass.push_back(clf.getMaxFinderPass());
	clus_pt.push_back(clf.getEnergyGeV());
	clus_pt_corr1D.push_back(ClusE1D);
	clus_pt_corr2D.push_back(ClusE);
	clus_eta.push_back(clf.getGlobalEtaRad(isect));
	clus_phi.push_back(clf.getGlobalPhiRad(isect));
	clus_sect.push_back(isect);
	clus_local_eta.push_back(clf.getEtaRad());
	clus_local_phi.push_back(clf.getLocalPhiRad());
	
	if(clf.getEnergyGeV()>1.0) nofClus1GeV++;
	if(clf.getEnergyGeV()>3.0) nofClus3GeV++;	
	if(clf.getEnergyGeV()>5.0) nofClus5GeV++;
	if(clf.getEnergyGeV()>10.0) nofClus10GeV++;

      }//cluster loop
      // if(isect==2 or isect==5){
      // 	hNClus1GeV->Fill(nofClus1GeV);
      // 	hNClus3GeV->Fill(nofClus3GeV);
      // 	hNClus5GeV->Fill(nofClus5GeV);
      // 	hNClus10GeV->Fill(nofClus10GeV);
      // }
    }//sector loop
    outputTree->Fill();
    
    // float maxpt = 0;
    // for(uint32_t ipt=0;ipt<pt_uncorr.size();ipt++)  {
    //   hBkgRateSingle_Uncorr->Fill(pt_uncorr.at(ipt));
    //   if(pt_uncorr.at(ipt)>maxpt) maxpt = pt_uncorr.at(ipt);
    // }
    // hBkgRateSingle_Uncorr_MaxPt->Fill(maxpt);

    // maxpt = 0;
    // for(uint32_t ipt=0;ipt<pt1d.size();ipt++) {
    //   hBkgRateSingle1D->Fill(pt1d.at(ipt));
    //   if(pt1d.at(ipt)>maxpt) maxpt = pt1d.at(ipt);
    // }
    // hBkgRateSingle1D_MaxPt->Fill(maxpt);

    // maxpt = 0;
    // for(uint32_t ipt=0;ipt<pt2d.size();ipt++) {
    //   hBkgRateSingle2D->Fill(pt2d.at(ipt));
    //   if(pt2d.at(ipt)>maxpt) maxpt = pt2d.at(ipt);
    // }
    // hBkgRateSingle2D_MaxPt->Fill(maxpt);
    
    // // for(uint32_t ipt=0;ipt<pt2d.size();ipt++) {
    // //   std::cout << "ievent : " << ievent  << ", size: " << pt2d.size() << ", ipt: " << ipt << ", val : " <<  pt2d.at(ipt) << std::endl;
    // //   for(int ir=pt2d.at(ipt);ir>=0;ir--){
    // // 	//hBkgRateSingle2D->Fill(ir);
    // // 	std::cout << "\t ir: " << ir  << std::endl;
    // //   }
    // // }
    
    // pt_uncorr.clear();
    // pt1d.clear();
    // pt2d.clear();
    ///////////////////=========== Fill Cluster only details ============= ///////////////////////
    
    hSelEvents->Fill(5);
    bool hasAnyClus = false;
    bool hasMatchClus = false;
    // ///////////////////=========== comapre genjet and Cluster ============= ///////////////////////
    // for(int ipjet=0; ipjet < jetlist.size() ; ipjet++ ){
    //   int ijet = jetlist.at(ipjet).index ;
    //   //hEt->Fill(genjet_pt->at(ijet));
    //   int minsect = (genjet_eta->at(ijet) < 0) ? 0 : 3;
    //   int maxsect = (genjet_eta->at(ijet) < 0) ? 3 : 6;
    //   double tcPtSum = (genjet_eta->at(ijet) < 0) ? totTCpt_negEta : totTCpt_posEta ;
    //   double avgXoZ = (genjet_eta->at(ijet) < 0) ? tcXoZnegEta/tcPtSum : tcXoZposEta/tcPtSum ;
    //   double avgYoZ = (genjet_eta->at(ijet) < 0) ? tcYoZnegEta/tcPtSum : tcYoZposEta/tcPtSum ;
    //   double avgPhi = (genjet_eta->at(ijet) < 0) ? tcPhi_negEta/tcPtSum : tcPhi_posEta/tcPtSum ;
    //   double avgEta = (genjet_eta->at(ijet) < 0) ? tcEta_negEta/tcPtSum : tcEta_posEta/tcPtSum ;
    //   double gjroz = 1/sinh(genjet_eta->at(ijet)) ; //tan(2*atan(exp(-1.0*genjet_eta->at(ijet))));
    //   double gjxoz = gjroz*cos(genjet_phi->at(ijet));
    //   double gjyoz = gjroz*sin(genjet_phi->at(ijet));
    //   double avgRoZ = sqrt((avgXoZ-gjxoz)*(avgXoZ-gjxoz) + (avgYoZ-gjyoz)*(avgYoZ-gjyoz));
    //   bool hasFound = false, hasFoundPt = false;
    //   int binPt, binEta;
    //   int binPtClus, binEtaClus;
    //   double cluspt = genjet_pt->at(ijet);
    //   int nofClus = 0;
    //   //hTrigGen->Fill(genjet_pt->at(ijet));
    //   double minRoz = 1., closestClusE= 0., closestClusE1D = 0., closestClusE_UnCorr= 0., closestClusEDiff= 100., totClusE = 0., totClusE_UnCorr = 0.;
    //   double clstdeltaR = 0., clstEta = 0., clstPhi = 0.;
    //   double clstXoZ = 0., clstYoZ = 0., clstRoZ = 0., clstZcm = 0.;
    //   int clstisect = -1;
    //   for (uint32_t isect = minsect ; isect < maxsect ; isect++ ){
    // 	for(TPGCluster const& clf : vCld[isect]){	  
    // 	  if(clf.getEnergyGeV()<3.0) continue;
    // 	  hasAnyClus = true;
    // 	  //hClusE->Fill(clf.getEnergyGeV());
    // 	  double deltaR = TMath::Sqrt((clf.getGlobalEtaRad(isect)-genjet_eta->at(ijet))*(clf.getGlobalEtaRad(isect)-genjet_eta->at(ijet))
    // 				      +
    // 				      (clf.getGlobalPhiRad(isect)-genjet_phi->at(ijet))*(clf.getGlobalPhiRad(isect)-genjet_phi->at(ijet)));

    // 	  double dClusGenRoz = sqrt( (clf.getGlobalXOverZF(isect) - gjxoz)*(clf.getGlobalXOverZF(isect) - gjxoz) + (clf.getGlobalYOverZF(isect) - gjyoz)*(clf.getGlobalYOverZF(isect) - gjyoz) );
    // 	  GetHistoBin(hJetEtaPtBin,clf.getGlobalEtaRad(isect), clf.getEnergyGeV(), binEtaClus, binPtClus);
    // 	  double ptcorr1D = fClusCorr->Eval(clf.getEnergyGeV()) ;
    // 	  double ptcorr2D = fClusEtaCorr[binEtaClus-1]->Eval(clf.getEnergyGeV()) ;
    // 	  ptcorr1D = (ptcorr1D>8)?8.:ptcorr1D; //to avoid very low energy over-correction
    // 	  ptcorr2D = (ptcorr2D>8)?8.:ptcorr2D; //to avoid very low energy over-correction
    // 	  if(dClusGenRoz<minRoz){
    // 	    minRoz = dClusGenRoz;
    // 	    closestClusE_UnCorr = clf.getEnergyGeV();	    
    // 	    closestClusE1D = clf.getEnergyGeV() * ptcorr1D;
    // 	    closestClusE = clf.getEnergyGeV() * ptcorr2D;
    // 	    closestClusEDiff = closestClusE - genjet_pt->at(ijet);
    // 	    clstdeltaR = deltaR;
    // 	    clstisect = isect;
    // 	    clstEta = clf.getGlobalEtaRad(isect);
    // 	    clstPhi = clf.getGlobalPhiRad(isect);
    // 	    clstXoZ = clf.getGlobalXOverZF(isect);
    // 	    clstYoZ = clf.getGlobalYOverZF(isect);
    // 	    clstRoZ = clf.getGlobalRhoOverZF(isect);
    // 	    clstZcm = clf.getGlobalZCm(isect);
    // 	  }//minRoZ condition
    // 	  totClusE_UnCorr += clf.getEnergyGeV() ;
    // 	  totClusE += ( clf.getEnergyGeV()*ptcorr2D ) ;
    // 	  nofClus++;
    // 	}//cluster loop
    //   }//sector loop      

      
    //   double minRoz_cl3d = 1., closestClusE_cl3d= 0., closestClusE_cl3d_UnCorr = 0., closestClusEta_cl3d= 0.;
    //   for(int icl3=0;icl3<cl3d_n;icl3++){	  
    // 	if(cl3d_pt->at(icl3)<3.0) continue;
    // 	double clusXoZ = cos(cl3d_phi->at(icl3))/sinh(cl3d_eta->at(icl3)) ; 
    // 	double clusYoZ = sin(cl3d_phi->at(icl3))/sinh(cl3d_eta->at(icl3)); 
    // 	double dClusGenRoz = sqrt( (clusXoZ - gjxoz)*(clusXoZ - gjxoz) + (clusYoZ - gjyoz)*(clusYoZ - gjyoz) );
    // 	if(dClusGenRoz<minRoz_cl3d){
    // 	  minRoz_cl3d = dClusGenRoz;
    // 	  closestClusE_cl3d_UnCorr = cl3d_pt->at(icl3) ;
    // 	  closestClusE_cl3d = cl3d_pt->at(icl3) * fcl3d->Eval(cl3d_pt->at(icl3)) ;
    // 	  closestClusEta_cl3d = cl3d_eta->at(icl3);
    // 	}//minRoZ condition
    //   }//cluster loop
      
    //   bool hasClosestFound = (minRoz<1.)?true:false;
    //   bool hasClosestFound_cl3d = (minRoz_cl3d<1.)?true:false;
    //   bool hasClosestFoundPt = (hasClosestFound and closestClusE>pt_clusThresh)?true:false;
    //   bool hasClosestFoundPt1D = (hasClosestFound and closestClusE1D>pt_clusThresh)?true:false;
    //   bool hasClosestFoundPt_cl3d = (hasClosestFound_cl3d and closestClusE_cl3d>pt_clusThresh)?true:false;
    //   bool hasClosestFoundPt_UnCorr = (hasClosestFound and closestClusE_UnCorr>pt_clusThresh)?true:false;
    //   bool hasClosestFoundPt_cl3d_UnCorr = (hasClosestFound_cl3d and closestClusE_cl3d_UnCorr>pt_clusThresh)?true:false;
    //   bool hasClosestFoundPtTot = (hasClosestFound and totClusE>(pt_clusThresh))?true:false;
    //   bool hasClosestFoundEta = (hasClosestFound and closestClusE>pt_clusThresh_etaeff)?true:false;
    //   bool hgcalInnerEta = (fabs(genjet_eta->at(ijet))>innerEtaMin and fabs(genjet_eta->at(ijet))<innerEtaMax) ? true : false;
    //   bool hgcalOuterEta = (fabs(genjet_eta->at(ijet))>etaMin and fabs(genjet_eta->at(ijet))<etaMax) ? true : false;

    //   hasMatchClus = (hgcalOuterEta and hasClosestFound);
    //   rateSingles1D->Fill((hgcalOuterEta and hasClosestFound), closestClusE1D);
    //   rateSingles2D->Fill((hgcalOuterEta and hasClosestFound), closestClusE);
    //   rateSingles_Uncorr->Fill((hgcalOuterEta and hasClosestFound), closestClusE_UnCorr);

    //   hTotClus->Fill(closestClusE_UnCorr);
    //   if(hasClosestFound) hClusPtNoCut->Fill(closestClusE_UnCorr);
      
    //   if(hgcalOuterEta and hasClosestFound) {

    // 	hClusPt1D->Fill(closestClusE1D);
    // 	hClusPt2D->Fill(closestClusE);
    // 	hClusPt_Uncorr->Fill(closestClusE_UnCorr);
	
    // 	if(hgcalInnerEta){
    // 	  effTrigGen->Fill(hasClosestFoundPt,genjet_pt->at(ijet));
    // 	  effTrigGenTDR->Fill(hasClosestFoundPt,genjet_pt->at(ijet));
    // 	  effTrigGenTDR1D->Fill(hasClosestFoundPt1D,genjet_pt->at(ijet));
    // 	  effTrigGenTDR_UnCorr->Fill(hasClosestFoundPt_UnCorr,genjet_pt->at(ijet));
	
    // 	  h2GenVsClusPt->Fill(genjet_pt->at(ijet), closestClusE);
    // 	  h2GenVsTCPt->Fill(genjet_pt->at(ijet), tcPtSum);
    // 	  h2TCVsClusPt->Fill( tcPtSum, closestClusE);
    // 	  if(hasClosestFound_cl3d){
    // 	    h2NewVsOldPtcl3d->Fill( closestClusE_cl3d , closestClusE);
    // 	  }
    // 	}	
    //   }
    //   if(genjet_pt->at(ijet)>pt_clusThresh_etaeff and hasClosestFound)
    // 	effTrigGenEta->Fill(hasClosestFoundPt,fabs(genjet_eta->at(ijet)));
      
    //   if(genjet_pt->at(ijet)>pt_clusThresh_etaeff and hgcalInnerEta)
    // 	effTrigGenPhi->Fill(hasClosestFoundPt,TMath::RadToDeg()*genjet_phi->at(ijet));

    //   if(hgcalInnerEta and hasClosestFound_cl3d) {
    // 	h2GenVsClusPtcl3d->Fill(genjet_pt->at(ijet),closestClusE_cl3d);
    // 	effTrigGenTDR_cl3d->Fill(hasClosestFoundPt_cl3d,genjet_pt->at(ijet));
    // 	effTrigGenTDR_cl3d_UnCorr->Fill(hasClosestFoundPt_cl3d_UnCorr,genjet_pt->at(ijet));
    //   }

    // }//jet from decay of pions from tau
    // ///////////////////=========== comapre genjet and Cluster ============= ///////////////////////
    if(hasAnyClus) hSelEvents->Fill(6);
    if(hasMatchClus) hSelEvents->Fill(7);
    
    for(uint16_t i=0;i<6;i++) {
      vCld[i].clear();
      vTcw[i].clear();
    }
    
    taudlist.clear();
    taugdlist.clear();
    partlist.clear();
    genlist.clear();
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
    tc_eta->clear();
    tc_phi->clear();
    tc_layer->clear();
    tc_x->clear();
    tc_y->clear();
    tc_z->clear();
    tc_energy->clear();
    tc_data->clear();
    tc_uncompressedCharge->clear();
    tc_subdet->clear();
    tc_wafertype->clear();
    // cl3d_pt->clear();
    // cl3d_phi->clear();
    // cl3d_eta->clear();

  }//event loop
  foutTree->Write();
  
  // std::string outname = "stage2SemiEmulator_" + outputfile_extn + "_" + index;
  
  // effTrigGen->SetStatisticOption(TEfficiency::kBJeffrey);
  // effTrigGenEta->SetStatisticOption(TEfficiency::kBJeffrey);
  // effTrigGenPhi->SetStatisticOption(TEfficiency::kBJeffrey);
  // //effTrigGenPhi->GetPaintedGraph()->SetMinimum(0.0);
 
  // hSelEvents->GetXaxis()->SetBinLabel(1,"NoCut");
  // hSelEvents->GetXaxis()->SetBinLabel(2,"NGJEvts");
  // hSelEvents->GetXaxis()->SetBinLabel(3,"NTrk=(1 or 2)");
  // hSelEvents->GetXaxis()->SetBinLabel(4,"Pre-emul");
  // hSelEvents->GetXaxis()->SetBinLabel(5,"Post-emul");
  // hSelEvents->GetXaxis()->SetBinLabel(6,"Post-ratecalc");
  // hSelEvents->GetXaxis()->SetBinLabel(7,"hasAnyClus");
  // hSelEvents->GetXaxis()->SetBinLabel(8,"hasMatchClus");
  
  
  // TFile *fout = new TFile(Form("%s.root",outname.c_str()),"recreate");
  // //////=== yardstick to match calibration ===////////
  // hSelEvents->Write();
  // h2GenVsClusPt->Write();
  // h2GenVsTCPt->Write();
  // h2TCVsClusPt->Write();
  // h2GenVsClusPtcl3d->Write();
  // h2NewVsOldPtcl3d->Write();
  // hNJets->Write();
  // hNJetsF1->Write();
  // //////////////////////////////////////////////  
  // hTotClus->Write();
  // hClusPtNoCut->Write();
  // hClusPt1D->Write();
  // hClusPt2D->Write();
  // hClusPt_Uncorr->Write();      
  // //////=== trigger efficiency histos ===////////
  // effTrigGen->Write();
  // effTrigGenTDR->Write();
  // effTrigGenTDR1D->Write();
  // effTrigGenTDR_UnCorr->Write();
  // effTrigGenTDR_cl3d->Write();
  // effTrigGenTDR_cl3d_UnCorr->Write();
  // //////////////////////////////////////////////  
  // effTrigGenEta->Write();
  // effTrigGenPhi->Write();
  // //////////////////////////////////////////////  
  // rateSingles1D->Write();
  // rateSingles2D->Write();
  // rateSingles_Uncorr->Write();
  // //////////////////////////////////////////////
  // hBkgRateSingle1D->Write();
  // hBkgRateSingle2D->Write();
  // hBkgRateSingle_Uncorr->Write();
  // hBkgRateSingle1D_MaxPt->Write();
  // hBkgRateSingle2D_MaxPt->Write();
  // hBkgRateSingle_Uncorr_MaxPt->Write();
  // fout->Close();
  // delete fout;
  
  // fin->Close();
  // delete fin;

  return true;
}
