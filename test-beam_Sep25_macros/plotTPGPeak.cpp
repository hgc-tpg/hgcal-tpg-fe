/**********************************************************************
 Created on : 28/09/2025
 Purpose    : Read the TPG data blocks and decode the energies
 Author     : Indranil Das, Research Associate, Imperial
 Email      : indranil.das@cern.ch | indra.ehep@gmail.com
**********************************************************************/
#include <iostream>

#include "TH1D.h"
#include "TFile.h"

#include "TFileHandlerLocal.h"
#include "FileReader.h"

#include "TPGFEDataformat.hh"
#include "TPGFEModuleEmulation.hh"
#include "Stage1IO.hh"

#include "SlinkBoe.h"
#include "SlinkEoe.h"
#include "TpgSubpacketHeader.h"
#include "OrbitReader.h"

// #include "OrbitDumpEvent.h"
// typedef OrbitDumpEvent OrbitCheckTypedef;

using namespace std;

///Source: offline/inc/OrbitCheck.hxx

class TPGEventReaderTB25 {
public:
  TPGEventReaderTB25(): iEvent(0) {    
  }

  bool runStart(uint32_t run, uint32_t sid) {
    return true;
  }

  bool orbit(const Hgcal10gLinkReceiver::OrbitHeader &oh) {
    //oh.print();
    std::cout << std::endl;

    return true;
  }

  bool event(const Hgcal10gLinkReceiver::OrbitReaderEvent &ore, TList *list) {

    //ore._ft->print();
    // std::cout << std::endl;
    
    unsigned n64(2*ore._ft->fragmentSize());
    uint64_t *p(ore._array);

    Hgcal10gLinkReceiver::SlinkBoe *b(nullptr);
    if(n64>=2) {
      b = ((Hgcal10gLinkReceiver::SlinkBoe*)p);
      // b->print();
      // std::cout << std::endl;
    }
    Hgcal10gLinkReceiver::SlinkEoe *e(nullptr);
    if(n64>=4) {
      e = ((Hgcal10gLinkReceiver::SlinkEoe*)(p+n64-2));
      // e->print();
      // std::cout << std::endl;
    }
    
    
    const Hgcal10gLinkReceiver::TpgSubpacketHeader *tsh(reinterpret_cast<const Hgcal10gLinkReceiver::TpgSubpacketHeader*>(b+1));
    const Hgcal10gLinkReceiver::TpgSubpacketHeader *tshEnd(reinterpret_cast<const Hgcal10gLinkReceiver::TpgSubpacketHeader*>(e-1));

    if(!tsh->validPattern()) {
      tsh->print();
      return false;
    }
    tsh=tsh->nextSubpacketHeader();

    int noffecafe = 0;
    bool done(false);
    while(tsh<=tshEnd && !done) {
      if(!tsh->validPattern()) {
	done=true;
	
      } else {
	if((tsh->channelId()%2)==0) { // RX only
	  //std::cout << "Processing event : " <<  iEvent << std::endl;
	  //tsh->print();
	  
	  unsigned emp_chan(tsh->channelId()/2);
	  //if(emp_chan==100 or emp_chan==102 or emp_chan==104 or emp_chan==106 or emp_chan==108 or emp_chan==110 or emp_chan==112 or emp_chan==114 or emp_chan==116 or emp_chan==118 or emp_chan==122){
	  if(emp_chan==100 or emp_chan==102 or emp_chan==104 or emp_chan==106 or emp_chan==108 or emp_chan==110 or emp_chan==112 or emp_chan==114 or emp_chan==116 or emp_chan==118 or emp_chan==120 or emp_chan==122){
	  //if(emp_chan==100){
	    uint wpspd = 0;
	    for(unsigned bx(0);bx<tsh->numberOfBxs();bx++) {
	      const uint64_t *el64packed((const uint64_t*)(tsh+1+bx*tsh->numberOfWordsPerBx()));
	      uint32_t elinks[7];
	      for(unsigned j(0);j<tsh->numberOfWordsPerBx();j++) {
		if(j<(tsh->numberOfWordsPerBx()-1)){
		  elinks[2*j] = el64packed[j] & 0xffffffff;
		  elinks[2*j+1] = (el64packed[j]>>32) & 0xffffffff;
		}else{
		  elinks[2*j] = el64packed[j] & 0xffffffff;
		}
		// std::cout << "Word " << std::setw(6) << wpspd++ << " = 0x"
		// 	  << std::hex << std::setfill('0')
		// 	  << std::setw(16) << el64packed[j]
		// 	  << std::dec << std::setfill(' ')
		// 	  << std::endl;	      
	      }
	      // for(unsigned iel(0);iel<7;iel++) {
	      // 	std::cout << "\t elink " << std::setw(3) << iel << " = 0x"
	      // 		  << std::hex << std::setfill('0')
	      // 		  << std::setw(8) << elinks[iel]
	      // 		  << std::dec << std::setfill(' ')
	      // 		  << std::endl;	      
	      // }
	      
	      int neTx = (emp_chan==120)?3:4;
	      uint32_t *el = new uint32_t[neTx];
	      el[0] = elinks[2];
	      el[1] = elinks[1];
	      el[2] = elinks[0];
	      if(neTx>3) el[3] = elinks[3];
	      // for(unsigned iel(0);iel<neTx;iel++){
	      // 	std::cout << "\t\t el " << std::setw(3) << iel << " = 0x"
	      // 		  << std::hex << std::setfill('0')
	      // 		  << std::setw(8) << el[iel]
	      // 		  << std::dec << std::setfill(' ')
	      // 		  << std::endl;	      		
	      // }
	      
	      TPGFEDataformat::TcRawDataPacket rdp;
	      if(neTx==4)
		TPGStage1Emulation::Stage1IO::convertElinksToTcRawData(TPGFEDataformat::BestC, 9, el, rdp);
	      else
		TPGStage1Emulation::Stage1IO::convertElinksToTcRawData(TPGFEDataformat::BestC, 6, el, rdp);
	      delete []el;
	      //rdp.print();
	      //std::cout<<"EmpCh : " << emp_chan << ", bx : "<< rdp.bx() << ", ModuleSum: " << TPGFEDataformat::TcRawData::Decode5E3M(rdp.moduleSum()) << std::endl;
	      uint64_t modsum = TPGFEDataformat::TcRawData::Decode5E3M(rdp.moduleSum());
	      TH1D *hE = (TH1D *) list->FindObject(Form("hEMS_empch%d",emp_chan)) ;
	      hE->SetBinContent(bx+1, (hE->GetBinContent(bx+1)+modsum) );
	    }
	  }
	  //std::cout<<std::endl; //after each lpGBT links
	}
	
	tsh=tsh->nextSubpacketHeader();
      }
      noffecafe++;
    }

    // for(unsigned j(0);j<n64;j++) {
    //   std::cout << "Word " << std::setw(6) << j << " = 0x"
    // 		<< std::hex << std::setfill('0')
    // 		<< std::setw(16) << p[j]
    // 		<< std::dec << std::setfill(' ')
    // 	      << std::endl;
    // }
    // std::cout << std::endl;
 
    //std::cout << "Processing event : " <<  iEvent++ << ", noffecafe: " << noffecafe << std::endl;
    
    return true;
  }

  bool runStop() {
    return true;
  }
private:
  uint64_t iEvent;
};

typedef TPGEventReaderTB25 OrbitCheckTypedef;

int main(int argc, char** argv){
  
  std::cout << "Nof arguments : " << argc << std::endl;
  if(argc < 3){
    std::cerr << argv[0] << ": no run numbers specified" << std::endl;
    return 1;
  }
  
  //Command line arg assignment
  //Assign relay and run numbers
  uint32_t runNumber(12600113);
  uint32_t sourceId(1260);
  uint32_t firstLs(1);
  
  unsigned dumpEvent(0);
  
  std::istringstream issRun(argv[1]);
  issRun >> runNumber;
  std::istringstream issLink(argv[2]);
  issLink >> sourceId;
  if(argc > 3){
    std::istringstream isfirstLs(argv[3]);
    isfirstLs >> firstLs;
  }
  
  /////=======================================================
  int emp_ch[12] = {100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120, 122} ;
  TFile *fout = new TFile(Form("output_run%u.root",runNumber),"recreate");
  TDirectory *dir_hist = fout->mkdir(Form("run%u",runNumber));
  //Create hists
  for(int iemp=0;iemp<12;iemp++){
    TH1D *hEMS = new TH1D(Form("hEMS_empch%d",emp_ch[iemp]),Form("Run:%u MS distn of layer %d (EMPch:%d)",runNumber,(iemp+1),emp_ch[iemp]),7,-3,3);
    hEMS->GetXaxis()->SetTitle("Bx");
    hEMS->GetYaxis()->SetTitle("Entries");
    hEMS->GetXaxis()->SetTitleSize(0.05);
    hEMS->GetYaxis()->SetTitleSize(0.05);
    hEMS->GetXaxis()->SetLabelSize(0.05);
    hEMS->GetYaxis()->SetLabelSize(0.05);
    hEMS->SetLineWidth(3);
    hEMS->SetLineColor(kBlue);
    hEMS->SetFillColor(kGreen);
    hEMS->SetDirectory(dir_hist);
  }
  TList *list = (TList *)dir_hist->GetList();
  /////========================================================
  
  std::vector<Hgcal10gLinkReceiver::OrbitReaderEvent> vEvents;
  
  OrbitCheckTypedef ct;
  assert(ct.runStart(runNumber,sourceId));
  
  std::string oDir("dat/");
  
  Hgcal10gLinkReceiver::OrbitReader oReader;
  oReader.setPrint(false);
  
  unsigned nEvents(0);
  
  const Hgcal10gLinkReceiver::OrbitHeader *oh;
  //Hgcal10gLinkReceiver::FragmentTrailer *ft;
  
  bool done(false);
  for(unsigned i(firstLs);!done;i++) {
    std::ostringstream oss;
    oss << oDir << std::setfill('0') << "run" << std::setw(6) << runNumber
	<< "/run" << std::setw(6) << runNumber << "_ls" << std::setw(4)
	<< i << "_index000000_source" << std::setw(4) << sourceId << ".raw";
    
    if(!oReader.open(oss.str())) {
      if(i==firstLs) std::cout << "Failed " << oss.str() << std::endl << std::endl;
      else           std::cout << "No more files to open" << std::endl << std::endl;
      done=true;
      continue;
    }
    std::cout << "Opened " << oss.str() << std::endl << std::endl;
    
    //while((oh=oReader.readOrbit(vEvents))!=nullptr and nEvents<1) {
    while((oh=oReader.readOrbit(vEvents))!=nullptr) {
      assert(ct.orbit(*oh));
      
      for(unsigned j(0);j<vEvents.size();j++) {
	assert(ct.event(vEvents[j],list));
	nEvents++;
      }
      
    }
    
    oReader.close();
  }

  std::cout << "Total number of events in run = " << nEvents << std::endl;

  assert(ct.runStop());
  /////========================================================
  TCanvas *c1 = new TCanvas(Form("Run%u_source%d",runNumber,sourceId),Form("Run%u_source%d",runNumber,sourceId),1600,800);
  c1->Divide(4,3);
  for(int iemp=0;iemp<12;iemp++){
    c1->cd(iemp+1);//->SetLogy();
    TH1D *hE = (TH1D *) list->FindObject(Form("hEMS_empch%d",emp_ch[iemp])) ;
    hE->Draw();
  }

  /////========================================================
  fout->cd();
  dir_hist->Write();
  c1->Write();
  fout->Close();
  delete fout;
  /////========================================================
  
  return 0;
}
