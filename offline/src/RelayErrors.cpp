/*
g++ -Wall -I. -Icommon/inc -Iexpert/inc src/RelayErrors.cpp -lyaml-cpp -o bin/RelayErrors.exe `root-config --libs --cflags`
*/

#include <iostream>
#include <iomanip>
#include <cassert>

#include <yaml-cpp/yaml.h>

#include "TH1D.h"
#include "TH2D.h"
#include "TFileHandler.h"

#include "RelayReader.h"
#include "EcondHeader.h"

bool printRelay1691487718(unsigned nEvents, unsigned statusBe, const Hgcal10gLinkReceiver::RecordRunning *rEvent) {
  return rEvent->payloadLength()==10 ||
    (rEvent->payloadLength()>6 && rEvent->payloadLength()<123 && statusBe==0) ||
    (rEvent->payloadLength()>6 && statusBe==4) ||
    nEvents==     152 ||
    nEvents==     154 ||
    nEvents==     902 ||
    nEvents==     904 ||
    nEvents==     905 ||
    nEvents==     920 ||
    nEvents==     923 ||
    nEvents==     924 ||
    nEvents==18437957 ||
    nEvents==18437958 ||
    nEvents==18437959 ||
    nEvents==18763176 ||
    nEvents==18763177 ||
    nEvents==18763178 ||
    nEvents==20998024 ||
    nEvents==20998025;
}


int main(int argc, char** argv) {
  if(argc<2) {
    std::cerr << argv[0] << ": No relay number specified" << std::endl;
    return 1;
  }

  uint32_t relayNumber(0);
  std::istringstream issRelay(argv[1]);
  issRelay >> relayNumber;

  if(relayNumber==0) {
    std::cerr << argv[0] << ": Relay number uninterpretable" << std::endl;
    return 2;
  }

  // Make the buffer space for the records and some useful casts for configuration and event records
  Hgcal10gLinkReceiver::RecordT<4*4095> *rxxx(new Hgcal10gLinkReceiver::RecordT<4*4095>);
  Hgcal10gLinkReceiver::RecordYaml *rCfg((Hgcal10gLinkReceiver::RecordYaml*)rxxx);
  Hgcal10gLinkReceiver::RecordRunning *rEvent((Hgcal10gLinkReceiver::RecordRunning*)rxxx);

  for(unsigned nLink(1);nLink<2;nLink++) {
  //for(unsigned nLink(2);nLink<3;nLink++) {
  
  std::ostringstream sout;
  sout << "RelayErrors_" << relayNumber << "_Link" << nLink;
  TFileHandler tfh(sout.str().c_str());

  TH1D *hEventTime,*hEventBx,*hEventBx6,*hEventBx10,*hEventBx106,*hSpillBx,*hPayloadSize,*hL1aDiff,*hL1aDiff6,*hL1aDiff10,*hL1aDiffS0,*hL1aDiffS4,*hL1aDiffS5,*hL1aDiffSX,
    *hStatusBe,*hRunCount,*hRunCount6,*hRunCount10,*hRunCountS5;
  TH2D *hPayloadSizeVsTime,*hPayloadSizeVsBx,*hPayloadSizeVsL1aDiff,*hPayloadSizeVsStatusBe;

  hSpillBx=new TH1D("SpillBx",";Time in run;Event rate (kHz)",300,20,26);
  hEventTime=new TH1D("EventTime",";Time in run;Event rate (kHz)",2500,0,2500);
  hRunCount=new TH1D("RunCount",";Run count;Number of L1As",100,0,100);
  hRunCount6=new TH1D("RunCount6",";Run count;Number of L1As",100,0,100);
  hRunCount10=new TH1D("RunCount10",";Run count;Number of L1As",100,0,100);
  hRunCountS5=new TH1D("RunCountS5",";Run count;Number of L1As",100,0,100);

  hStatusBe=new TH1D("StatusBe",";ECON-D BE status;Number of L1As",8,0,8);

  hEventBx=new TH1D("EventBx",";Bunch crossing;Number of L1As",3600,0,3600);
  hEventBx6=new TH1D("EventBx6",";Bunch crossing;Number of L1As",3600,0,3600);
  hEventBx10=new TH1D("EventBx10",";Bunch crossing;Number of L1As",3600,0,3600);
  hEventBx106=new TH1D("EventBx106",";Bunch crossing;Number of L1As",3600,0,3600);

  hPayloadSize=new TH1D("PayloadSize",";Payload size (64-bit words)",250,0,250);  
  hPayloadSizeVsTime=new TH2D("PayloadSizeVsTime",";Time in run;Payload size (64-bit words)",2500,0,2500,250,0,250);
  hPayloadSizeVsBx=new TH2D("PayloadSizeVsBx",";Bunch crossing;Payload size (64-bit words)",3600,0,3600,250,0,250);
  hPayloadSizeVsL1aDiff=new TH2D("PayloadSizeVsL1ADiff",";Difference of L1As (BX);Payload size (64-bit words)",1000,0,1000,250,0,250);
  hPayloadSizeVsStatusBe=new TH2D("PayloadSizeVsStatusBe",";ECON-D BE status;Payload size (64-bit words)",8,0,8,250,0,250);

  hL1aDiff=new TH1D("L1aDiff",";Difference of L1As (BX);Number of L1As",1000,0,1000);
  hL1aDiff6=new TH1D("L1aDiff6",";Difference of L1As (BX);Number of L1As",1000,0,1000);
  hL1aDiff10=new TH1D("L1aDiff10",";Difference of L1As (BX);Number of L1As",1000,0,1000);
  hL1aDiffSX=new TH1D("L1aDiffSX",";Difference of L1As (BX);Number of L1As",1000,0,1000);
  hL1aDiffS0=new TH1D("L1aDiffS0",";Difference of L1As (BX);Number of L1As",1000,0,1000);
  hL1aDiffS4=new TH1D("L1aDiffS4",";Difference of L1As (BX);Number of L1As",1000,0,1000);
  hL1aDiffS5=new TH1D("L1aDiffS5",";Difference of L1As (BX);Number of L1As",1000,0,1000);

  // Create the file reader
  Hgcal10gLinkReceiver::RelayReader _relayReader;
  _relayReader.enableLink(0,false);
  _relayReader.enableLink(1,nLink==1);
  _relayReader.enableLink(2,nLink==2);
  _relayReader.openRelay(relayNumber);
  
  unsigned nEvents(0);
  uint32_t runCount(0xffffffff);
  uint64_t relayUtc(0);
  uint64_t runUtc;
  uint64_t oldBx(0);

  unsigned nEvS4(0);
  unsigned nEvS5(0);
  unsigned oldStatusBe(0);
  
  while(_relayReader.read(rxxx)) {
    if(rxxx->state()!=Hgcal10gLinkReceiver::FsmState::Running) {
      
      if(rxxx->state()==Hgcal10gLinkReceiver::FsmState::Starting) {
	if(relayUtc==0) relayUtc=rxxx->utc();
	runUtc=rxxx->utc();
	oldBx=0;
	runCount++;
      }

      if(rxxx->state()==Hgcal10gLinkReceiver::FsmState::Configuration) {
	YAML::Node n(YAML::Load(rCfg->string()));
	
	if(n["Source"].as<std::string>()=="Serenity") {
	}
	if(n["Source"].as<std::string>()=="TCDS2") {
	}
      }
      
      if(rxxx->state()==Hgcal10gLinkReceiver::FsmState::Status) {
	YAML::Node n(YAML::Load(rCfg->string()));
	
	if(n["Source"].as<std::string>()=="Serenity") {
	}
	if(n["Source"].as<std::string>()=="TCDS2") {
	  //std::cout << n << std::endl;
	}
      }
      
    } else {
      nEvents++;

      bool printIt(false);
      
      // Access the Slink header ("begin-of-event")
      // This should always be present; check pattern is correct
      const Hgcal10gLinkReceiver::SlinkBoe *b(rEvent->slinkBoe());
      assert(b!=nullptr);
      //if(!b->validPattern()) b->print();
      
      // Access the Slink trailer ("end-of-event")
      // This should always be present; check pattern is correct
      const Hgcal10gLinkReceiver::SlinkEoe *e(rEvent->slinkEoe());
      assert(e!=nullptr);
      //if(!e->validPattern()) e->print();
      
      const Hgcal10gLinkReceiver::BePacketHeader *bph(rEvent->bePacketHeader());
      assert(bph!=nullptr);
	
      uint64_t tBx(e->totalBx());
      double timeBx(tBx*25.0e-9+runUtc-relayUtc);
      unsigned statusBe(bph->econdStatus(0));

      if(statusBe==4) nEvS4++;
      if(statusBe==5) nEvS5++;
      
      hSpillBx->Fill(50.0*timeBx);
      hEventTime->Fill(timeBx);
      hStatusBe->Fill(statusBe);
      hRunCount->Fill(runCount);
      if(rEvent->payloadLength()==  6) hRunCount6->Fill(runCount);
      if(rEvent->payloadLength()== 10) hRunCount10->Fill(runCount);
      if(statusBe==5) hRunCountS5->Fill(runCount);

      hEventBx->Fill(e->bxId());
      if(rEvent->payloadLength()==  6) hEventBx6->Fill(e->bxId());
      if(rEvent->payloadLength()== 10) hEventBx10->Fill(e->bxId());
      if(rEvent->payloadLength()==106) hEventBx106->Fill(e->bxId());

      hL1aDiff->Fill(tBx-oldBx);
      if(rEvent->payloadLength()== 6) hL1aDiff6->Fill(tBx-oldBx);
      if(rEvent->payloadLength()==10) hL1aDiff10->Fill(tBx-oldBx);
      if(oldStatusBe==0               ) hL1aDiffSX->Fill(tBx-oldBx);
      if(oldStatusBe==0 && statusBe==0) hL1aDiffS0->Fill(tBx-oldBx);
      if(oldStatusBe==0 && statusBe==4) hL1aDiffS4->Fill(tBx-oldBx);
      if(oldStatusBe==0 && statusBe==5) hL1aDiffS5->Fill(tBx-oldBx);
      
      hPayloadSize->Fill(rEvent->payloadLength());
      hPayloadSizeVsTime->Fill(timeBx,rEvent->payloadLength());
      hPayloadSizeVsBx->Fill(e->bxId(),rEvent->payloadLength());
      hPayloadSizeVsL1aDiff->Fill(tBx-oldBx,rEvent->payloadLength());
      hPayloadSizeVsStatusBe->Fill(statusBe,rEvent->payloadLength());

      oldBx=tBx;
      oldStatusBe=statusBe;

      // Set pointer to the beginning of the packet
      const uint64_t *p64(((const uint64_t*)rEvent)+1);

      // Check for odd missing packet padding
      if(rEvent->payloadLength()==6) {
	if(p64[3]!=0) {
	  std::cout << "Event " << nEvents << " has non-zero padding" << std::endl;
	  printIt=true;
	}
      }

      // Check for idles and/or repetition
      unsigned nRep(0);

      for(unsigned i(0);i<rEvent->payloadLength();i++) {
	if(i>0) {
	  if(p64[i]==p64[i-1]) nRep++;
	}

	if((p64[i]>>40)==0xaaaaff || ((p64[i]>>8)&0xffffff)==0xaaaaff) {
	  std::cout << "Event " << nEvents << " has idle for word " << i << std::endl;
	  printIt=true;
	}
      }
      
      if(nRep>1) {
	std::cout << "Event " << nEvents << " has " << nRep << " repetitions" << std::endl;
	printIt=true;
      }
      
      //if(nEvents>30000000) {
      //std::cout << "Total BX = " << tBx << ", time = " << tBx*25.0e-9 << std::endl;
      //}
	    
      //if(nEvents<=10 || rEvent->payloadLength()==10) {
      if((relayNumber==1691487718 && printRelay1691487718(nEvents,statusBe,rEvent)) ||
	 (statusBe==4 && nEvS4<10) ||
	 (statusBe==5 && nEvS5<10) ||
	 printIt) {
	
	std::cout << "Event " << nEvents << std::endl;
	b->print();
	if(relayNumber==1691487718 && nEvents==20998026) e=(const Hgcal10gLinkReceiver::SlinkEoe*)(p64+4); // Repair!
	e->print();				
	bph->print();

	if(rEvent->payloadLength()>6) {
	  const Hgcal10gLinkReceiver::EcondHeader *edh((const Hgcal10gLinkReceiver::EcondHeader*)(p64+3));
	  edh->print();
	}
	
	for(unsigned i(0);i<rEvent->payloadLength();i++) {
	  std::cout << "Word " << std::setw(3) << i << " ";
	  std::cout << std::hex << std::setfill('0');
	  std::cout << "0x" << std::setw(16) << p64[i] << std::endl;
	  std::cout << std::dec << std::setfill(' ');
	}
	std::cout << std::endl;
      }
    }
  }

  std::cout << "Total number of events found = " << nEvents << std::endl;
  std::cout << "Number of events with status 4 found = " << nEvS4 << std::endl;
  std::cout << "Number of events with status 5 found = " << nEvS5 << std::endl;
  }
  
  delete rxxx;
  
  return 0;
}
