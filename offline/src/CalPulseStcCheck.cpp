/*
g++ -Wall -Icommon/inc offline/src/CalPulseStcCheck.cpp -lyaml-cpp -o bin/CalPulseStcCheck.exe
*/

#include <iostream>
#include <iomanip>
#include <cassert>

#include <yaml-cpp/yaml.h>

#include "TH1D.h"
#include "TH2D.h"
#include "TFileHandler.h"

#include "RelayReader.h"
#include "EcontEnergies.h"



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

  std::ostringstream sout;
  sout << "CalPulseStcCheck_" << relayNumber;
  TFileHandler tfh(sout.str().c_str());

  TH2D *hStcVsRun=new TH2D("StcVsRun",";Run within relay;STC number;Total energy",
			   100,0,100,24,0,24);
  TH2D *hStcVsRunM=new TH2D("StcVsRunM",";Run within relay;STC number;Total energy",
			   100,0,100,24,0,24);
  TH2D *hStcVsRunP=new TH2D("StcVsRunP",";Run within relay;STC number;Total energy",
			   100,0,100,24,0,24);

  
  // Create the file reader
  Hgcal10gLinkReceiver::RelayReader _relayReader;
  _relayReader.enableLink(0,true);
  _relayReader.enableLink(1,false);
  _relayReader.enableLink(2,false);
  _relayReader.openRelay(relayNumber);
  
  // Make the buffer space for the records and some useful casts for configuration and event records
  Hgcal10gLinkReceiver::RecordT<16*4095> *rxxx(new Hgcal10gLinkReceiver::RecordT<16*4095>);
  Hgcal10gLinkReceiver::RecordYaml *rCfg((Hgcal10gLinkReceiver::RecordYaml*)rxxx);
  Hgcal10gLinkReceiver::RecordRunning *rEvent((Hgcal10gLinkReceiver::RecordRunning*)rxxx);

  unsigned nEvents(0);
  unsigned nEventsPerRun(0);
  unsigned nRuns(0);
  uint64_t relayUtc(0);
  uint64_t runUtc;


  std::vector<uint16_t> vStc[3];
  
  while(_relayReader.read(rxxx)) {

    if(rxxx->state()!=Hgcal10gLinkReceiver::FsmState::Running) {
      rxxx->RecordHeader::print();
      
      if(rxxx->state()==Hgcal10gLinkReceiver::FsmState::Starting) {
	nRuns++;
	nEventsPerRun=0;
	if(relayUtc==0) relayUtc=rxxx->utc();
	runUtc=rxxx->utc();
      }

      if(rxxx->state()==Hgcal10gLinkReceiver::FsmState::Configuration) {
	/*
	YAML::Node n(YAML::Load(rCfg->string()));
	
	if(n["Source"].as<std::string>()=="Serenity") {
	}
	if(n["Source"].as<std::string>()=="TCDS2") {
	}
	*/
      }
      
      if(rxxx->state()==Hgcal10gLinkReceiver::FsmState::Status) {
	/*
	YAML::Node n(YAML::Load(rCfg->string()));
	
	if(n["Source"].as<std::string>()=="Serenity") {
	}
	if(n["Source"].as<std::string>()=="TCDS2") {
	  std::cout << n << std::endl;
	}
	*/
      }
      
      if(rxxx->state()==Hgcal10gLinkReceiver::FsmState::Stopping) {
      }
      
    } else {
      nEvents++;
      
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
      
      // Set pointer to the beginning of the packet
      const uint64_t *p64(((const uint64_t*)rEvent)+1);

      Hgcal10gLinkReceiver::unpackerEnergies(p64+101,vStc[0]);
      Hgcal10gLinkReceiver::unpackerEnergies(p64+107,vStc[1]);
      Hgcal10gLinkReceiver::unpackerEnergies(p64+113,vStc[2]);

      for(unsigned i(0);i<24;i++) {
	//hStcVsRunM->Fill(nRuns-0.5,i+0.5,pow(2.0,vStc[0][i]/8.0));
	//hStcVsRun ->Fill(nRuns-0.5,i+0.5,pow(2.0,vStc[1][i]/8.0));
	//hStcVsRunP->Fill(nRuns-0.5,i+0.5,pow(2.0,vStc[2][i]/8.0));
	hStcVsRunM->Fill(nRuns-0.5,i+0.5,vStc[0][i]);
	hStcVsRun ->Fill(nRuns-0.5,i+0.5,vStc[1][i]);
	hStcVsRunP->Fill(nRuns-0.5,i+0.5,vStc[2][i]);
      }
      
      if(nEventsPerRun<=10 && nRuns==90) {
	for(unsigned i(0);i<rEvent->payloadLength();i++) {
	  std::cout << "Word " << std::setw(3) << i << " ";
	  std::cout << std::hex << std::setfill('0');
	  std::cout << "0x" << std::setw(16) << p64[i] << std::endl;
	  std::cout << std::dec << std::setfill(' ');
	}
	
	std::cout << std::endl << "STC energies" << std::endl;
	for(unsigned j(0);j<3;j++) {
	  for(unsigned i(0);i<24;i++) {
	    std::cout << " " << std::setw(5) << vStc[j][i];
	    if((i%12)==11) std::cout << std::endl;
	  }
	  std::cout << std::endl;
	}
      }
    }
  }

  std::cout << "Total number of events found = " << nEvents << std::endl;
  
  delete rxxx;
  
  return 0;
}
