/*
  g++ -Wall -Icommon/inc src/RelayTcds2Check.cpp -lyaml-cpp -o bin/RelayTcds2Check.exe
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
#include "RelayTcds2Check.h"

#define NumberOfLinks 2

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
  sout << "RelayTcds2Check_" << relayNumber;
  TFileHandler tfh(sout.str().c_str());

  
  // Create the file reader
  /*
    Hgcal10gLinkReceiver::RelayReader _relayReader;
    _relayReader.enableLink(0,true);
    _relayReader.enableLink(1,false);
    _relayReader.enableLink(2,false);
  */
  
  Hgcal10gLinkReceiver::FileReader _relayReader;
  Hgcal10gLinkReceiver::FileReader runReader[NumberOfLinks];

  _relayReader.openRelay(relayNumber);

  std::ostringstream srun;
  srun << "dat/Relay" << relayNumber;

  runReader[0].setDirectory(srun.str());
  runReader[1].setDirectory(srun.str());
  runReader[2].setDirectory(srun.str());
  
  Hgcal10gLinkReceiver::RelayTcds2Check rtc;
  
  // Make the buffer space for the records and some useful casts for configuration and event records
  Hgcal10gLinkReceiver::RecordT<32*1024-1> *rRun(new Hgcal10gLinkReceiver::RecordT<32*1024-1>[NumberOfLinks]);

  //std::vector< Hgcal10gLinkReceiver::RecordT<4095>* > vRunning;
  std::vector<const Hgcal10gLinkReceiver::Record*> vRunning;
  vRunning.resize(NumberOfLinks);
  //vRunning.push_back(new Hgcal10gLinkReceiver::RecordT<4095>);
  //vRunning.push_back(new Hgcal10gLinkReceiver::RecordT<4095>);
  //vRunning.push_back(new Hgcal10gLinkReceiver::RecordT<4095>);
  
  Hgcal10gLinkReceiver::RecordT<32*1024-1> *rxxx(new Hgcal10gLinkReceiver::RecordT<32*1024-1>);
  assert(rxxx!=nullptr);
  
  unsigned nEvents(0);
  uint64_t relayUtc(0);
  uint64_t runUtc;

  uint32_t nRunningRecords_[NumberOfLinks];
  std::memset(nRunningRecords_,0,4*NumberOfLinks);

  bool repeated[NumberOfLinks];
  std::memset(repeated,0,sizeof(bool)*NumberOfLinks);

  bool stopped[NumberOfLinks];
  std::memset(stopped,0,sizeof(bool)*NumberOfLinks);
  
  Hgcal10gLinkReceiver::RecordHeader rh[NumberOfLinks],rhFirst[NumberOfLinks];
  Hgcal10gLinkReceiver::SlinkBoe boe_[NumberOfLinks];
  Hgcal10gLinkReceiver::SlinkEoe eoe_[NumberOfLinks];

  bool firstEventInRun(true);
  
  while(_relayReader.read(rxxx)) {
    std::cout << "Relay read successful, record state = " << Hgcal10gLinkReceiver::FsmState::stateName(rxxx->state()) << std::endl;
    rxxx->RecordHeader::print();
    
    if(rxxx->state()==Hgcal10gLinkReceiver::FsmState::Starting) {
      Hgcal10gLinkReceiver::RecordYaml *rCfg((Hgcal10gLinkReceiver::RecordYaml*)rxxx);
      YAML::Node n(YAML::Load(rCfg->string()));
      std::cout << n << std::endl << std::endl;
      
      unsigned runNumber(n["RunNumber"].as<unsigned>());
      for(unsigned l(0);l<NumberOfLinks;l++) {
	runReader[l].openRun(runNumber,l);
	runReader[l].read(rRun+l);
      }
      
      firstEventInRun=true;
    }
    
    if(rxxx->state()==Hgcal10gLinkReceiver::FsmState::Stopping) {
      bool done(false);
      
      while(!done) {
	done=true;
	uint32_t eventId(0xffffffff);

	if(firstEventInRun) {
	  for(unsigned l(0);l<NumberOfLinks;l++) {
	    assert(runReader[l].read(rRun+l));
	    rhFirst[l]=*((Hgcal10gLinkReceiver::RecordHeader*)(rRun+l));
	    std::cout << "INFO: first running record for link " << l << std::endl;
	    rhFirst[l].print();
	  }
	  firstEventInRun=false;
	}

	for(unsigned l(0);l<NumberOfLinks;l++) {
	  if(!stopped[l] && !repeated[l]) {
	    if(runReader[l].read(rRun+l)) {
	      if(rRun[l].state()!=Hgcal10gLinkReceiver::FsmState::Running) stopped[l]=true;

	      if(!stopped[l]) {
		done=false;

		vRunning[l]=(Hgcal10gLinkReceiver::Record*)(rRun+l);
		
		Hgcal10gLinkReceiver::RecordRunning *rEvt((Hgcal10gLinkReceiver::RecordRunning*)(rRun+l));
		const Hgcal10gLinkReceiver::SlinkBoe *b(rEvt->slinkBoe());
		assert(b!=nullptr);
		if(eventId>b->eventId()) eventId=b->eventId();

	      } else {
		std::cout << "WARNING: read stopped for link " << l << std::endl;
		vRunning[l]=nullptr;
	      }
	      
	    } else {
	      std::cout << "WARNING: read failed for link " << l << std::endl;
	      vRunning[l]=nullptr;
	    }
	  }
	}

	for(unsigned l(0);l<NumberOfLinks;l++) {
	  if(!stopped[l]) {
	    Hgcal10gLinkReceiver::RecordRunning *rEvt((Hgcal10gLinkReceiver::RecordRunning*)(rRun+l));
	    const Hgcal10gLinkReceiver::SlinkBoe *b(rEvt->slinkBoe());
	    assert(b!=nullptr);

	    const Hgcal10gLinkReceiver::SlinkEoe *e(rEvt->slinkEoe());
	    assert(e!=nullptr);

	    if(b->eventId()>eventId) {
	      if(!repeated[l]) {
		std::cout << "ERROR: running records for link " << l << " = " << std::setw(10) << nRunningRecords_[l] << std::endl;
		rh[l].RecordHeader::print();
		rEvt->RecordHeader::print();

		boe_[l].print();
		b->print();
		
		uint32_t diffEventId(b->eventId()-eventId);
		uint32_t diffSequence(rEvt->sequenceCounter()-rh[l].sequenceCounter());
		assert(diffEventId==diffSequence);
	      }
	      repeated[l]=true;

	    } else {
	      repeated[l]=false;

	      // Save for next event
	      rh[l]=*((Hgcal10gLinkReceiver::RecordHeader*)(rRun+l));
	      boe_[l]=*b;
	      eoe_[l]=*e;
	    }
	  }
	}
	
	//rtc.runningRecord(vRunning);
	/*
	unsigned eventId(0xffffffff);
      
	for(unsigned l(0);l<r.size();l++) {
	  if(r[l]!=nullptr) {
	    Hgcal10gLinkReceiver::RecordRunning *rEvt((Hgcal10gLinkReceiver::RecordRunning*)r[l]);
	    nRunningRecords_[l]++;

	    const Hgcal10gLinkReceiver::SlinkBoe *b(rEvt->slinkBoe());
	    assert(b!=nullptr);
	    if(nRunningRecords_[l]>200000000) {
	      std::cout << "INFO: running records for link " << l << " = " << std::setw(10) << nRunningRecords_[l] << std::endl;
	      b->print();
	    }
	  
	    const Hgcal10gLinkReceiver::SlinkEoe *e(rEvt->slinkEoe());
	    assert(e!=nullptr);
	    //e->print();
	  
	    if(eventId==0xffffffff) {
	      eventId=b->eventId();
	      if(nRunningRecords_[0]>200000000) {
		std::cout << "INFO: running records for link " << l << " = " << std::setw(10) << nRunningRecords_[l]
			  << ", event id set to " << eventId << std::endl;
	      }
	    } else {
	      if(nRunningRecords_[0]>200000000) {
		std::cout << "INFO: running records for link " << l << " = " << std::setw(10) << nRunningRecords_[l]
			  << ", event id check " << eventId << " vs " << b->eventId() << std::endl;
	      }
	      if(eventId!=b->eventId()) {
		std::cout << "WARNING: running records for link " << l << " = " << std::setw(10) << nRunningRecords_[l]
			  << ", event id check failed " << eventId << " vs " << b->eventId() << std::endl;
		rh[l].print();
		boe_[l].print();
		r[l]->RecordHeader::print();
		b->print();
	      }
	    }

	    // Save for next event
	    rh[l]=*((Hgcal10gLinkReceiver::RecordHeader*)r[l]);
	    boe_[l]=*b;
	    eoe_[l]=*e;

	  } else {
	    std::cout << "WARNING: nullptr for running records for link " << l << " = " << std::setw(10) << nRunningRecords_[l] << std::endl;
	  }
	}
	*/

      }
    }
    
    if(rxxx->state()==Hgcal10gLinkReceiver::FsmState::Status) {
      Hgcal10gLinkReceiver::RecordYaml *rCfg((Hgcal10gLinkReceiver::RecordYaml*)rxxx);
      YAML::Node n(YAML::Load(rCfg->string()));

      //if(n["Source"].as<std::string>()=="TCDS2") {
      std::cout <<  std::endl << n << std::endl << std::endl;
	// L1A total
	//}
    }
    if(rxxx->state()==Hgcal10gLinkReceiver::FsmState::Configuration) {
      Hgcal10gLinkReceiver::RecordYaml *rCfg((Hgcal10gLinkReceiver::RecordYaml*)rxxx);
      std::cout << rCfg->string() << std::endl;
      //rCfg->print();
      //YAML::Node n(YAML::Load(rCfg->string()));

      //if(n["Source"].as<std::string>()=="TCDS2") {
      //std::cout << std::endl << n << std::endl << std::endl;
	// L1A total
	//}
    }
  }
    
  //rtc.nonRunningRecord(rxxx);
  
  std::cout << "Total number of events found = " << nEvents << std::endl;
  
  delete rxxx;
  
  return 0;
}
