#ifndef Hgcal10gLinkReceiver_RelayTcds2Check_h
#define Hgcal10gLinkReceiver_RelayTcds2Check_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <vector>
#include <string>
#include <cstring>

#include "RecordYaml.h"


namespace Hgcal10gLinkReceiver {

  class RelayTcds2Check {

  public:
  
    RelayTcds2Check() { 
      nRunningRecords_[0]=0;
      nRunningRecords_[1]=0;
      nRunningRecords_[2]=0;
    }
    
    virtual ~RelayTcds2Check() {
      for(unsigned l(0);l<3;l++) {
	std::cout << "Running records for link " << l << " = " << std::setw(10) << nRunningRecords_[l] << std::endl;
      }
    }
    
    void nonRunningRecord(const Record *r) {
      if(r->state()==Hgcal10gLinkReceiver::FsmState::Starting) {
        //if(relayUtc==0) relayUtc=rxxx->utc();
        //runUtc=rxxx->utc();
        //oldBx=0;
        //runCount++;
      }

      if(r->state()==Hgcal10gLinkReceiver::FsmState::Configuration) {
	Hgcal10gLinkReceiver::RecordYaml *rCfg((Hgcal10gLinkReceiver::RecordYaml*)r);
        YAML::Node n(YAML::Load(rCfg->string()));
        
        if(n["Source"].as<std::string>()=="Serenity") {
          std::cout << n << std::endl;
        }
        if(n["Source"].as<std::string>()=="TCDS2") {
          std::cout << n << std::endl;
        }
      }
      
      if(r->state()==Hgcal10gLinkReceiver::FsmState::Status) {
	Hgcal10gLinkReceiver::RecordYaml *rCfg((Hgcal10gLinkReceiver::RecordYaml*)r);
        YAML::Node n(YAML::Load(rCfg->string()));
        
        if(n["Source"].as<std::string>()=="Serenity") {
          std::cout << n << std::endl;
        }
        if(n["Source"].as<std::string>()=="TCDS2") {
          std::cout << n << std::endl;
        }
      }



    }
    
    void runningRecord(const std::vector<const Record*> &r) {

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
	      rh_[l].print();
	      boe_[l].print();
	      r[l]->RecordHeader::print();
	      b->print();
	    }
	  }

	  // Save for next event
	  rh_[l]=*((Hgcal10gLinkReceiver::RecordHeader*)r[l]);
	  boe_[l]=*b;
	  eoe_[l]=*e;

	} else {
	  std::cout << "WARNING: nullptr for running records for link " << l << " = " << std::setw(10) << nRunningRecords_[l] << std::endl;
	}
      }
    }

  protected:
    unsigned nRunningRecords_[3];
    RecordHeader rh_[3];
    Hgcal10gLinkReceiver::SlinkBoe boe_[3];
    Hgcal10gLinkReceiver::SlinkEoe eoe_[3];
  };

}

#endif
