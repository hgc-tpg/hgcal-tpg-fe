#ifndef Hgcal10gLinkReceiver_CounterCheck_h
#define Hgcal10gLinkReceiver_CounterCheck_h

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <thread>

#include "TH1D.h"
#include "TFileHandler.h"

#include "SlinkBoe.h"
#include "SlinkEoe.h"
#include "OrbitReader.h"

class CounterCheck {
public:
  CounterCheck() {
  }

  bool runStart(uint32_t nRun, uint32_t sid) {
    std::cout << "CounterCheck: called runStart for run " << nRun
	      <<", source id " << sid << std::endl;
    
    _sourceId=sid;

    _tfh.initialise(std::string("CounterCheckRun")+std::to_string(nRun)+"Sid"+std::to_string(sid));
    
    _hSlinkHeaderBx=new TH1D("SlinkHeaderBx",";Slink Header BX;Number of events",3600,0,3600);
    
    return true;
  }

  bool event(unsigned n64, const uint64_t *p) {
    std::cout << "CounterCheck: called event" << std::endl;
    
    if(n64>=2) {
      Hgcal10gLinkReceiver::SlinkBoe *b((Hgcal10gLinkReceiver::SlinkBoe*)p);
      b->print();
      std::cout << std::endl;
    }
    if(n64>=4) {
      Hgcal10gLinkReceiver::SlinkEoe *e((Hgcal10gLinkReceiver::SlinkEoe*)(p+n64-2));
      e->print();
      std::cout << std::endl;

      _hSlinkHeaderBx->Fill(e->bxId());
    }
    
    for(unsigned j(0);j<n64;j++) {
      std::cout << "Word " << std::setw(6) << j << " = 0x"
		<< std::hex << std::setfill('0')
		<< std::setw(16) << p[j]
		<< std::dec << std::setfill(' ')
	      << std::endl;
    }
    std::cout << std::endl;
    
    return true;
  }

  bool runStop() {
    std::cout << "CounterCheck: called runStop" << std::endl;
    return true;
  }
  
private:
  uint32_t _sourceId;

  TFileHandler _tfh;
  
  TH1D *_hSlinkHeaderBx;
};

#endif
