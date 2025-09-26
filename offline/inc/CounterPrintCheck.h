#ifndef Hgcal10gLinkReceiver_CounterPrintCheck_h
#define Hgcal10gLinkReceiver_CounterPrintCheck_h

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

#include "SlinkBoe.h"
#include "SlinkEoe.h"
#include "BePacketHeader.h"
#include "OrbitReader.h"

class CounterPrintCheck {
public:
  CounterPrintCheck() {
  }

  bool runStart(uint32_t nRun, uint32_t sid) {
    std::cout << "CounterPrintCheck: called runStart for run " << nRun
	      <<", source id " << sid << std::endl;
    
    _sourceId=sid;
    _nEvents=0;
    _initialOrbit=0;
    
     return true;
  }

  bool event(unsigned n64, const uint64_t *p) {
    //std::cout << "CounterPrintCheck: called event" << std::endl;
    
    Hgcal10gLinkReceiver::SlinkBoe *b(nullptr);
    if(n64>=2) {
      b=(Hgcal10gLinkReceiver::SlinkBoe*)p;
      //b->print();
      //std::cout << std::endl;
    }
    
    Hgcal10gLinkReceiver::SlinkEoe *e(nullptr);
    if(n64>=4) {
      e=(Hgcal10gLinkReceiver::SlinkEoe*)(p+n64-2);
      //e->print();
      //std::cout << std::endl;
    }

    if(_nEvents==0) {
      _initialOrbit=e->orbitId();
    }
    
    std::array<unsigned,20> array;
    array[0]=b->eventId();
    array[1]=e->orbitId();
    array[2]=e->orbitId()-_initialOrbit;
    array[3]=e->bxId();
    array[4]=b->eventId()%64;
    array[5]=e->orbitId()%8;
   
    Hgcal10gLinkReceiver::BePacketHeader *bph((Hgcal10gLinkReceiver::BePacketHeader*)(p+2));
    array[6]=bph->bunchCounter();
    array[7]=bph->eventCounter();
    array[8]=bph->orbitCounter();

    _nEvents++;
    std::cout << "Event " << std::setw(10) << _nEvents << ":";
    for(unsigned i(0);i<9;i++) {
      std::cout << " " << std::setw(6) << array[i];
      if(i==5) std::cout << "     | ";
    }
    if(array[3]!=array[6] || array[4]!=array[7] || array[5]!=array[8]) std::cout << "   <<===";
    std::cout << std::endl;
    
    
    return true;
  }

  bool runStop() {
    std::cout << "CounterPrintCheck: called runStop" << std::endl;
    return true;
  }
  
private:
  uint32_t _sourceId;
  uint32_t _nEvents;
  uint32_t _initialOrbit;
};

#endif
