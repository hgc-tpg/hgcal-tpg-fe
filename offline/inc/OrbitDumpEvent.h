#ifndef Hgcal10gLinkReceiver_OrbitDumpEvent_h
#define Hgcal10gLinkReceiver_OrbitDumpEvent_h

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
#include "OrbitReader.h"

class OrbitDumpEvent {
public:
  OrbitDumpEvent() {
  }

  bool runStart(uint32_t run, uint32_t sid) {
    return true;
  }

  bool orbit(const Hgcal10gLinkReceiver::OrbitHeader &oh) {
    oh.print();
    std::cout << std::endl;

    return true;
  }

  bool event(const Hgcal10gLinkReceiver::OrbitReaderEvent &ore) {
    ore._ft->print();
    std::cout << std::endl;
    
    unsigned n64(2*ore._ft->fragmentSize());
    uint64_t *p(ore._array);

    if(n64>=2) {
      Hgcal10gLinkReceiver::SlinkBoe *b((Hgcal10gLinkReceiver::SlinkBoe*)p);
      b->print();
      std::cout << std::endl;
    }
    if(n64>=4) {
      Hgcal10gLinkReceiver::SlinkEoe *e((Hgcal10gLinkReceiver::SlinkEoe*)(p+n64-2));
      e->print();
      std::cout << std::endl;
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
    return true;
  }
  
};

#endif
