// g++ -std=c++11 -I/opt/local/include -I ../hgcal10glinkreceiver OrbitReader.cpp -L /opt/local/lib -l yaml-cpp

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

#include "yaml-cpp/yaml.h"

#include "OrbitReader.h"
#include "SlinkBoe.h"
#include "SlinkEoe.h"


void printEvent(uint32_t n, const uint64_t *p) {
  //std::cout << "printEvent: n = " << n << ", p = "
  //	    << (p==nullptr?"nullptr":std::to_string((uint64_t)p)) << std::endl;

  if(n>=2) {
    Hgcal10gLinkReceiver::SlinkBoe *b((Hgcal10gLinkReceiver::SlinkBoe*)p);
    b->print();
    std::cout << std::endl;
  }
  if(n>=4) {
    Hgcal10gLinkReceiver::SlinkEoe *e((Hgcal10gLinkReceiver::SlinkEoe*)(p+n-2));
    e->print();
    std::cout << std::endl;
  }

  for(unsigned j(0);j<n;j++) {
    std::cout << "Word " << std::setw(6) << j << " = 0x"
	      << std::hex << std::setfill('0')
	      << std::setw(16) << p[j]
	      << std::dec << std::setfill(' ')
	      << std::endl;
  }
  std::cout << std::endl;
}

int main(int argc, char *argv[]) {

  std::ifstream fout;
  //fout.open("output-stream-1hz.raw",std::ios::binary);
  fout.open("dth_p2_test.raw_sourceid01230_index000.raw");
  
  uint64_t array[16*1024];

  Hgcal10gLinkReceiver::OrbitHeader *orbitHeader;
  Hgcal10gLinkReceiver::FragmentTrailer *fragmentTrailer;

  Hgcal10gLinkReceiver::SlinkBoe *slinkBoe;
  Hgcal10gLinkReceiver::SlinkEoe *slinkEoe;

  std::vector< std::pair<unsigned,uint64_t*> > vEvents;
  
  for(unsigned i(0);i<10;i++) {
    fout.read((char*)array,8*4);
    orbitHeader=(Hgcal10gLinkReceiver::OrbitHeader*)(array);

    orbitHeader->print();
    std::cout << std::endl;

    vEvents.resize(orbitHeader->eventCount());
    
    fout.read((char*)(array+4),16*(orbitHeader->packetWordCount()-2));
    
    unsigned nEnd(2*orbitHeader->packetWordCount()-2);
    
    for(unsigned j(0);j<orbitHeader->eventCount();j++) {
      fragmentTrailer=(Hgcal10gLinkReceiver::FragmentTrailer*)(array+nEnd);
      fragmentTrailer->print();

      std::cout << std::endl;

      nEnd-=2;
      slinkEoe=(Hgcal10gLinkReceiver::SlinkEoe*)(array+nEnd);
      slinkEoe->print();

      std::cout << std::endl;

      nEnd-=2*slinkEoe->eventLength()-2;
      
      slinkBoe=(Hgcal10gLinkReceiver::SlinkBoe*)(array+nEnd);
      slinkBoe->print();

      vEvents[vEvents.size()-j-1].first=2*slinkEoe->eventLength();
      vEvents[vEvents.size()-j-1].second=(uint64_t*)slinkBoe;
      std::cout << std::endl;

      nEnd-=2;
    }
  }

  for(unsigned i(0);i<vEvents.size() && i<1;i++) {
    uint64_t *a(vEvents[i].second);
    for(unsigned j(0);j<vEvents[i].first;j++) {
      std::cout << "Word " << std::setw(6) << j << " = 0x"
		<< std::hex << std::setfill('0')
		<< std::setw(16) << a[j]
		<< std::dec << std::setfill(' ')
		<< std::endl;
    }
    std::cout << std::endl;
  }
  

  fout.close();
  
  /*  
  orbitHeader=(Hgcal10gLinkReceiver::OrbitHeader*)(array+506);
  orbitHeader->print();

  slinkBoe=(Hgcal10gLinkReceiver::SlinkBoe*)(array+4+506);
  slinkBoe->print();

  slinkEoe=(Hgcal10gLinkReceiver::SlinkEoe*)(array+502+506);
  slinkEoe->print();

  fragmentTrailer=(Hgcal10gLinkReceiver::FragmentTrailer*)(array+504+506);
  fragmentTrailer->print();
  */

  Hgcal10gLinkReceiver::OrbitReader oReader;
  unsigned nEvents(0);
  
  for(unsigned i(0);i<3 && oReader.open(std::string("dth_p2_test.raw_sourceid01230_index00")+std::to_string(i)+".raw");i++) {
    std::cout << std::endl << "FILE " << i << std::endl;
    
    //oReader.open("output-stream-1hz.raw");
    //oReader.open(std::string("dth_p2_test.raw_sourceid01230_index00")+std::to_string(i)+".raw");
    
    while(oReader.readOrbit(vEvents)) {
      std::cout << "Number of events in orbit = " << vEvents.size() << std::endl;
      
      for(unsigned i(0);i<vEvents.size();i++) {
	nEvents++;
	printEvent(vEvents[i].first,vEvents[i].second);
      }
    }
    
    oReader.close();
  }

  std::cout << "Number of events in run = " << nEvents << std::endl;

  return 0;
}
