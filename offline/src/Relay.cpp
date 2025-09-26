/*
g++ -Wall -Icommon/inc src/Relay.cpp -lyaml-cpp -o bin/Relay.exe
*/

#include <iostream>
#include <iomanip>
#include <cassert>

#include <yaml-cpp/yaml.h>

#include "RelayReader.h"

void econtEnergies(const uint64_t *p, std::vector<uint16_t> &v) {
  v.resize(24);

  uint64_t e0,e1;
  
  e0=p[0]<<32|(p[1]&0xffffffff);
  e1=p[1]<<48|(p[2]&0xffffffff)<<16|(p[3]&0xffffffff)>>16;
  //std::cout << "e0,e1 = " << std::hex << e0 << "," << e1 << std::endl;
  
  v[ 0]=(e0>> 1)&0x7f;
  v[ 1]=(e0>> 8)&0x7f;
  v[ 2]=(e0>>15)&0x7f;
  v[ 3]=(e0>>22)&0x7f;
  v[ 4]=(e0>>29)&0x7f;
  
  v[ 5]=(e1    )&0x7f;
  v[ 6]=(e1>> 7)&0x7f;
  v[ 7]=(e1>>14)&0x7f;
  v[ 8]=(e1>>21)&0x7f;
  v[ 9]=(e1>>28)&0x7f;
  v[10]=(e1>>35)&0x7f;
  v[11]=(e1>>42)&0x7f;
  
  e0=(p[0]&0xffffffff00000000)    |(p[1]&0xffffffff00000000)>>32;
  e1=(p[1]&0xffffffff00000000)<<16|(p[2]&0xffffffff00000000)>>16|(p[3]&0xffffffff00000000)>>48;
  //std::cout << "e0,e1 = " << std::hex << e0 << "," << e1 << std::endl;
  
  v[12]=(e0>> 1)&0x7f;
  v[13]=(e0>> 8)&0x7f;
  v[14]=(e0>>15)&0x7f;
  v[15]=(e0>>22)&0x7f;
  v[16]=(e0>>29)&0x7f;
  
  v[17]=(e1    )&0x7f;
  v[18]=(e1>> 7)&0x7f;
  v[19]=(e1>>14)&0x7f;
  v[20]=(e1>>21)&0x7f;
  v[21]=(e1>>28)&0x7f;
  v[22]=(e1>>35)&0x7f;
  v[23]=(e1>>42)&0x7f;
}

void unpackerEnergies(const uint64_t *p, std::vector<uint16_t> &v) {
  v.resize(24);

  for(unsigned i(0);i<6;i++) {
    v[2*i   ]=(p[i]    )&0x7f;
    v[2*i+ 1]=(p[i]>>13)&0x7f;
    v[2*i+12]=(p[i]>>32)&0x7f;
    v[2*i+13]=(p[i]>>45)&0x7f;
  }
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

  // Create the file reader
  Hgcal10gLinkReceiver::RelayReader _relayReader;
  _relayReader.enableLink(0,true);
  _relayReader.enableLink(1,false);
  _relayReader.openRelay(relayNumber);
  
  // Make the buffer space for the records and some useful casts for configuration and event records
  Hgcal10gLinkReceiver::RecordT<4*4095> *rxxx(new Hgcal10gLinkReceiver::RecordT<4*4095>);
  Hgcal10gLinkReceiver::RecordYaml *rCfg((Hgcal10gLinkReceiver::RecordYaml*)rxxx);
  Hgcal10gLinkReceiver::RecordRunning *rEvent((Hgcal10gLinkReceiver::RecordRunning*)rxxx);

  unsigned nEvents(0);
  
  while(_relayReader.read(rxxx)) {
    if(rxxx->state()!=Hgcal10gLinkReceiver::FsmState::Running) {
      
      if(rxxx->state()==Hgcal10gLinkReceiver::FsmState::Configuration) {
	YAML::Node n(YAML::Load(rCfg->string()));
	
	if(n["Source"].as<std::string>()=="Serenity") {
	}
	if(n["Source"].as<std::string>()=="TCDS2") {
	}
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
      
      if(nEvents<=10) {
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
  
  delete rxxx;
  
  return 0;
}
