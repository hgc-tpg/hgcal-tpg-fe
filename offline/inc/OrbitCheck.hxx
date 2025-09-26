#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <thread>

#include "OrbitReader.h"

int main(int argc, char *argv[]) {
  uint32_t runNumber(0xffffffff);
  uint32_t sourceId(0xffffffff);
  uint32_t firstLs(1);
  if(argc<3) return 1;

  std::istringstream issrn(argv[1]);
  issrn >> runNumber;
  std::istringstream issid(argv[2]);
  issid >> sourceId;
  if(argc>=4) {
    std::istringstream issls(argv[3]);
    issls >> firstLs;
  }
  
  std::vector<Hgcal10gLinkReceiver::OrbitReaderEvent> vEvents;

  OrbitCheckTypedef ct;
  assert(ct.runStart(runNumber,sourceId));
  
  std::string oDir("OrbitData/");

  Hgcal10gLinkReceiver::OrbitReader oReader;
  if(argc>3) oReader.setPrint(true);
  
  unsigned nEvents(0);

  const Hgcal10gLinkReceiver::OrbitHeader *oh;
  Hgcal10gLinkReceiver::FragmentTrailer *ft;
  
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

    while((oh=oReader.readOrbit(vEvents))!=nullptr) {
      assert(ct.orbit(*oh));
      
      for(unsigned j(0);j<vEvents.size();j++) {
	assert(ct.event(vEvents[j]));
	nEvents++;
      }
    }
    
    oReader.close();
  }

  std::cout << "Total number of events in run = " << nEvents << std::endl;

  assert(ct.runStop());

  return 0;
}
