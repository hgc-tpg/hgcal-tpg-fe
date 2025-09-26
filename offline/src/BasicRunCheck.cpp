#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <queue>
#include <stdlib.h>
#include <cassert>
#include <sys/types.h>
#include <getopt.h>

#include "FileReader.h"
#include "RecordPrinter.h"

using namespace Hgcal10gLinkReceiver;

int main(int argc, char** argv) {
  if(argc<3) {
    std::cerr << argv[0] << ": no run and/or link number specified" << std::endl;
    return 1;
  }

  bool doSrun(false);

  unsigned runNumber(0);
  unsigned linkNumber(0);
  /*  
  unsigned fileNumber(0);
  unsigned packetNumber(0);
  unsigned i(0),j(0);
  
  uint64_t initialLF(0);
  const uint64_t maskLF(0xff00ffffffffffff);
  */
  
  // Handle run number
  std::istringstream issRun(argv[1]);
  issRun >> runNumber;
  std::istringstream issLink(argv[2]);
  issLink >> linkNumber;

  FsmState::State state(doSrun?FsmState::Halted:FsmState::ConfiguredB);
  FsmState::State previousState(FsmState::EndOfStateEnum);

  FileReader _fileReader;
  RecordT<1024> h;

  unsigned recordCount(0);
  unsigned recordQueueDepth(5);
  std::queue< RecordT<1024> > recordQueue;
  
  _fileReader.open(runNumber,linkNumber,doSrun);
  //_fileReader.open(argv[1]);

  while(_fileReader.read(&h)) {
    recordCount++;
    
    recordQueue.push(h);
    if(recordQueue.size()>recordQueueDepth) recordQueue.pop();

    if(previousState!=state) {
      std::cout << std::endl << "***** Previous state = "
		<< FsmState::stateName(state)
		<< " *****" << std::endl << std::endl;
      previousState=state;
    }
    
    if(h.state()!=FsmState::Running) RecordPrinter(&h);

    if(!h.validPattern() || !h.validState()) {
      std::cout << std::endl << "BAD PATTERN OR STATE AT RECORD COUNT " << recordCount
		<< ", DUMPING LAST " << recordQueueDepth
		<< " RECORDS" << std::endl << std::endl;
      
      while(recordQueue.size()>0) {
	RecordPrinter(recordQueue.front());
	recordQueue.pop();
      }
      //assert(false);
    }
    /*    
    if(FsmState::staticState(state)) {
      if(state!=h.state()) {
	assert(state==FsmState::staticStateBeforeTransient(h.state()));
	state=FsmState::staticStateAfterTransient(h.state());
      }
    } else {
      assert(h.state()==FsmState::staticStateBeforeTransient(state));
      state=h.state();      
    }
    */
  }

  std::cout << "Total record count = " << recordCount << std::endl;
  return 0;
}
