#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <cassert>
#include <sys/types.h>
#include <getopt.h>

#include "TFileHandler.h"
#include "TH1F.h"

#include "FileReader.h"
#include "RecordPrinter.h"

using namespace Hgcal10gLinkReceiver;

int main(int argc, char** argv) {
  if(argc<2) {
    std::cerr << argv[0] << ": no run number specified" << std::endl;
    return 1;
  }

  TFileHandler tfh(std::string("BasicSuperRunCheck")+argv[1]);
  TH1I *h0=new TH1I("EventsVsTime0",";Time;Events/sec",50000,0,50000);
  TH1F *h1=new TH1F("EventsVsTime1",";Time;Events/sec", 5000,0,50000);
  
  bool doSrun(true);

  unsigned superRunNumber(0);
  
  // Handle run number
  std::istringstream issRun(argv[1]);
  issRun >> superRunNumber;

  FsmState::State state(doSrun?FsmState::Halted:FsmState::ConfiguredB);
  
  FileReader _fileReader;
  RecordT<1024> h;

  RecordConfiguringA &rca((RecordConfiguringA&)h);
  RecordConfiguringB &rcb((RecordConfiguringB&)h);
  RecordStarting     &rst((RecordStarting&    )h);
  RecordStopping     &rsp((RecordStopping&    )h);
  RecordHaltingB     &rhb((RecordHaltingB&    )h);
  RecordHaltingA     &rha((RecordHaltingA&    )h);
  
  _fileReader.openRelay(superRunNumber);

  unsigned nCfgA(0);
  unsigned nCfgB(0);
  unsigned nSt(0);
  unsigned nPs(0);
  unsigned nRs(0);
  unsigned nSp(0);

  unsigned nRunSt(0);
  unsigned nRunPs(0);
  unsigned nRunEv(0);
  unsigned nRunRs(0);
  unsigned nRunSp(0);
  unsigned nSeqOff;
  
  unsigned nCfgInSrun(0);
  unsigned nRunInCfg(0);
  unsigned nRunInSrun(0);
  unsigned nEvtInRun(0);
  unsigned nEvtInCfg(0);
  unsigned nEvtInSrun(0);
      
  
  RecordHeader previous;
  
  while(_fileReader.read(&h)) {
    std::cout << std::endl << "***** Previous state = "
	      << FsmState::stateName(state)
	      << " *****" << std::endl;
    
    if(h.state()!=FsmState::Running) RecordPrinter(&h);
  
    if(FsmState::staticState(state)) {
      if(state!=h.state()) {
	assert(state==FsmState::staticStateBeforeTransient(h.state()));
	state=FsmState::staticStateAfterTransient(h.state());
      }
    } else {
      assert(h.state()==FsmState::staticStateBeforeTransient(state));
      state=h.state();      
    }

    if(h.state()==FsmState::ConfiguringA) {
      assert(rca.valid());
      assert(rca.utc()==superRunNumber);
      assert(rca.superRunNumber()==superRunNumber);

      assert(nCfgA==0);
      nCfgA++;
    }

    if(h.state()==FsmState::ConfiguringB) {
      //assert(rcb.valid());
      assert(rcb.utc()>=previous.utc());
      assert(rcb.superRunNumber()==superRunNumber);

      nCfgB++;
      assert(rcb.configurationCounter()==nCfgB);
      nSt=0;
      nSp=0;

      nCfgInSrun++;
      nRunInCfg=0;
      nEvtInCfg=0;
    }
    
    if(h.state()==FsmState::Starting) {
      //assert(rst.valid());
      assert(rst.utc()>=previous.utc());
      assert(rst.runNumber()==rst.utc());

      nSt++;
      nPs=0;
      nRs=0;
      
      nRunInSrun++;
      nRunInCfg++;
      nEvtInRun=0;

      RecordT<1024> hRun;
      RecordStarting     &runRst((RecordStarting&    )hRun);
      RecordPausing      &runRps((RecordPausing&     )hRun);
      RecordRunning      &runRrn((RecordRunning&     )hRun);
      RecordResuming     &runRrs((RecordResuming&    )hRun);
      RecordStopping     &runRsp((RecordStopping&    )hRun);

      FileReader _fileReaderRun;
      _fileReaderRun.openRun(rst.runNumber(),0);
      //if(_fileReaderRun.closed()) return 1;
	 
      nRunSt=0;
      nRunPs=0;
      nRunEv=0;
      nRunRs=0;
      nRunSp=0;
      
      std::cout << "RUN FILE!" << std::endl;
      while(_fileReaderRun.read(&hRun)) {
	if(hRun.state()!=FsmState::Running) {
	  RecordPrinter(&hRun);

	  if(hRun.state()==FsmState::Starting) {
	    uint64_t *a((uint64_t*)&rst);
	    uint64_t *b((uint64_t*)&runRst);
	    for(unsigned i(0);i<1+rst.payloadLength();i++) {
	      std::cout << i << ", " << a[i] << ", " << b[i] << std::endl;
	      assert(a[i]==b[i]);
	    }
	    assert(nRunSt==0);
	    nRunSt++;

	    
	  } else if(hRun.state()==FsmState::Pausing) {
	    assert(nRunPs==nRunRs);
	    nRunPs++;
	    
	  } else if(hRun.state()==FsmState::Resuming) {
	    nRunRs++;
	    assert(nRunPs==nRunRs);
	    
	  } else if(hRun.state()==FsmState::Stopping) {
	    //assert(runRsp.numberOfEvents()<=rst.maxEvents());
	    assert(runRsp.numberOfSeconds()<=rst.maxSeconds());
	    assert(runRsp.numberOfSpills()<=rst.maxSpills());

	    assert(nRunPs==runRsp.numberOfPauses());
	    assert(nRunRs==runRsp.numberOfPauses());


	    assert(nRunEv<=runRsp.numberOfEvents());
	    if(nRunEv<runRsp.numberOfEvents()) {
	      std::cerr << "Events dropped = "
			<< runRsp.numberOfEvents() << " - "
			<< nRunEv << std::endl;
	    }
	    
	    assert(nRunSp==0);
	    nRunSp++;

	  } else {
	    assert(false);
	  }
	} else {
	  if(nRunEv==0) nSeqOff=runRrn.sequenceCounter();
	  if(nRunEv<10) runRrn.print();

	  if(runRrn.sequenceCounter()!=nSeqOff+nRunEv) {
	    std::cerr << "Skip in event sequence counter from "
		      << nSeqOff+nRunEv << " expected to "
		      << runRrn.sequenceCounter() << " seen" << std::endl;
	    nSeqOff=runRrn.sequenceCounter()-nRunEv;
	  }

	  h0->Fill(runRrn.slinkEoe()->orbitId()-superRunNumber);
	  h1->Fill(runRrn.slinkEoe()->orbitId()-superRunNumber,0.1);
	  nRunEv++;

	  nEvtInSrun++;
	  nEvtInCfg++;
	  nEvtInRun++;
	}
      }
      std::cout << "nRunSt = " << nRunSt << ", nRunSp = " << nRunSp
		<< ", nRunPs = " << nRunPs << ", nRunRs = " << nRunRs
		<< ", nRunEv = " << nRunEv << std::endl;
      std::cout << "RUN FILE DONE!" << std::endl;

      _fileReaderRun.close();
    }

    if(h.state()==FsmState::Stopping) {
      //assert(rsp.valid());
      assert(rsp.utc()>=previous.utc());
      assert(rsp.runNumber()==rst.runNumber());

      assert(rsp.numberOfEvents()<=rst.maxEvents());
      assert(rsp.numberOfSeconds()<=rst.maxSeconds());
      assert(rsp.numberOfSpills()<=rst.maxSpills());

      assert(nRunPs==rsp.numberOfPauses());
      assert(nRunRs==rsp.numberOfPauses());

      //assert(rsp.numberOfEvents()+1==nRunEv);
      if(nRunEv<rsp.numberOfEvents()) {
	std::cerr << "Events dropped = "
		  << rsp.numberOfEvents() << " - "
		  << nRunEv << " = " << rsp.numberOfEvents()+1-nRunEv
		  << std::endl;
	    }

      nSp++;
      assert(nSp==nSt);
      assert(nPs==nRs);
    }

    if(h.state()==FsmState::HaltingB) {
      assert(rhb.valid());
      assert(rhb.utc()>=previous.utc());
      assert(rhb.superRunNumber()==superRunNumber);

      assert(rhb.configurationNumber()==nCfgInSrun);
      assert(rhb.numberOfRuns()==nRunInCfg);
      //assert(rhb.numberOfEvents()==nEvtInCfg);
      assert(rhb.numberOfRuns()==nSt);
      assert(rhb.numberOfRuns()==nSp);
    }


    if(h.state()==FsmState::HaltingA) {
      assert(rha.valid());
      assert(rha.utc()>=previous.utc());
      assert(rha.superRunNumber()==superRunNumber);
      assert(rha.numberOfConfigurations()==nCfgInSrun);
      assert(rha.numberOfRuns()==nRunInSrun);
      assert(rha.numberOfEvents()==nEvtInSrun);
    }

    previous=h;
  }
  return 0;
}

#ifdef JUNK
  /*
  
  std::ostringstream ossRun;
  ossRun << std::setfill('0') << std::setw(9) << runNumber;
  
  std::string sRun("Run"); // FIX!!!
  sRun+=ossRun.str()+"_File";
  




  bool doPrint(false);
  */

  uint64_t p64[512*16];

  uint64_t predictedPH(0);
  uint64_t maskPH(0xff00ffffffffffff);

  bool endOfRun(false);
  
  while(!endOfRun) {
    /*
    // Handle file number
    std::ostringstream ossFile;
    ossFile << std::setfill('0') << std::setw(9) << fileNumber;
    */
    //std::cout << "Directly reading " << sRun+ossFile.str()+".bin" << std::endl;

    // Open file
    std::ifstream fin;

    // tellg, seekg?
    
    unsigned p(0);
    unsigned length;
    bool endOfFile(false);

    // Read starting packet
    bool startOfFile(true);
    
    while(!endOfFile && !endOfRun) {
      if(startOfFile) {
	fin.open(RunFileName(runNumber,linkNumber,fileNumber).c_str(),std::ios::binary);
	if(!fin) {
	  std::cerr << argv[0] << ": cannot open file " << RunFileName(runNumber,0,fileNumber) << std::endl;
	  return 2;
	}
      }
      
      assert(fin);

      fin.read((char*)p64,8);
      if(startOfFile) {
	if(fileNumber==0) assert((p64[0]>>56)==0xbb);
	else              assert((p64[0]>>56)==0xdd);
      }
      
      if((p64[0]>>56)!=0xac) {
	std::cout << "  Word " << std::setw(10) << j << ", data = 0x"
		  << std::hex << std::setfill('0')
		  << std::setw(16) << p64[0]
		  << std::dec << std::setfill(' ')
		  << std::endl;
      } else {
	if(predictedPH==0) {
	  predictedPH=p64[0]+1;
	  
	} else {
	  if((p64[0]&maskPH)!=(predictedPH&maskPH)) {
	    std::cout << "  Packet header " << std::setw(10) << j << ", data = 0x"
		      << std::hex << std::setfill('0')
		      << std::setw(16) << p64[0] << ", predicted = 0x"
		      << std::setw(16) << predictedPH
		      << std::dec << std::setfill(' ');

	    if((p64[0]&maskPH)>(predictedPH&maskPH)) {
	      std::cout << ", difference = " << (p64[0]&maskPH)-(predictedPH&maskPH) << std::endl;
	    } else {
	      std::cout << ", difference = -" << (predictedPH&maskPH)-(p64[0]&maskPH) << std::endl;
	    }
	  }
	  predictedPH=p64[0]+1;
	}
      }
      
      j++;
      length=(p64[0]>>48)&0xff;
      fin.read((char*)(p64+1),8*length);
      j+=length;

      startOfFile=false;

      if((p64[0]>>56)==0xcc) {
	fin.close();

	startOfFile=true;
	fileNumber++;
      } else {
      }
      
      if((p64[0]>>56)==0xee) {
	fin.close();
	endOfRun=true;
      }

    }

  }


/*
  if(fileNumber==0) assert((p64[0]>>56)==0xbb);
    else              assert((p64[0]>>56)==0xcc);
    
    
    fin.read((char*)(p64+1),8);
    j++;
    
    std::cout << "Run " << (p64[1]>>32)
	      << ", UTC " << (p64[1]&0xffffffff)
	      << ", file number = " << fileNumber
	      << std::endl;
    
    // Events
    
    fin.read((char*)p64,8);
    j++;
    
    if(fileNumber==0) initialLF=p64[0];
    assert((p64[0]>>56)==0xac);
    
    while(fin && !endOfFile && !endOfRun) {
      if(doPrint) std::cout << "New  Word " << std::setw(6) << 0 << ", data = 0x"
			    << std::hex << std::setfill('0')
			    << std::setw(16) << p64[0]
			    << std::dec << std::setfill(' ')
			    << std::endl;
      
      assert((p64[0]>>56)==0xac);
      
      if((p64[0]&maskLF)!=((initialLF+packetNumber)&maskLF)) {
	std::cout << "Bad  Word " << std::setw(6) << 0 << ", data = 0x"
		  << std::hex << std::setfill('0')
		  << std::setw(16) << p64[0] << ", expected = 0x"
		  << std::setw(16) << initialLF+packetNumber
		  << std::dec << std::setfill(' ')
		  << ", difference = " << (p64[0]&maskLF)-((initialLF+packetNumber)&maskLF)
		  << std::endl;
	
	initialLF=p64[0]-packetNumber;
      }
      
      //assert((p64[0]&maskLF)==((initialLF+p)&maskLF));
      packetNumber++;

      unsigned length((p64[0]>>48)&0xff);
  
      std::cout << "Length = " << length << " words, "
		<< "  Word " << std::setw(6) << 0 << ", data = 0x"
		<< std::hex << std::setfill('0')
		<< std::setw(16) << p64[0]
		<< std::dec << std::setfill(' ')
		<< std::endl;
      
      fin.read((char*)(p64+1),8*length);
      j+=length;
      
      if(doPrint || j< 500) {
	for(unsigned i(1);i<=length;i++) {
	  std::cout << "  Word " << std::setw(6) << i << ", data = 0x"
		    << std::hex << std::setfill('0')
		    << std::setw(16) << p64[i]
		    << std::dec << std::setfill(' ')
		    << std::endl;
	}
      }
      
      fin.read((char*)p64,8);
      j++;
      if((p64[0]>>56)==0xcc) endOfFile=true;
      if((p64[0]>>56)==0xee) endOfRun=true;
    }
    
    length=(p64[0]>>48)&0xff;
    assert(length==1);
    
    fin.read((char*)(p64+1),8);
    j++;
    
    std::cout << "Run " << (p64[1]>>32)
	      << ", UTC " << (p64[1]&0xffffffff)
	      << std::endl;
    
    fin.close();
    fileNumber++;
  }

  std::cout << "Directly read " << j << " words" << std::endl;
  */
#endif
