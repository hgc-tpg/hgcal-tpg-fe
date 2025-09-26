#include <iostream>
#include <iomanip>
#include <cassert>

#include "FileReader.h"

int main(int argc, char** argv) {
  if(argc<2) {
    std::cerr << argv[0] << ": no relay and/or run numbers specified" << std::endl;
    return 1;
  }

  // Handle run number
  unsigned relayNumber(0);
  std::istringstream issRelay(argv[1]);
  issRelay >> relayNumber;

  unsigned runNumber(relayNumber);
  if(argc>2) {
    std::istringstream issRun(argv[2]);
    issRun >> runNumber;
  }
  
  if(relayNumber==0 || runNumber==0) {
    std::cerr << argv[0] << ": relay and/or run numbers uninterpretable" << std::endl;
    return 2;
  }

  // Create the file reader
  Hgcal10gLinkReceiver::FileReader _fileReader;

  // Make the buffer space for the records
  Hgcal10gLinkReceiver::RecordT<4095> *r(new Hgcal10gLinkReceiver::RecordT<4095>);

  // Set up specific records to interpet the formats
  const Hgcal10gLinkReceiver::RecordStarting *rStart((Hgcal10gLinkReceiver::RecordStarting*)r);
  const Hgcal10gLinkReceiver::RecordStopping *rStop ((Hgcal10gLinkReceiver::RecordStopping*)r);
  const Hgcal10gLinkReceiver::RecordRunning  *rEvent((Hgcal10gLinkReceiver::RecordRunning*) r);

  // Defaults to the files being in directory "dat"
  // Can call setDirectory("blah") to change this
  //_fileReader.setDirectory("somewhere/else");
  _fileReader.setDirectory(std::string("dat/Relay")+argv[1]);
  //_fileReader.openRun(runNumber,0);
  _fileReader.openRun(runNumber,1);

  std::ostringstream fName;
  fName << "Run" << runNumber << "_Counters.txt";
  std::ofstream fout;
  fout.open(fName.str().c_str());
  
  //bool anyError(false);
  unsigned nEvents(0);//,nErrors[10]={0,0,0,0,0,0,0,0,0,0};
  /*
  unsigned bc(0),oc(0),ec(0),seq(0);
  unsigned ocOffset(0);

  bool noCheck(false);
  uint32_t initialSeq;
  
  bool doubleEvents(false);
  unsigned nEventsPerOrbit(1);
  */
  while(_fileReader.read(r)) {
    if(       r->state()==Hgcal10gLinkReceiver::FsmState::Starting) {
      rStart->print();
      std::cout << std::endl;

    } else if(r->state()==Hgcal10gLinkReceiver::FsmState::Stopping) {
      rStop->print();
      std::cout << std::endl;

    } else {

	// We have an event record
      nEvents++;
      bool print(nEvents<=100);// || nEvents>2000);
      //anyError=false;

      // Check id is correct
      
      // Access the Slink header ("begin-of-event")
      // This should always be present; check pattern is correct
      const Hgcal10gLinkReceiver::SlinkBoe *b(rEvent->slinkBoe());
      assert(b!=nullptr);
      //if(b->l1aType()==0x001c) b->print();
      fout << b->eventId() << ",";
      
      
      // Access the Slink trailer ("end-of-event")
      // This should always be present; check pattern is correct
      const Hgcal10gLinkReceiver::SlinkEoe *e(rEvent->slinkEoe());
      assert(e!=nullptr);
      fout << e->bxId() << "," << e->orbitId() << ",";
      
      // Access the BE packet header
      const Hgcal10gLinkReceiver::BePacketHeader *bph(rEvent->bePacketHeader());
      assert(bph!=nullptr);

      fout << unsigned(bph->eventCounter()) << ","
	   << unsigned(bph->bunchCounter()) << ","
	   << unsigned(bph->orbitCounter()) << ",";

      if(rEvent->payloadLength()>6) {
	const uint64_t *p64(((const uint64_t*)rEvent)+4);

	unsigned bcEcond((p64[4]>>20)&0xfff);
	unsigned ecEcond((p64[4]>>14)&0x3f);
	unsigned ocEcond((p64[4]>>11)&0x7);
	

	fout << ecEcond << ","
	     << bcEcond << ","
	     << ocEcond << std::endl;

      } else {
	fout << "9999,9999,9999" << std::endl;
      }
    }
      
  }

  std::cout << "Total number of event records seen = "
	    << nEvents << std::endl;
  /*
  for(unsigned i(0);i<10;i++) {
    std::cout << "Total number of event records with error " << i << " = "
	      << nErrors[i] << std::endl;
  }
  std::cout << "Final OC offset = " << ocOffset << std::endl;
  */

  delete r;

  fout.close();
  
  return 0;
}
