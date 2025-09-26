#include <iostream>
#include <iomanip>
#include <cassert>

#include "TFileHandler.h"

#include "FileReader.h"

int main(int argc, char** argv) {
  if(argc<3) {
    std::cerr << argv[0] << ": no relay and/or run numbers specified" << std::endl;
    return 1;
  }

  // Handle run number
  unsigned relayNumber(0);
  unsigned runNumber(0);

  std::istringstream issRelay(argv[1]);
  issRelay >> relayNumber;

  std::istringstream issRun(argv[2]);
  issRun >> runNumber;

  if(relayNumber==0 || runNumber==0) {
    std::cerr << argv[0] << ": relay and/or run numbers uninterpretable" << std::endl;
    return 2;
  }

  TFileHandler tfh("Hacked");
  
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
  _fileReader.openRun(runNumber,1);
  
  bool anyError(false);
  unsigned nEvents(0),nErrors[10]={0,0,0,0,0,0,0,0,0,0};
  unsigned bc(0),oc(0),ec(0),seq(0);
  unsigned ocOffset(0);
  
  bool doubleEvents(false);
  unsigned nEventsPerOrbit(1);
  
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
      bool print(nEvents<=3);
      anyError=false;
      /*
      if(print) {
	rEvent->print();
	std::cout << std::endl;
      }
      */

      if(nEvents==1) seq=rEvent->sequenceCounter();

      if(seq!=rEvent->sequenceCounter()) {
	std::cout << "Event " << nEvents << " Sequence error; seen = " << rEvent->sequenceCounter()
		  << ", expected = " << seq << ", difference = "
		  << rEvent->sequenceCounter()-seq << std::endl;
	  seq=rEvent->sequenceCounter();
	  print=true;
	  anyError=true;
	  nErrors[8]++;
      }
      seq++;

      // Check id is correct
      if(!rEvent->valid()) rEvent->print();

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
      
      // Access the BE packet header
      const Hgcal10gLinkReceiver::BePacketHeader *bph(rEvent->bePacketHeader());
      //if(bph!=nullptr && print) bph->print();
      
      // Access ECON-D packet as an array of 32-bit words
      const uint32_t *pEcond(rEvent->econdPayload());
      
      // Check this is not an empty event
      if(pEcond!=nullptr) {

	const uint64_t *p64(((const uint64_t*)rEvent)+1);
	bph=(const Hgcal10gLinkReceiver::BePacketHeader*)p64;

	if(nEvents==1) {
	  bc=bph->bunchCounter();
	  oc=bph->orbitCounter();
	  ec=bph->eventCounter();
	} 

	if(bph->bunchCounter()!=(bc+200*((nEvents-1)%nEventsPerOrbit))) {
	    std::cout << "Event " << nEvents << " BC error; seen = " << bph->bunchCounter()
		      << ", expected = " << bc << std::endl;
	    bc=bph->bunchCounter();
	    print=true;
	    anyError=true;
	    nErrors[0]++;
	  }
	
	  if(bph->orbitCounter()!=(oc%8)) {
	    std::cout << "Event " << nEvents << " OC error; seen = " << unsigned(bph->orbitCounter())
		      << ", expected = " << unsigned(oc%8) << std::endl;
	    oc=bph->orbitCounter();
	    print=true;
	    anyError=true;
	    nErrors[1]++;
	  }
	  
	  if(bph->eventCounter()!=(ec%64)) {
	    std::cout << "Event " << nEvents << " EC error; seen = " << unsigned(bph->eventCounter())
		      << ", expected = " << unsigned(ec%64) << std::endl;
	    ec=bph->eventCounter();
	    print=true;
	    anyError=true;
	    nErrors[2]++;
	  }
	  
	if(rEvent->payloadLength()>2) {
	  unsigned bcEcond((p64[doubleEvents?3:2]>>20)&0xfff);
	  unsigned ecEcond((p64[doubleEvents?3:2]>>14)&0x3f);
	  unsigned ocEcond(((p64[doubleEvents?3:2]>>11)+ocOffset)&0x7);
	  
	  if(bcEcond!=bc) {
	    std::cout << "Event " << nEvents << " ECOND BC error; seen = " << bcEcond
		      << ", expected = " << bc << std::endl;
	    //bc=bcEcond;
	    print=true;
	    anyError=true;
	    nErrors[3]++;
	  }
	  
	  if(ocEcond!=(oc%8)) {
	    std::cout << "Event " << nEvents << " ECOND OC offset " << ocOffset << ", error; seen = " << ocEcond
		      << ", expected = " << unsigned(oc%8) << std::endl;
	    if(ocEcond>(oc%8)) ocOffset+=8-ocEcond+(oc%8);
	    else ocOffset+=(oc%8)-ocEcond;
	    ocOffset=(ocOffset%8);
	    print=true;
	    anyError=true;
	    nErrors[4]++;
	  }
	  
	  if(ecEcond!=(ec%64)) {
	    std::cout << "Event " << nEvents << " ECOND EC error; seen = " << ecEcond
		      << ", expected = " << unsigned(ec%64) << std::endl;
	    //ec=ecEcond;
	    print=true;
	    anyError=true;
	    nErrors[5]++;
	  }
	}
	
	if(print) bph->print();

	


	if((nEvents%nEventsPerOrbit)==0) oc++;
	ec++;

	if(doubleEvents) {
	bph=(const Hgcal10gLinkReceiver::BePacketHeader*)(p64+122);
	for(unsigned i(0);i<121;i++) {
	  if(p64[i]==0xcdcdcdcdcdcdcdcd) {
	    std::cout << "Event " << nEvents << " ERROR divider found with i = " << i << std::endl;
	    bph=(const Hgcal10gLinkReceiver::BePacketHeader*)(p64+i+1);
	    print=true;
	    anyError=true;
	  }
	}
	
        if(bph->bunchCounter()!=bc) {
          std::cout << "Event " << nEvents << " BC error; seen = " << bph->bunchCounter()
                    << ", expected = " << bc << std::endl;
          bc=bph->bunchCounter();
          print=true;
	  anyError=true;
        }

	if(bph->orbitCounter()!=(oc%8)) {
	  std::cout << "Event " << nEvents << " OC error; seen = " << unsigned(bph->orbitCounter())
		    << ", expected = " << unsigned(oc%8) << std::endl;
	  oc=bph->orbitCounter();
	  print=true;
	  anyError=true;
	}
	
	if(bph->eventCounter()!=(ec%64)) {
	  std::cout << "Event " << nEvents << " EC error; seen = " << unsigned(bph->eventCounter())
		    << ", expected = " << unsigned(ec%64) << std::endl;
          ec=bph->eventCounter();
          print=true;
	  anyError=true;
        }

        if(print) bph->print();

        oc++;
        ec++;
	}
	if(anyError) nErrors[9]++;

	
	if(print && nEvents<100) {
	  std::cout << "Header and words of ECON-D packet " << nEvents << std::endl;
	  rEvent->RecordHeader::print();

	  for(unsigned i(0);i<rEvent->payloadLength();i++) {
	    std::cout << "Word " << std::setw(3) << i << " ";
	    std::cout << std::hex << std::setfill('0');
	    std::cout << "0x" << std::setw(16) << p64[i] << std::endl;
	    std::cout << std::dec << std::setfill(' ');
	  }
	  std::cout << std::endl;
	}

	// Do other stuff here

      }

      
    }
  }

  std::cout << "Total number of event records seen = "
	    << nEvents << std::endl;
  for(unsigned i(0);i<10;i++) {
    std::cout << "Total number of event records with error " << i << " = "
	      << nErrors[i] << std::endl;
  }
  std::cout << "Final OC offset = " << ocOffset << std::endl;


  delete r;

  return 0;
}
