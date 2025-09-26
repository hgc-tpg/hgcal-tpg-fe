#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>

#include "FileReader.h"
#include "EcontEnergies.h"
#include "TpgSubpacketHeader.h"

/*
void econtEnergies(const uint64_t *p, std::vector<uint16_t> &v) {
  v.resize(24);

  uint64_t e0,e1;
  
  e0=(p[0]<<32)|(p[1]&0xffffffff);
  e1=(p[1]<<48)|(p[2]&0xffffffff)<<16|(p[3]&0xffffffff)>>16;
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

unsigned unpackerEnergies(const uint64_t *p, std::vector<uint16_t> &v) {
  v.resize(24);
  for(unsigned i(0);i<6;i++) {
    v[2*i   ]=(p[i]    )&0x7f;
    v[2*i+ 1]=(p[i]>>13)&0x7f;
    v[2*i+12]=(p[i]>>32)&0x7f;
    v[2*i+13]=(p[i]>>45)&0x7f;
  }

  return (p[0]>>26)&0xf;
}
*/

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

  std::vector<uint16_t> v;

  // Create the file reader
  Hgcal10gLinkReceiver::FileReader _fileReader;

  // Make the buffer space for the records
  Hgcal10gLinkReceiver::RecordT<8*4095> *r(new Hgcal10gLinkReceiver::RecordT<8*4095>);

  // Set up specific records to interpet the formats
  const Hgcal10gLinkReceiver::RecordStarting *rStart((Hgcal10gLinkReceiver::RecordStarting*)r);
  const Hgcal10gLinkReceiver::RecordStopping *rStop ((Hgcal10gLinkReceiver::RecordStopping*)r);
  const Hgcal10gLinkReceiver::RecordRunning  *rEvent((Hgcal10gLinkReceiver::RecordRunning*) r);

  // Defaults to the files being in directory "dat"
  // Can call setDirectory("blah") to change this
  //_fileReader.setDirectory("somewhere/else");
  _fileReader.setDirectory(std::string("dat/Relay")+argv[1]);
  _fileReader.openRun(runNumber,0);
  //_fileReader.openRun(runNumber,1);
  
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
  uint32_t eventIdOld(0);

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
      bool print(nEvents<=100000);// || nEvents>2000);
      //anyError=false;

      // Check id is correct

	// Access the Slink header ("begin-of-event")
	// This should always be present; check pattern is correct
	const Hgcal10gLinkReceiver::SlinkBoe *b(rEvent->slinkBoe());
	assert(b!=nullptr);
	//if(b->l1aType()==0x001c) b->print();
      
	// Access the Slink trailer ("end-of-event")
	// This should always be present; check pattern is correct
	const Hgcal10gLinkReceiver::SlinkEoe *e(rEvent->slinkEoe());
	assert(e!=nullptr);

	if(b->eventId()>=0x441cea && b->eventId()<0x441cea+10) print=true;
	
	const uint64_t *p64(((const uint64_t*)rEvent)+1);
	unsigned nEcontWords(p64[3]&0xff);
	unsigned nUnpackerWords(p64[4+nEcontWords]&0xff);
	
	unsigned bx(Hgcal10gLinkReceiver::unpackerEnergies(p64+5+nEcontWords+6*((nUnpackerWords/6)/2),v));

	if(b->eventId()==eventIdOld || b->eventId()!=eventIdOld+1) {
	  std::cout << "ERROR: Old, new BX id = " << eventIdOld << ", " << b->eventId() << std::endl;
	}
	eventIdOld=b->eventId();
	
	if(print || bx==15 || e->bxId()==55) {
	std::cout << std::endl;
	rEvent->print();

	if(false) {
	for(unsigned u(0);u<(nEcontWords/4);u++) {
	  Hgcal10gLinkReceiver::econtEnergies(p64+4+4*u,v);
	  std::cout << "ECON-T energies " << u << std::endl;
	  for(unsigned j(0);j<12;j++) {
	    std::cout << " " << std::setw(3) << v[j];
	  }
	  std::cout << std::endl;
	  for(unsigned j(12);j<24;j++) {
	    std::cout << " " << std::setw(3) << v[j];
	  }
	  std::cout << std::endl;
	}
	for(unsigned u(0);u<(nUnpackerWords/6);u++) {
	  bx=Hgcal10gLinkReceiver::unpackerEnergies(p64+5+nEcontWords+6*u,v);
	  std::cout << "Unpacker energies " << u << ", bx = " << bx << std::endl;
	  for(unsigned j(0);j<12;j++) {
	    std::cout << " " << std::setw(3) << v[j];
	  }
	  std::cout << std::endl;
	  for(unsigned j(12);j<24;j++) {
	    std::cout << " " << std::setw(3) << v[j];
	  }
	  std::cout << std::endl;
	}
	}
	/*
	b->print();
	e->print();
      
	// Access the BE packet header
	const Hgcal10gLinkReceiver::BePacketHeader *bph(rEvent->bePacketHeader());
	//if(bph!=nullptr && print)
	bph->print();
	*/



	/*      
      // Access ECON-D packet as an array of 32-bit words
      const uint32_t *pEcond(rEvent->econdPayload());
      
      // Check this is not an empty event
      if(pEcond!=nullptr) {

	const uint64_t *p64(((const uint64_t*)rEvent)+1);
	bph=(const Hgcal10gLinkReceiver::BePacketHeader*)(p64+2);

	const uint32_t *p32(rEvent->daqPayload());

	
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
	  //unsigned bcEcond((p64[doubleEvents?3:2]>>20)&0xfff);
	  //unsigned ecEcond((p64[doubleEvents?3:2]>>14)&0x3f);
	  //unsigned ocEcond(((p64[doubleEvents?3:2]>>11)+ocOffset)&0x7);
	  unsigned bcEcond((p64[4]>>20)&0xfff);
	  unsigned ecEcond((p64[4]>>14)&0x3f);
	  unsigned ocEcond(((p64[4]>>11)+ocOffset)&0x7);
	  
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

	unsigned nSubPackets(0);
	for(unsigned i(0);i<2*rEvent->payloadLength()-4U;i++) {
	  if(p32[i]==0xe000001f) nSubPackets++;
	}
	hSubpacketCount->Fill(nSubPackets);
	hSubpacketCountVsPl->Fill(rEvent->payloadLength(),nSubPackets);
	
	if(p32[7]!=0xe000001f) {
	  hSubpackets->Fill(0);
	  nErrors[7]++;
	  print=true;
	  anyError=true;
	  
	} else if(p32[44]!=0xe000001f) {
	  hSubpackets->Fill(1);
	  nErrors[7]++;
	  print=true;
	  anyError=true;

	} else if(p32[85]!=0xe000001f) {
	  hSubpackets->Fill(2);
	  nErrors[7]++;
	  print=true;
	  anyError=true;

	} else if(p32[122]!=0xe000001f) {
	  hSubpackets->Fill(3);
	  nErrors[7]++;
	  print=true;
	  anyError=true;

	} else if(p32[163]!=0xe000001f) {
	  hSubpackets->Fill(4);
	  nErrors[7]++;
	  print=true;
	  anyError=true;

	} else if(p32[200]!=0xe000001f) {
	  hSubpackets->Fill(5);
	  nErrors[7]++;
	  print=true;
	  anyError=true;
	}

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

	*/	
	if(print) {
	  for(unsigned i(0);i<rEvent->payloadLength();i++) {
	    if((p64[i]>>32)==0xcafecafe) {
	      std::cout << std::endl;
	      const Hgcal10gLinkReceiver::TpgSubpacketHeader *t((Hgcal10gLinkReceiver::TpgSubpacketHeader*)(p64+i));
	      t->print();
	    }
	  }

	  std::cout << "Record " << nEvents << std::endl;
	  //rEvent->RecordHeader::print();

	  for(unsigned i(0);i<rEvent->payloadLength();i++) {
	    std::cout << "Word " << std::setw(3) << i << " ";
	    std::cout << std::hex << std::setfill('0');
	    std::cout << "0x" << std::setw(16) << p64[i] << std::endl;
	    std::cout << std::dec << std::setfill(' ');
	  }
	  std::cout << std::endl;
	  /*
	  for(unsigned i(0);i<2*rEvent->payloadLength()-4U;i++) {
	    std::cout << "Word " << std::setw(3) << i << " ";
	    std::cout << std::hex << std::setfill('0');
	    std::cout << "0x" << std::setw(8) << p32[i] << std::endl;
	    std::cout << std::dec << std::setfill(' ');
	  }

	  std::cout << std::endl;
	  */
	}

	// Do other stuff here

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

  return 0;
}
