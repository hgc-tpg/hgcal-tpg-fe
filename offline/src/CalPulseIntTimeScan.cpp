#include <iostream>
#include <iomanip>
#include <cassert>

#include <yaml-cpp/yaml.h>

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

#include "TFileHandler.h"

#include "RelayReader.h"

void econtEnergies(const uint64_t *p, std::vector<uint16_t> &v) {
  v.resize(24);

  uint64_t e0,e1;
  
  e0=p[0]<<32|p[1]&0xffffffff;
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
  
  e0= p[0]&0xffffffff00000000     |(p[1]&0xffffffff00000000)>>32;
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

//void packetDump(const RecordRunning *rEvent) {
//}


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

  bool beamScintTiming(true);
  unsigned linkNumber(0);

  
  std::ostringstream sout;
  sout << "CalPulseIntTimeScan_Relay" << relayNumber << "_Link" << linkNumber;
  TFileHandler tfh(sout.str().c_str());

  TH1D *hPayloadLength,*hSequence,*hSubpackets,*hSubpacketCount,*hEvents,*hScintLength,*hScintLengthS,*hScintLengthE,*hNtrains;
  TH2D *hPayloadLengthVsSeq,*hSubpacketCountVsPl,*hSequenceVsNe;
  TH2D *hStcVsCpd[2][12],*hToaVsCpd,*hUStcVsCpd[2][12],*hUStcVsCpdM2[2][12],*hScint,*hTpVsCpd[2],*hScintFirst;
  TProfile *pUStcVsCpd[2];
  hEvents=new TH1D("Events",";Scintillator trigger delay;Number of events",
		   200,0,200);

  hSequence=new TH1D("Sequence",";L1A BX number;Number of events",
		     4000,0,4000);

  hPayloadLength=new TH1D("PayloadLength",";Payload length (8-byte words);Number of events",
			  150,0,300);  
  hPayloadLengthVsSeq=new TH2D("PayloadLengthVsSeq",";L1A BX number;Payload length (8-byte words);Number of events",
			  4000,0,4000,150,0,300);

  hStcVsCpd[0][0]=new TH2D("StcVsCpd_E0_STC00",";CalComing-to-L1A delay;STC energy;Number of events",
			   200,0,200,128,0,128);
  hStcVsCpd[1][0]=new TH2D("StcVsCpd_E1_STC00",";CalComing-to-L1A delay;STC energy;Number of events",
			   200,0,200,128,0,128);

  hToaVsCpd=new TH2D("ToaVsCpd",";CalComing-to-L1A delay;TOA value;Number of events",
		     200,0,200,1024,0,1024);
  hTpVsCpd[0]=new TH2D("TpVsCpd_0",";L1A delay;TP 2-bit values;Number of events",
		     200,0,200,4,0,4);
  hTpVsCpd[1]=new TH2D("TpVsCpd_1",";L1A delay;TP 2-bit values;Number of events",
		     200,0,200,4,0,4);

  pUStcVsCpd[0]=new TProfile("PUStcVsCpd_E0_STC00",";CalComing-to-L1A delay;STC energy",
			   200,0,200);
  pUStcVsCpd[1]=new TProfile("PUStcVsCpd_E1_STC00",";CalComing-to-L1A delay;STC energy",
			   200,0,200);
  hUStcVsCpd[0][0]=new TH2D("UStcVsCpd_E0_STC00",";CalComing-to-L1A delay;STC energy;Number of events",
			   200,0,200,128,0,128);
  hUStcVsCpd[1][0]=new TH2D("UStcVsCpd_E1_STC00",";CalComing-to-L1A delay;STC energy;Number of events",
			   200,0,200,128,0,128);
  hUStcVsCpdM2[0][0]=new TH2D("UStcVsCpdM2_E0_STC00",";CalComing-to-L1A delay - 2;STC energy;Number of events",
			      200,0,200,128,0,128);
  hUStcVsCpdM2[1][0]=new TH2D("UStcVsCpdM2_E1_STC00",";CalComing-to-L1A delay - 2;STC energy;Number of events",
			      200,0,200,128,0,128);


  hScint=new TH2D("Scint",";Scintillator trigger delay;Scintillator bits;Number of events",
			      200,0,200,32,0,32);
  hScintFirst=new TH2D("ScintFirst",";Scintillator trigger delay;Scintillator rising edge;Number of events",
			      200,0,200,32*9,0,32*9);

  hNtrains=new TH1D("Ntrains",";Scintillator trigger number of trains;Number of events",
		    100,0,100);
  hScintLength=new TH1D("ScintLength",";Scintillator trigger train length (units of 0.8ps);Number of events",
			9*32,0,9*32);
  hScintLengthS=new TH1D("ScintLengthS",";Starting scintillator trigger train length (units of 0.8ps);Number of events",
			9*32,0,9*32);
  hScintLengthE=new TH1D("ScintLengthE",";Ending scintillator trigger train length (units of 0.8ps);Number of events",
			9*32,0,9*32);

 
  hSubpacketCount=new TH1D("SubpacketCount",";Number of half-HGCROC subpackets;Number of events",
		       20,0,20);
  hSubpacketCountVsPl=new TH2D("SubpacketCountVsPl",";Payload length (8-byte words);Number of half-HGCROC subpackets;Number of events",
			       150,0,300,20,0,20);
  
  hSubpackets=new TH1D("Subpackets",";Half-HGCROC empty subpacket;Number of events",
		       20,0,20);  

  
  // Create the file reader
  Hgcal10gLinkReceiver::RelayReader _relayReader;
  _relayReader.enableLink(0,false); // TEMP
  _relayReader.enableLink(1,false); // TEMP
  _relayReader.enableLink(linkNumber,true); // TEMP
  
  // Make the buffer space for the records
  Hgcal10gLinkReceiver::RecordT<4*4095> *rxxx(new Hgcal10gLinkReceiver::RecordT<4*4095>);
  Hgcal10gLinkReceiver::RecordYaml *r((Hgcal10gLinkReceiver::RecordYaml*)rxxx);

  // Set up specific records to interpet the formats
  const Hgcal10gLinkReceiver::RecordStarting *rStart((Hgcal10gLinkReceiver::RecordStarting*)r);
  const Hgcal10gLinkReceiver::RecordStopping *rStop ((Hgcal10gLinkReceiver::RecordStopping*)r);
  const Hgcal10gLinkReceiver::RecordRunning  *rEvent((Hgcal10gLinkReceiver::RecordRunning*) r);

  // Defaults to the files being in directory "dat"
  // Can call setDirectory("blah") to change this
  //_fileReader.setDirectory("somewhere/else");
  //_relayReader.setDirectory(std::string("dat/Relay")+argv[1]);
  _relayReader.openRelay(relayNumber);


  uint32_t allRecordCount(0);
  uint32_t recordCount[Hgcal10gLinkReceiver::FsmState::EndOfStateEnum+1];
  std::memset(recordCount,0,4*(Hgcal10gLinkReceiver::FsmState::EndOfStateEnum+1));
  
  bool anyError(false);
  unsigned nEvents(0),nErrors[10]={0,0,0,0,0,0,0,0,0,0};
  unsigned bc(0),oc(0),ec(0),seq(0);
  unsigned ocOffset(0);

  bool noCheck(true);
  uint32_t initialSeq;
  
  bool doubleEvents(false);
  unsigned nEventsPerOrbit(1);

  unsigned nEventsPerRun(0);
  unsigned nRuns(0);
  
  std::vector<uint16_t> vEcontEnergy;
  std::vector<uint16_t> vUnpackerEnergy;
  
  uint32_t cpd(0);
  unsigned runCount(0);
  
  while(_relayReader.read(r)) {
    allRecordCount++;
    recordCount[r->state()]++;
    /*
    if(       r->state()==Hgcal10gLinkReceiver::FsmState::Starting) {
      rStart->print();
      std::cout << std::endl;

    } else if(r->state()==Hgcal10gLinkReceiver::FsmState::Stopping) {
      rStop->print();
      std::cout << std::endl;
    */
    if(r->state()!=Hgcal10gLinkReceiver::FsmState::Running) {
      std::cout << std::endl;

      if(r->state()==Hgcal10gLinkReceiver::FsmState::Starting) {
	nRuns++;
	nEventsPerRun=0;
      }
      
      if(r->state()==Hgcal10gLinkReceiver::FsmState::Configuration) {
	YAML::Node n(YAML::Load(r->string()));
	//std::cout << n << std::endl;

	if(beamScintTiming) {
	  if(n["Source"].as<std::string>()=="Serenity") {
	    std::cout << n << std::endl;
	    cpd=n["TriggerDaq"]["daq_ro.DAQro2.latency"].as<uint32_t>()/n["TriggerDaq"]["daq_ro.DAQro2.event_size"].as<uint32_t>();
	    std::cout << "DAQro2 latency = " << n["TriggerDaq"]["daq_ro.DAQro2.latency"].as<uint32_t>() << ", size = " << n["TriggerDaq"]["daq_ro.DAQro2.event_size"].as<uint32_t>() << std::endl;

	  }
	} else {
	  if(n["Source"].as<std::string>()=="TCDS2") {
	    //cpd=n["Configuration"]["ctrl_stat.ctrl.calpulse_delay"].as<uint32_t>();
	    cpd=n["Configuration"]["reg_320.ctrl1.ext_trigger_delay"].as<uint32_t>();
	  }
	}
	
	std::cout << "CalComing-to-L1A delay set to " << cpd << std::endl;
	runCount++;
      }
      
    } else {

	// We have an event record
      nEvents++;
      nEventsPerRun++;
      
      bool print(nEventsPerRun==5);
      anyError=false;

      hEvents->Fill(cpd);
      /*
      if(print) {
	rEvent->print();
	std::cout << std::endl;
      }
      */

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

      const uint64_t *p64(((const uint64_t*)rEvent)+1);

      if(linkNumber==0) {
	//std::cout << "P64[20] = " << std::hex << p64[20] << std::endl;
	/*
	uint64_t e0,e1;

	e0=p64[20]<<32|p64[21]&0xffffffff;
	e1=p64[21]<<48|(p64[22]&0xffffffff)<<16|(p64[23]&0xffffffff)>>16;
	//std::cout << "e0,e1 = " << std::hex << e0 << "," << e1 << std::endl;

	hStcVsCpd[0][ 0]->Fill(cpd,(e0>> 1)&0x7f);
	hStcVsCpd[0][ 0]->Fill(cpd,(e0>> 8)&0x7f);
	hStcVsCpd[0][ 0]->Fill(cpd,(e0>>15)&0x7f);
	hStcVsCpd[0][ 0]->Fill(cpd,(e0>>22)&0x7f);
	hStcVsCpd[0][ 0]->Fill(cpd,(e0>>29)&0x7f);
	hStcVsCpd[0][ 0]->Fill(cpd,(e1    )&0x7f);
	hStcVsCpd[0][ 0]->Fill(cpd,(e1>> 7)&0x7f);
	hStcVsCpd[0][ 0]->Fill(cpd,(e1>>14)&0x7f);
	hStcVsCpd[0][ 0]->Fill(cpd,(e1>>21)&0x7f);
	hStcVsCpd[0][ 0]->Fill(cpd,(e1>>28)&0x7f);
	hStcVsCpd[0][ 0]->Fill(cpd,(e1>>35)&0x7f);
	hStcVsCpd[0][ 0]->Fill(cpd,(e1>>42)&0x7f);

	e0= p64[20]&0xffffffff00000000     |(p64[21]&0xffffffff00000000)>>32;
	e1=(p64[21]&0xffffffff00000000)<<16|(p64[22]&0xffffffff00000000)>>16|(p64[23]&0xffffffff00000000)>>48;
	//std::cout << "e0,e1 = " << std::hex << e0 << "," << e1 << std::endl;

	hStcVsCpd[1][ 0]->Fill(cpd,(e0>> 1)&0x7f);
	hStcVsCpd[1][ 0]->Fill(cpd,(e0>> 8)&0x7f);
	hStcVsCpd[1][ 0]->Fill(cpd,(e0>>15)&0x7f);
	hStcVsCpd[1][ 0]->Fill(cpd,(e0>>22)&0x7f);
	hStcVsCpd[1][ 0]->Fill(cpd,(e0>>29)&0x7f);
	hStcVsCpd[1][ 0]->Fill(cpd,(e1    )&0x7f);
	hStcVsCpd[1][ 0]->Fill(cpd,(e1>> 7)&0x7f);
	hStcVsCpd[1][ 0]->Fill(cpd,(e1>>14)&0x7f);
	hStcVsCpd[1][ 0]->Fill(cpd,(e1>>21)&0x7f);
	hStcVsCpd[1][ 0]->Fill(cpd,(e1>>28)&0x7f);
	hStcVsCpd[1][ 0]->Fill(cpd,(e1>>35)&0x7f);
	hStcVsCpd[1][ 0]->Fill(cpd,(e1>>42)&0x7f);


	for(unsigned i(0);i<6;i++) {
	  // Should be 2*i, 2*i+1
	  hUStcVsCpd[0][0]->Fill(cpd,(p64[69+i]    )&0x7f);
	  hUStcVsCpd[0][0]->Fill(cpd,(p64[69+i]>>13)&0x7f);
	  hUStcVsCpd[1][0]->Fill(cpd,(p64[69+i]>>32)&0x7f);
	  hUStcVsCpd[1][0]->Fill(cpd,(p64[69+i]>>45)&0x7f);

	  hUStcVsCpdM2[0][0]->Fill(cpd,(p64[57+i]    )&0x7f);
	  hUStcVsCpdM2[0][0]->Fill(cpd,(p64[57+i]>>13)&0x7f);
	  hUStcVsCpdM2[1][0]->Fill(cpd,(p64[57+i]>>32)&0x7f);
	  hUStcVsCpdM2[1][0]->Fill(cpd,(p64[57+i]>>45)&0x7f);
	}
	*/
	econtEnergies(p64+20,vEcontEnergy);

	if(nEvents<10) {
	  for(unsigned p(0);p<9;p++) {
	    std::cout << "Unpacker trigger packet " << p << std::endl;
	    unpackerEnergies(p64+69-6*4+6*p,vUnpackerEnergy);

	    for(unsigned i(0);i<24;i++) {
	      std::cout << "Compare STC and USTC " << std::setw(2) << " = "
			<< vEcontEnergy[i] << " ?= " << vUnpackerEnergy[i] << std::endl;
	    }
	  }
	}
	
	unpackerEnergies(p64+69,vUnpackerEnergy);

	for(unsigned i(0);i<24;i++) {
	  hStcVsCpd[ i/12][0]->Fill(cpd,vEcontEnergy[i]);
	  hUStcVsCpd[i/12][0]->Fill(cpd,vUnpackerEnergy[i]);
	  pUStcVsCpd[i/12]->Fill(cpd,vUnpackerEnergy[i]);
	}


	/////////////////////////////////////////////////////////
	// SCINTILLATOR
	/////////////////////////////////////////////////////////

	if(nRuns>1) {

	  unsigned nCentre(116);
	  if(cpd>24) nCentre=118;
	  uint32_t sWord(p64[nCentre]);

	//bool nonZero(false);
	/*
	if(sWord!=0) {
	  std::cout << "Delay = " << std::setw(3) << cpd
		    << ", Scint word = 0x" << std::hex << std::setfill('0')
		    << std::setw(8) << sWord << std::dec << std::setfill(' ') << std::endl;
	  for(int i(-4);i<5;i++) {
	    std::cout  << " Packet scint word = 0x" << std::hex << std::setfill('0')
		       << std::setw(8) << p64[117+5*i] << std::dec << std::setfill(' ') << std::endl;
	  }

	}
	*/
	for(unsigned i(0);i<32;i++) {
	  if((sWord&(1<<i))!=0) hScint->Fill(cpd,i);
	}

	bool inTrain(false);
	unsigned first,last;
	unsigned nTrains(0);
	
	for(unsigned i(0);i<9;i++) {
	  //std::cout << "Scint Word " << i << " = " << std::hex << std::setw(16) << p64[117-20+5*i] << std::dec << std::endl;
	  uint64_t pTemp=p64[nCentre-20+5*i]&0xffffffff;
	  //uint64_t pTemp=p64[nCentre-20+5*i];
	  
	  for(unsigned j(0);j<32;j++) {
	    //std::cout << "Scint Word in loop" << i << " = " << std::hex << std::setw(16) << pTemp << std::dec << std::endl;
	    
	    bool b((pTemp&uint64_t(1<<(31-j)))!=0);

	    //std::cout << "Bit " << 31-j << " = " << std::setw(2) << (b?"True":"False") << ", mask "
	    //	      << std::hex << std::setw(16) << (1<<(31-j)) << std::dec << std::endl;
	    
	    if(b && !inTrain) {
	      first=32*i+j;
	      inTrain=true;
	      nTrains++;
	    }

	    if(b) last=32*i+j;

	    if(!b && inTrain) {
	      inTrain=false;

	      //if(last==first) std::cout << "F=L word = " << std::hex << std::setw(16) << p64[nCentre-20+5*i] << std::dec << std::endl;

	      if(first==0) hScintLengthS->Fill(last-first+1);
	      else {
		hScintFirst->Fill(cpd,first);
		hScintLength->Fill(last-first+1);
	      }
	    }

	    //std::cout << "InTrain = " << (inTrain?"True":"False") << std::endl;
	  }
	}

	if(inTrain) hScintLengthE->Fill(last-first+1);
	//if(nTrains>0) std::cout << "Number of trains = " << nTrains << std::endl;
	hNtrains->Fill(nTrains);
	/*
	if(nTrains>99) {
	  std::cout << "Number of trains = " << std::setw(3) << nTrains
		    << ", Scint word = 0x" << std::hex << std::setfill('0')
		    << std::setw(8) << sWord << std::dec << std::setfill(' ') << std::endl;
	  for(int i(-4);i<5;i++) {
	    std::cout  << " Packet scint word = 0x" << std::hex << std::setfill('0')
		       << std::setw(8) << p64[nCentre+5*i] << std::dec << std::setfill(' ') << std::endl;
	  }
	  rEvent->print();
	}
	*/
	}
	
      } else {
	if(rEvent->payloadLength()>13) {
	  for(unsigned i(10);i<rEvent->payloadLength()-3;i++) {
	    hToaVsCpd->Fill(cpd,(p64[i]    )&0x3ff);
	    hToaVsCpd->Fill(cpd,(p64[i]>>32)&0x3ff);
	    hTpVsCpd[0]->Fill(cpd,(p64[i]>>30)&0x3);
	    hTpVsCpd[1]->Fill(cpd,(p64[i]>>62)&0x3);
	  }
	}
      }
	
      if(!noCheck) {
      if(nEvents==1) {
	seq=rEvent->sequenceCounter();
	initialSeq=seq;
      }
      
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

      if(e->bxId()==0 ||e->bxId()==3564) {
	e->print();
	b->print();
	bph->print();
	print=true;

	  for(unsigned i(0);i<rEvent->payloadLength();i++) {
	    std::cout << "Word " << std::setw(3) << i << " ";
	    std::cout << std::hex << std::setfill('0');
	    std::cout << "0x" << std::setw(16) << *((const uint64_t*)(b+1)+i) << std::endl;
	    std::cout << std::dec << std::setfill(' ');
	  }
	  std::cout << std::endl;


      }
      
      hSequence->Fill(e->bxId());
      hPayloadLength->Fill(rEvent->payloadLength());      
      hPayloadLengthVsSeq->Fill(e->bxId(),rEvent->payloadLength());
      hSequenceVsNe->Fill(nEvents,e->bxId());
      /*      
	if(print && nEvents<1000) {
	  b->print();
	  e->print();
	  
	}
      */
      

      // Access ECON-D packet as an array of 32-bit words
      const uint32_t *pEcond(rEvent->econdPayload());
      
      // Check this is not an empty event
      if(pEcond!=nullptr) {

	bph=(const Hgcal10gLinkReceiver::BePacketHeader*)(p64+2);

	const uint32_t *p32(rEvent->daqPayload());
      /*      
	
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
      }
 }
      if(print) {// && nEvents<10000) {
	  b->print();
	  e->print();
	  
	  std::cout << "Header and words of ECON-D packet " << nEvents << std::endl;
	  rEvent->RecordHeader::print();

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
    }

  }

  std::cout << "Count of all records = " << allRecordCount << std::endl;
  unsigned nCheck(0);
  for(unsigned i(0);i<=Hgcal10gLinkReceiver::FsmState::EndOfStateEnum;i++) {
    std::cout << "Count of record state " << std::setw(3) << i << " = "
	      << Hgcal10gLinkReceiver::FsmState::stateName(Hgcal10gLinkReceiver::FsmState::State(i))
	      << " = " << std::setw(10) << recordCount[i] << std::endl;
    nCheck+=recordCount[i];
  }
  std::cout << "Count of all record states = " << nCheck << std::endl << std::endl;
  
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
