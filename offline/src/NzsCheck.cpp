#include <iostream>
#include <iomanip>
#include <cassert>

#include <yaml-cpp/yaml.h>

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

#include "TFileHandler.h"
#include "FileReader.h"
#include "EcontEnergies.h"

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

  std::ostringstream sout;
  sout << "NzsCheck_Relay" << relayNumber;
  TFileHandler tfh(sout.str().c_str());
  
  TH1D *hPacketSize;
  TH2D *hPacketSizeVsEvent,*hPacketSizeVsBx,*hPacketSizeVsDbx;
  
  hPacketSize=new TH1D("PacketSize",";Packet size (64-bit words);Number of events",
		       60,0,120);

  hPacketSizeVsEvent=new TH2D("PacketSizeVsEvent",";Event counter;Packet size (64-bit words);Number of events",
			 600,0,600000,60,0,120);
  hPacketSizeVsBx=new TH2D("PacketSizeVsBx",";Bunch counter;Packet size (64-bit words);Number of events",
			 3600,0,3600,60,0,120);
  hPacketSizeVsDbx=new TH2D("PacketSizeVsDbx",";Difference of bunch counter;Packet size (64-bit words);Number of events",
			 1000,0,1000,60,0,120);
  /*
  TH1D *hScintFirstWord,*hScintNwords,*hStcEnergy[12],*hStcEnergy0[12];
  TH2D *hScintWord;

  TH1D *hPayloadLength,*hSequence,*hSubpackets,*hSubpacketCount,*hEvents,*hScintLength,*hScintLengthS,*hScintLengthE,*hNtrains,*hScintFirst,*hScintFirst1,*hScintFirst2;
  TH2D *hPayloadLengthVsSeq,*hSubpacketCountVsPl,*hSequenceVsNe,*hEcontVsBc,*hScintFirstStc[12];
  TH2D *hStcVsCpd[2][12],*hToaVsCpd,*hUStcVsCpd[2][12],*hUStcVsCpdM2[2][12],*hScint,*hTpVsCpd[2],*hScintLengthVsFirst;
  TProfile *pUStcVsCpd[2];





  
  hEcontVsBc=new TH2D("EcontVsBc",";BC value;ECON-T counter value;Number of events",
		      3600,0,3600,16,0,16);

  
  hScintNwords=new TH1D("ScintNwords",";Scintillator number of words;Number of events",
			256,0,256);
  hScintWord=new TH2D("ScintWord",";Scintillator word type;Number of events",
		      256,0,256,4,0,4);


  hScint=new TH2D("Scint",";Scintillator trigger delay;Scintillator bits;Number of events",
                              200,0,200,32,0,32);
  hScintFirst=new TH1D("ScintFirst",";Scintillator rising edge;Number of events",
                              32*9,0,32*9);
  hScintFirst1=new TH1D("ScintFirst1",";Scintillator rising edge;Number of events",
			32*9,0,32*9);
  hScintFirst2=new TH1D("ScintFirst2",";Scintillator rising edge;Number of events",
			32*9,0,32*9);

  for(unsigned i(0);i<12;i++) {
    std::ostringstream sout;
    sout << "StcEnergy" << std::setfill('0') << std::setw(2) << i;
    hStcEnergy[i]=new TH1D(sout.str().c_str(),";STC energy;Number of events",
			   128,0,128);

    hStcEnergy0[i]=new TH1D((sout.str()+"0").c_str(),";STC energy;Number of events",
			    128,0,128);

    std::ostringstream sout2;
    sout2 << "ScintFirstStc" << std::setfill('0') << std::setw(2) << i;
    hScintFirstStc[i]=new TH2D(sout2.str().c_str(),";Scintillator rising edge;STC energy;Number of events",
			       32*9,0,32*9,128,0,128);
  }
  
  hNtrains=new TH1D("Ntrains",";Scintillator trigger number of trains;Number of events",
                    100,0,100);
  hScintLength=new TH1D("ScintLength",";Scintillator trigger train length (units of 0.8ps);Number of events",
                        9*32,0,9*32);
  hScintLengthS=new TH1D("ScintLengthS",";Starting scintillator trigger train length (units of 0.8ps);Number of events",
                        9*32,0,9*32);
  hScintLengthE=new TH1D("ScintLengthE",";Ending scintillator trigger train length (units of 0.8ps);Number of events",
                        9*32,0,9*32);

  hScintLengthVsFirst=new TH2D("hScintLengthVsFirst",";Scintillator train rising edge;Scintillator train length;Number of events",
			       32*9,0,32*9,200,0,200);

  std::vector<uint16_t> v;
  */
  
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
  //_fileReader.openRun(runNumber,2);
  
  //bool anyError(false);
  unsigned nEvents(0);//,nErrors[10]={0,0,0,0,0,0,0,0,0,0};
  unsigned bxOld;
  /*
  unsigned bc(0),oc(0),ec(0),seq(0);
  unsigned ocOffset(0);

  bool noCheck(false);
  uint32_t initialSeq;

  
  bool doubleEvents(false);
  unsigned nEventsPerOrbit(1);
  */

  //std::ofstream fout("Aidan.csv");
  
  
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
      bool print(nEvents<=10 || rEvent->payloadLength()<90);
      //anyError=false;

      const uint64_t *p64(((const uint64_t*)rEvent)+1);

      const Hgcal10gLinkReceiver::SlinkBoe *boe(rEvent->slinkBoe());
      assert(boe!=nullptr);
      const Hgcal10gLinkReceiver::SlinkEoe *eoe(rEvent->slinkEoe());
      assert(eoe!=nullptr);

      hPacketSize->Fill(rEvent->payloadLength());
      hPacketSizeVsEvent->Fill(boe->eventId(),rEvent->payloadLength());
      hPacketSizeVsBx->Fill(eoe->bxId(),rEvent->payloadLength());
      if(nEvents>1) hPacketSizeVsDbx->Fill(eoe->totalBx()-bxOld,rEvent->payloadLength());

      bxOld=eoe->totalBx();


      if(print) {
	std::cout << std::endl;
	rEvent->print();

	// Access the Slink header ("begin-of-event")
	// This should always be present; check pattern is correct
	const Hgcal10gLinkReceiver::SlinkBoe *b(rEvent->slinkBoe());
	assert(b!=nullptr);
	//if(b->l1aType()==0x001c) b->print();
      
	// Access the Slink trailer ("end-of-event")
	// This should always be present; check pattern is correct
	const Hgcal10gLinkReceiver::SlinkEoe *e(rEvent->slinkEoe());
	assert(e!=nullptr);
	

	b->print();
	e->print();
      
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

  //fout.close();
  delete r;

  return 0;
}
