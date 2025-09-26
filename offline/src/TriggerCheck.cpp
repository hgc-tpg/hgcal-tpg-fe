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
  sout << "TriggerCheck_Relay" << relayNumber;
  TFileHandler tfh(sout.str().c_str());

  TH1D *hScintFirstWord,*hScintNwords,*hStcEnergy[12],*hStcEnergy0[12];
  TH2D *hScintWord;

  TH1D *hPayloadLength,*hSequence,*hSubpackets,*hSubpacketCount,*hEvents,*hScintLength,*hScintLengthS,*hScintLengthE,*hNtrains,*hScintFirst,*hScintFirst1,*hScintFirst2;
  TH2D *hPayloadLengthVsSeq,*hSubpacketCountVsPl,*hSequenceVsNe,*hEcontVsBc,*hScintFirstStc[12];
  TH2D *hStcVsCpd[2][12],*hToaVsCpd,*hUStcVsCpd[2][12],*hUStcVsCpdM2[2][12],*hScint,*hTpVsCpd[2],*hScintLengthVsFirst;
  TProfile *pUStcVsCpd[2];


  hEcontVsBc=new TH2D("EcontVsBc",";BC value;ECON-T counter value;Number of events",
		      3600,0,3600,16,0,16);

  
  hScintFirstWord=new TH1D("ScintFirstWord",";Scintillator first word location;Number of events",
			   256,0,256);
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
      bool print(nEvents<=1000);// || nEvents>2000);
      //anyError=false;

      const uint64_t *p64(((const uint64_t*)rEvent)+1);

      /*
      // Run 1691537835 only
      uint64_t RaghuWord(p64[ 96]);
      uint64_t AidanWord(p64[145]);
      */
      
      const Hgcal10gLinkReceiver::SlinkBoe *boe(rEvent->slinkBoe());
      assert(boe!=nullptr);
      const Hgcal10gLinkReceiver::SlinkEoe *eoe(rEvent->slinkEoe());
      assert(eoe!=nullptr);

      /*
      fout << std::setw(4) << eoe->bxId() << "," << std::setw(10) << eoe->orbitId()
		<< "," << std::setw(10) << boe->eventId() << ",0x"
		<< std::hex << std::setfill('0')
		<< std::setw(8) << (AidanWord>>32) << ",0x"
		<< std::setw(8) << (RaghuWord&0xffffffff)
	   << std::dec << std::setfill(' ') << std::endl;
      */
      
      
      assert((p64[3]&0xffffffffffffff00)==0xfecafecafecafe00);
      unsigned nEcontWords(p64[3]&0xff);
      assert((nEcontWords%4)==0);
      unsigned nEcontBxs(nEcontWords/4);
      unsigned nEcontStart(3+1);
      unsigned nEcontCentre(nEcontStart+4*(nEcontBxs/2));
      //for(unsigned i(nEcontStart);i<nEcondStart+nEcontWords;i+=4) {
      //}
      
      assert((p64[nEcontStart+nEcontWords]&0xffffffffffffff00)==0xfecafecafecafe00);
      unsigned nUnpackerWords(p64[nEcontStart+nEcontWords]&0xff);
      assert((nUnpackerWords%6)==0);
      unsigned nUnpackerBxs(nUnpackerWords/6);
      unsigned nUnpackerStart(nEcontStart+nEcontWords+1);
      unsigned nUnpackerCentre(nUnpackerStart+6*(nUnpackerBxs/2));
      //for(unsigned i(nUnpackerStart);i<nUnpackerStart+nUnpackerWords;i+=6) {
      //}
      
      unsigned unpackerBx(Hgcal10gLinkReceiver::unpackerEnergies(p64+nUnpackerCentre,v));
      hEcontVsBc->Fill(eoe->bxId(),unpackerBx);
      if((boe->l1aSubType()&0x1)!=0) {
	for(unsigned stc(0);stc<12;stc++) {
	  hStcEnergy[stc]->Fill(v[stc]);
	}
      }

      // Scintillator
      unsigned nWord(0),nFirstWord(0),nLastWord(0);
      for(unsigned i(0);i<rEvent->payloadLength();i++) {
	if((p64[i]>>32)==0xabcdfeed) {
	  if(nWord==0) {
	    nFirstWord=i;
	    hScintFirstWord->Fill(nFirstWord);
	  }
	  
	  uint32_t d(p64[i]&0xffffffff);
	  if(d==0x000000000)  hScintWord->Fill(i,1);
	  else if(d==0xffffffff)  hScintWord->Fill(i,3);
	  else hScintWord->Fill(i,2);

	  nWord++;
	}
      }
      hScintNwords->Fill(nWord);


      unsigned nCentre(nFirstWord-2+nWord/2);
      nCentre=176; // HACK
	
      //std::cout << "nCentre = " << nCentre << std::endl;
      if(nCentre<300) {
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
            //        << std::hex << std::setw(16) << (1<<(31-j)) << std::dec << std::endl;
            
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
                hScintFirst->Fill(first);
                if((boe->l1aSubType()&0x1)!=0) {
		  hScintFirst1->Fill(first);
		  for(unsigned stc(0);stc<12;stc++) {
		    hScintFirstStc[stc]->Fill(first,v[stc]);
		  }
		}
		if((boe->l1aSubType()&0x2)!=0) hScintFirst2->Fill(first);

                hScintLength->Fill(last-first+1);
		hScintLengthVsFirst->Fill(first,last-first+1);
	      }
            }

            //std::cout << "InTrain = " << (inTrain?"True":"False") << std::endl;
          }
        }

        if(inTrain) hScintLengthE->Fill(last-first+1);
        //if(nTrains>0) std::cout << "Number of trains = " << nTrains << std::endl;
        hNtrains->Fill(nTrains);

	if(nTrains==0) {
	  for(unsigned stc(0);stc<12;stc++) {
	    hStcEnergy0[stc]->Fill(v[stc]);
	  }
	}
      }

	
      
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
