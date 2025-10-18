/**********************************************************************
 Created on : 28/09/2025
 Purpose    : Read the TPG data blocks and decode the energies
 Author     : Indranil Das, Research Associate, Imperial
 Email      : indranil.das@cern.ch | indra.ehep@gmail.com
**********************************************************************/
#include <iostream>

#include "TFileHandlerLocal.h"
#include "FileReader.h"

#include "TPGFEDataformat.hh"
#include "TPGFEModuleEmulation.hh"
#include "Stage1IO.hh"

#include "SlinkBoe.h"
#include "SlinkEoe.h"
#include "TpgSubpacketHeader.h"
#include "OrbitReader.h"

// #include "OrbitDumpEvent.h"
// typedef OrbitDumpEvent OrbitCheckTypedef;

using namespace std;

///Source: offline/inc/OrbitCheck.hxx

class TPGEventReaderTB25 {
public:
  TPGEventReaderTB25(): iEvent(0) {    
  }

  bool runStart(uint32_t run, uint32_t sid) {
    return true;
  }

  bool orbit(const Hgcal10gLinkReceiver::OrbitHeader &oh) {
    oh.print();
    std::cout << std::endl;

    return true;
  }

  bool event(const Hgcal10gLinkReceiver::OrbitReaderEvent &ore) {

    ore._ft->print();
    std::cout << std::endl;
    
    unsigned n64(2*ore._ft->fragmentSize());
    uint64_t *p(ore._array);

    Hgcal10gLinkReceiver::SlinkBoe *b(nullptr);
    if(n64>=2) {
      b = ((Hgcal10gLinkReceiver::SlinkBoe*)p);
      b->print();
      std::cout << std::endl;
    }
    Hgcal10gLinkReceiver::SlinkEoe *e(nullptr);
    if(n64>=4) {
      e = ((Hgcal10gLinkReceiver::SlinkEoe*)(p+n64-2));
      e->print();
      std::cout << std::endl;
    }
    
    
    const Hgcal10gLinkReceiver::TpgSubpacketHeader *tsh(reinterpret_cast<const Hgcal10gLinkReceiver::TpgSubpacketHeader*>(b+1));
    const Hgcal10gLinkReceiver::TpgSubpacketHeader *tshEnd(reinterpret_cast<const Hgcal10gLinkReceiver::TpgSubpacketHeader*>(e-1));

    if(!tsh->validPattern()) {
      tsh->print();
      return false;
    }
    tsh=tsh->nextSubpacketHeader();

    int noffecafe = 0;
    bool done(false);
    while(tsh<=tshEnd && !done) {
      if(!tsh->validPattern()) {
	done=true;
	
      } else {
	if((tsh->channelId()%2)==0) { // RX only
	  std::cout << "Processing event : " <<  iEvent << std::endl;
	  tsh->print();
	  
	  unsigned emp_chan(tsh->channelId()/2);
	  //if(emp_chan==100 or emp_chan==102 or emp_chan==104 or emp_chan==106 or emp_chan==108 or emp_chan==110 or emp_chan==112 or emp_chan==114 or emp_chan==116 or emp_chan==118 or emp_chan==122){
	  //if(emp_chan==100 or emp_chan==102 or emp_chan==104 or emp_chan==106 or emp_chan==108 or emp_chan==110 or emp_chan==112 or emp_chan==114 or emp_chan==116 or emp_chan==118 or emp_chan==120 or emp_chan==122){
	  //if(emp_chan==100 or emp_chan==102 or emp_chan==104 or emp_chan==106 or emp_chan==123){
	  if(emp_chan==100){
	    uint wpspd = 0;
	    for(unsigned bx(0);bx<tsh->numberOfBxs();bx++) {
	      const uint64_t *el64packed((const uint64_t*)(tsh+1+bx*tsh->numberOfWordsPerBx()));
	      uint32_t elinks[8];
	      for(unsigned j(0);j<tsh->numberOfWordsPerBx();j++) {
		// if(j<(tsh->numberOfWordsPerBx()-1)){
		//   elinks[2*j] = el64packed[j] & 0xffffffff;
		//   elinks[2*j+1] = (el64packed[j]>>32) & 0xffffffff;
		// }else{
		//   elinks[2*j] = el64packed[j] & 0xffffffff;
		// }
		elinks[2*j] = el64packed[j] & 0xffffffff;
		elinks[2*j+1] = (el64packed[j]>>32) & 0xffffffff;
		std::cout << "Word " << std::setw(6) << wpspd++ << " = 0x"
			  << std::hex << std::setfill('0')
			  << std::setw(16) << el64packed[j]
			  << std::dec << std::setfill(' ')
			  << std::endl;	      
	      }
	      for(unsigned iel(0);iel<8;iel++) {
		std::cout << "\t elink " << std::setw(3) << iel << " = 0x"
			  << std::hex << std::setfill('0')
			  << std::setw(8) << elinks[iel]
			  << std::dec << std::setfill(' ')
			  << std::endl;	      
	      }

	      if(emp_chan!=123){
		// /////////////////////////// Si ////////////////////////////
		const int neTx = 4;
		uint32_t el[neTx];
		el[0] = elinks[2];
		el[1] = elinks[1];
		el[2] = elinks[0];
		if(neTx>3) el[3] = elinks[3];

		// el[0] = elinks[0];
		// el[1] = elinks[1];
		// el[2] = elinks[2];
		// if(neTx>3) el[3] = elinks[3];
	      
		// //Run 111137 and 111138
		// el[0] = elinks[6];
		// el[1] = elinks[5];
		// el[2] = elinks[4];
		// if(neTx>3) el[3] = elinks[3];
	      
		for(unsigned iel(0);iel<neTx;iel++){
		  std::cout << "\t\t el " << std::setw(3) << iel << " = 0x"
			    << std::hex << std::setfill('0')
			    << std::setw(8) << el[iel]
			    << std::dec << std::setfill(' ')
			    << std::endl;	      		
		}
	      
		TPGFEDataformat::TcRawDataPacket rdp;
		if(neTx==4){
		  TPGStage1Emulation::Stage1IO::convertElinksToTcRawData(TPGFEDataformat::BestC, 9, el, rdp);
		  //TPGStage1Emulation::Stage1IO::convertElinksToTcRawData(TPGFEDataformat::STC4A, 12, el, rdp);
		}else
		  TPGStage1Emulation::Stage1IO::convertElinksToTcRawData(TPGFEDataformat::BestC, 6, el, rdp);
		rdp.print();
		// /////////////////////////// Si ////////////////////////////
	      }else{
		/////////////////////////// Sci ////////////////////////////
		const int neTx1 = 4;
		const int neTx2 = 3;
		uint32_t *el1 = new uint32_t[neTx1];
		uint32_t *el2 = new uint32_t[neTx2];
		el1[0] = elinks[1];
		el1[1] = elinks[0];
		el1[2] = elinks[2];
		el1[3] = elinks[3];
	      
		el2[0] = elinks[4];
		el2[1] = elinks[5];
		el2[2] = elinks[6];
	      
		//Run 111137 and 111138
		// el[0] = elinks[6];
		// el[1] = elinks[5];
		// el[2] = elinks[4];
		// if(neTx>3) el[3] = elinks[3];
	      
		for(unsigned iel(0);iel<neTx1;iel++){
		  std::cout << "\t\t el1 " << std::setw(3) << iel << " = 0x"
			    << std::hex << std::setfill('0')
			    << std::setw(8) << el1[iel]
			    << std::dec << std::setfill(' ')
			    << std::endl;	      		
		}
	      
		for(unsigned iel(0);iel<neTx2;iel++){
		  std::cout << "\t\t el1 " << std::setw(3) << iel << " = 0x"
			    << std::hex << std::setfill('0')
			    << std::setw(8) << el2[iel]
			    << std::dec << std::setfill(' ')
			    << std::endl;	      		
		}
	      
		TPGFEDataformat::TcRawDataPacket rdp1, rdp2;
		TPGStage1Emulation::Stage1IO::convertElinksToTcRawData(TPGFEDataformat::STC4A, 12, el1, rdp1);
		rdp1.print();
		TPGStage1Emulation::Stage1IO::convertElinksToTcRawData(TPGFEDataformat::STC4A, 10, el2, rdp2);
		rdp2.print();

		uint64_t tot0 = 0, tot1 = 0;
		for(const auto& itc: rdp1.getTcData()) tot0 += itc.decodedE(rdp1.type());
		for(const auto& itc: rdp2.getTcData()) tot1 += itc.decodedE(rdp2.type());
		std::cout << "tot0: " << tot0 << ", tot1: " << tot1 << std::endl;
	      
		delete []el1;
		delete []el2;
		/////////////////////////// Sci ////////////////////////////
	      }
	    }
	  }
	}
	
	tsh=tsh->nextSubpacketHeader();
      }
      noffecafe++;
    }

    for(unsigned j(0);j<n64;j++) {
      std::cout << "Word " << std::setw(6) << j << " = 0x"
		<< std::hex << std::setfill('0')
		<< std::setw(16) << p[j]
		<< std::dec << std::setfill(' ')
	      << std::endl;
    }
    std::cout << std::endl;
 
    std::cout << "Processing event : " <<  iEvent++ << ", noffecafe: " << noffecafe << std::endl;
    return true;
  }

  bool runStop() {
    return true;
  }
private:
  uint64_t iEvent;
};

typedef TPGEventReaderTB25 OrbitCheckTypedef;

int main(int argc, char** argv){
  
  std::cout << "Nof arguments : " << argc << std::endl;
  if(argc < 3){
    std::cerr << argv[0] << ": no run numbers specified" << std::endl;
    return 1;
  }
  
  //Command line arg assignment
  //Assign relay and run numbers
  uint32_t runNumber(12600113);
  uint32_t sourceId(1260);
  uint32_t firstLs(1);
  
  unsigned dumpEvent(0);
  
  std::istringstream issRun(argv[1]);
  issRun >> runNumber;
  std::istringstream issLink(argv[2]);
  issLink >> sourceId;
  if(argc > 3){
    std::istringstream isfirstLs(argv[3]);
    isfirstLs >> firstLs;
  }
  
  /////========================================================
  std::vector<Hgcal10gLinkReceiver::OrbitReaderEvent> vEvents;
  
  OrbitCheckTypedef ct;
  assert(ct.runStart(runNumber,sourceId));
  
  std::string oDir("dat/");
  
  Hgcal10gLinkReceiver::OrbitReader oReader;
  oReader.setPrint(false);
  
  unsigned nEvents(0);
  
  const Hgcal10gLinkReceiver::OrbitHeader *oh;
  //Hgcal10gLinkReceiver::FragmentTrailer *ft;
  
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
    
    while((oh=oReader.readOrbit(vEvents))!=nullptr and nEvents<4) {
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
  /////========================================================
  
  return 0;
}

