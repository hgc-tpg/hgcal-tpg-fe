/**********************************************************************
 Created on : 06/10/2025
 Purpose    : Read the DAQ data blocks
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

using namespace std;

namespace Hgcal10gLinkReceiver {
  
  class ErxPassThSubpacket {
    
  public:
    ErxPassThSubpacket() {
      reset();
    }
    virtual ~ErxPassThSubpacket() {
      reset();
    }
    
    void reset() {
      if(!_data32) delete []_data32;
      _data32 = 0x0;
      ffidxlst.resize(0);
      n32 = 0;
    }

    bool validNeRxs(uint32_t neRx) const {      
      return (ffidxlst.size()==neRx);
    }
    
    bool setData(uint64_t *p, uint32_t gwpos, uint32_t size) {
      //uint32_t n32 = EconDSubpacket::getPayloadSize() + 1;
      //uint32_t neRx = (EconDSubpacket::getPayloadSize() - 1)/39;
      n32 = size;
      _data32 = new uint32_t[n32];
      for(unsigned j(0);j<n32/2;j++) {
	_data32[2*j] = p[gwpos+j]>>32;
	_data32[2*j+1] = p[gwpos+j];
      }
      for(unsigned j(0);j<n32;j++) if(_data32[j]==0xffffffff) ffidxlst.push_back(j);
      uint32_t neRx = ((n32-1) - 1)/39;
      //Validate that all channels are present for all halfrocs connected to ECON-D
      if(!validNeRxs(neRx)) {reset(); return false;} 
      return true;
    }

    uint32_t getNeRx() const {return ffidxlst.size();}
    uint32_t getChData(uint32_t ihroc, uint32_t ich) const {
      assert(ihroc<getNeRx() and ich<37);
      return _data32[(ffidxlst[ihroc] + 1 + ich)];
    }
    
    void print(std::ostream &o=std::cout, std::string s="") {
      if(s.find("v")!=std::string::npos or s.find("V")!=std::string::npos){ //verbose
	for(unsigned j(0);j<n32;j++) {
	  o << s << "ErxPassThSubpacket::print() Word32 "
	    << std::setw(6) << j << " = 0x"
	    << std::hex << std::setfill('0')
	    << std::setw(8) << _data32[j]
	    << std::dec << std::setfill(' ')
	    << std::endl;
	}
      }else{
	o << s << "ErxPassThSubpacket::print() nofeRxs : " << getNeRx() << std::endl;
	for(unsigned j(0);j<getNeRx();j++) {
	  uint32_t idata = ffidxlst.at(j) + 1;
	  for(unsigned k(0);k<37;k++) {
	    o << s << "ErxPassThSubpacket::print() "
	      << "\t iData ith-half-roc: " << std::setfill('0') << std::setw(2) << j
	      <<", ch: " << std::setfill('0') << std::setw(2) << k
	      << " = 0x"
	      << std::hex << std::setfill('0')
	      << std::setw(8) << getChData(j,k)
	      << std::dec << std::setfill(' ')
	      << std::endl;
	  }//channel loop
	  o << s << std::endl;
	}//eRx loop
      }
    }
    
  private:
    uint32_t* _data32;
    uint32_t n32;
    std::vector<int> ffidxlst;
    
  };

  class EconDSubpacket {
    
  public:
    EconDSubpacket() {
      reset();
    }
    
    void reset() {
      _data=0x0;
      gwpos = 2; //first two 64-bit words are slink headers
    }
    
    void setWordPos(uint32_t iw)	{wpos	= iw;}    
    void setEventWordPos(uint32_t iw)	{gwpos	= iw;}
    void resetWordPos()			{setWordPos(0);}    
    void resetEventWordPos()		{setEventWordPos(0);}
    
    void setData(uint64_t *p) {
      _data=p;
      wpos = 0;
    }
    
    uint64_t getWordPos() const		{return wpos;}
    uint64_t getEventWordPos() const	{return gwpos;}
    uint64_t getEconDHeaderWord() const {return uint64_t(*(_data));}
    
    uint16_t numberOfValidEcons() const {
      uint16_t nofvalidecons = 0;
      for(int iecon=0;iecon<12;iecon++) if( ((getEconDHeaderWord()>>(3*iecon))&0x3)==0x0) nofvalidecons++;      
      return nofvalidecons;
    }
    
    uint64_t getWord(int iw) const	{return uint64_t(*(_data+iw));}    
    uint64_t getEventWord() const	{return getWord(wpos);}    
    uint32_t getPayloadSize() const	{return ((getEventWord()>>46)&0x1ff);}
    bool isPassTh() const		{return (((getEventWord()>>45)&0x1)==1)?true:false;}
    
    int nextEcon(int iecon){
      if(iecon==0)
	wpos = wpos + 1;
      else 
	wpos = wpos + ((getPayloadSize() + 1)/2 + 1);
      return wpos;
    }
    
    bool nextSubpacketHeader(uint32_t maxNof64Words) {
      resetWordPos();
      for(uint32_t iecon=0;iecon<=numberOfValidEcons();iecon++) nextEcon(iecon);
      gwpos += wpos;
      if(gwpos>=maxNof64Words) return false;      
      setData(_data+wpos);      
      return true;
    }
    
    void print(std::ostream &o=std::cout, std::string s="") {
      resetWordPos();
      o << s << "EconDSubpacket::print()  EconDHeaderWord = 0x"
        << std::hex << std::setfill('0')
        << std::setw(16) << getEconDHeaderWord()
        << std::dec << std::setfill(' ')
	<< std::endl;
      for(uint32_t iecon=0;iecon<numberOfValidEcons();iecon++) {
	nextEcon(iecon);
	o << s << " EconD event header Word = 0x"
	  << std::hex << std::setfill('0')
	  << std::setw(8) << getEventWord()
	  << std::dec << std::setfill(' ')
	  << ", payload size: "
	  << std::setw(4) << getPayloadSize()       
	  << ", word pos: "
	  << std::setw(4) << getWordPos()
	  << ", global word pos: "
	  << std::setw(4) << getEventWordPos()
	  << std::endl;
      }

      o << s << " Number of Valid Econs = "
	<< std::setw(3) << unsigned(numberOfValidEcons())
	<< std::endl;
      o << s << std::endl;
    }
    
  private:
    uint64_t* _data;
    uint32_t wpos;
    uint32_t gwpos;
  };

}

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
    
    Hgcal10gLinkReceiver::EconDSubpacket *econdsp = new Hgcal10gLinkReceiver::EconDSubpacket();
    econdsp->setData(p+2);
    do{
      econdsp->print();
      econdsp->resetWordPos();
      for(uint32_t iecon=0;iecon<econdsp->numberOfValidEcons();iecon++) {
	econdsp->nextEcon(iecon);
	if(econdsp->isPassTh()){
	  ///////////////// The following details are replaced by class ///////////////////////
	  // uint32_t neRx = (econdsp->getPayloadSize() - 1)/39;
	  // std::cout << "iecon: " << iecon << ", neRx: " << neRx << ", globpos: " << econdsp->getEventWordPos() << ", locpos: " << econdsp->getWordPos() << std::endl;
	  // uint32_t n32 = econdsp->getPayloadSize() + 1;
	  // uint32_t *p32 = new uint32_t[n32];
	  // unsigned fpos = econdsp->getEventWordPos()+econdsp->getWordPos()+1;
	  // for(unsigned j(0);j<n32/2;j++) {
	  //   p32[2*j] = p[fpos+j]>>32;
	  //   p32[2*j+1] = p[fpos+j];
	  // }
	  // std::vector<int> ffidxlst;
	  // for(unsigned j(0);j<n32;j++) if(p32[j]==0xffffffff) ffidxlst.push_back(j);
	  // std::cout << "Nof ffidxlst: " << ffidxlst.size() << std::endl;
	  // assert(ffidxlst.size()==neRx);
	  // for(unsigned j(0);j<ffidxlst.size();j++) {
	  //   uint32_t idata = ffidxlst.at(j) + 1;
	  //   for(unsigned k(0);k<37;k++) {
	  //     std::cout << "\t\tiData j: " << std::setw(6) << j
	  // 		<<", idata: " << idata
	  // 		<<", k: " << k
	  // 		<< " = 0x"
	  // 		<< std::hex << std::setfill('0')
	  // 		<< std::setw(8) << p32[idata+k]
	  // 		<< std::dec << std::setfill(' ')
	  // 		<< std::endl;
	  //   }
	  // }
	  // ffidxlst.clear();
	  // for(unsigned j(0);j<n32;j++) {
	  //   std::cout << "\tWord " << std::setw(6) << j << " = 0x"
	  // 	      << std::hex << std::setfill('0')
	  // 	      << std::setw(8) << p32[j]
	  // 	      << std::dec << std::setfill(' ')
	  // 	      << std::endl;
	  // }
	  // std::cout << std::endl;
	  // delete []p32;
	  ///////////////// The above details are replaced by class ///////////////////////
	  uint32_t n32 = econdsp->getPayloadSize() + 1;
	  unsigned fpos = econdsp->getEventWordPos()+econdsp->getWordPos()+1;
	  Hgcal10gLinkReceiver::ErxPassThSubpacket *desp = new Hgcal10gLinkReceiver::ErxPassThSubpacket();
	  if(!desp->setData(p,fpos,n32)) continue;
	  desp->print(std::cout,"verbose::");
	}	
      }//econ loop;
    }while(econdsp->nextSubpacketHeader(n64-2));
    delete econdsp;

    // for(unsigned j(0);j<n64;j++) {
    //   std::cout << "Word " << std::setw(6) << j << " = 0x"
    // 		<< std::hex << std::setfill('0')
    // 		<< std::setw(16) << p[j]
    // 		<< std::dec << std::setfill(' ')
    // 	      << std::endl;
    // }
    // std::cout << std::endl;
 
    std::cout << "Processing event : " <<  iEvent++
      //<< ", noffecafe: " << noffecafe
	      << std::endl;
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

