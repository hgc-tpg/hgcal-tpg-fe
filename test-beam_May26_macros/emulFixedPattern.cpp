/**********************************************************************
 Created on : 07/10/2025
 Purpose    : Emulation with fixed pattern data for test-beam of September 2025
 Author     : Indranil Das, Research Associate
 Email      : indranil.das@cern.ch | indra.ehep@gmail.com
**********************************************************************/
#include <iostream>

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"

#include "TFileHandlerLocal.h"
#include "FileReader.h"

#include "TPGFEDataformat.hh"
#include "TPGFEConfiguration.hh"
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
    uint16_t getBx() const		{return ((getEventWord()>>20)&0xfff) ;}
    
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
      o << s << " Number of Valid Econs = "
	<< std::setw(3) << unsigned(numberOfValidEcons())
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
	  << ", isPassThrough: "
	  << std::setw(4) << isPassTh()
	  << ", getBx: "
	  << std::setw(4) << getBx()
	  << std::endl;
      }
      o << s << std::endl;
    }
    
  private:
    uint64_t* _data;
    uint32_t wpos;
    uint32_t gwpos;
  };

}

class DAQEventReaderTB25 {
public:
  DAQEventReaderTB25(): iEvent(0) {    
  }

  bool runStart(uint32_t run, uint32_t sid) {
    return true;
  }

  bool orbit(const Hgcal10gLinkReceiver::OrbitHeader &oh) {
    oh.print();
    std::cout << std::endl;

    return true;
  }

  bool event(const Hgcal10gLinkReceiver::OrbitReaderEvent &ore, uint64_t eventId, uint32_t modId, std::map<uint64_t,TPGFEDataformat::HRocarray>& rocarr) {

    //============================
    //rocarr.clear();
    TPGFEConfiguration::TPGFEIdPacking pck;
    pck.setModId(modId);
    
    uint32_t zside = pck.getZside();
    uint32_t sector = pck.getSector();
    uint32_t link = pck.getLink();
    uint32_t det = pck.getDetType();
    uint32_t econd = pck.getEconN();
    uint32_t econt = pck.getEconN();
    uint32_t selTC4 = pck.getSelTC4();
    uint32_t module = pck.getModule();
    uint32_t iroc = 0;
    uint32_t daqlink = 0;
    //============================
    
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
	  uint32_t n32 = econdsp->getPayloadSize() + 1;
	  unsigned fpos = econdsp->getEventWordPos()+econdsp->getWordPos()+1;
	  Hgcal10gLinkReceiver::ErxPassThSubpacket *desp = new Hgcal10gLinkReceiver::ErxPassThSubpacket();
	  if(!desp->setData(p,fpos,n32)) { delete desp ; continue;}
	  //desp->print(std::cout,"verbose::");
	  //========================================================
	  if(iecon==econd and link==daqlink){
	    std::cout << "Econ: " << iecon << ", link: " << link << std::endl;
	    uint32_t rocorder[6] = {4,5,2,3,0,1} ;
	    for(unsigned ihrc(0);ihrc<desp->getNeRx();ihrc++) {
	      TPGFEDataformat::HalfHgcrocData hrocdata;
	      int ich = 0;
	      for(unsigned iseq(0);iseq<37;iseq++) {
		if(iseq==18) continue;
		TPGFEDataformat::HalfHgcrocChannelData chdata;
		//uint32_t rawdata = desp->getChData(ihrc,iseq);
		uint32_t rawdata = desp->getChData(rocorder[ihrc],iseq);
		const uint16_t trigflag = (rawdata>>30) & 0x3;
		if(trigflag==3){
		  uint16_t ttot = (rawdata>>10) & 0x3FF;
		  chdata.setTot(ttot,trigflag);
		}else if(trigflag==2){
		  uint16_t ttot = (rawdata>>10) & 0x3FF;
		  chdata.setTot(ttot,trigflag);
		}else if(trigflag==1){
		  uint16_t tadc = (rawdata>>10) & 0x3FF;
		  chdata.setAdc(0,trigflag);
		}else if(trigflag==0){
		  uint16_t tadc = (rawdata>>10) & 0x3FF;
		  chdata.setAdc(tadc,trigflag);
		}else{
		  chdata.setZero();
		}
		hrocdata.setChannel(ich, chdata);
		ich++;
	      }//channel loop
	      hrocdata.setBx(econdsp->getBx());
	      hrocdata.setSlinkBx(econdsp->getBx());
	      pck.setZero();
	      int ihalf = (ihrc%2==0)?0:1;
	      iroc = ihrc/2;
	      rocarr[eventId].push_back(std::make_pair(pck.packRocId(zside, sector, link, det, econd, selTC4, module, iroc, ihalf),hrocdata));
	    }//eRx loop
	    
	  }
	  //========================================================
	  delete desp;
	}	
      }//econ loop;
      daqlink++;
    }while(econdsp->nextSubpacketHeader(n64-2));
    delete econdsp;
    
    std::cout << "Processing event : " <<  iEvent++
	      << std::endl;
    return true;
  }

  bool runStop() {
    return true;
  }
private:
  uint64_t iEvent;
};

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

  bool event(const Hgcal10gLinkReceiver::OrbitReaderEvent &ore, uint64_t eventId, uint32_t modId, std::map<uint64_t,TPGFEDataformat::TcModuleBxPackets>& tpgarray) {

    //================================
    //tpgarray.clear();
    TPGFEConfiguration::TPGFEIdPacking pck;
    pck.setModId(modId);
    
    uint32_t zside = pck.getZside();
    uint32_t sector = pck.getSector();
    uint32_t link = pck.getLink();
    uint32_t det = pck.getDetType();
    uint32_t econd = pck.getEconN();
    uint32_t econt = pck.getEconN();
    uint32_t selTC4 = pck.getSelTC4();
    uint32_t module = pck.getModule();
    uint32_t iroc = 0;    
    //================================
    
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
	  // // std::cout << "Processing event : " <<  iEvent << std::endl;
	  // // tsh->print();
	  
	  unsigned emp_chan(tsh->channelId()/2);
	  //if(emp_chan==100 or emp_chan==102 or emp_chan==104 or emp_chan==106 or emp_chan==108 or emp_chan==110 or emp_chan==112 or emp_chan==114 or emp_chan==116 or emp_chan==118 or emp_chan==122){
	  //if(emp_chan==100 or emp_chan==102 or emp_chan==104 or emp_chan==106 or emp_chan==108 or emp_chan==110 or emp_chan==112 or emp_chan==114 or emp_chan==116 or emp_chan==118 or emp_chan==120 or emp_chan==122){
	  if(emp_chan==104){ 
	    uint wpspd = 0;
	    std::vector<TPGFEDataformat::TcRawDataPacket> tcpktarr;
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
		// el[0] = elinks[0];
		// el[1] = elinks[1];
		// el[2] = elinks[2];
		// if(neTx>3) el[3] = elinks[3];
	      
		//Run 111137 and 111138
		if(emp_chan==106){
		  el[0] = elinks[2];
		  el[1] = elinks[1];
		  el[2] = elinks[0];
		  el[3] = elinks[7];		  
		}else{
		  el[0] = elinks[6];
		  el[1] = elinks[5];
		  el[2] = elinks[4];
		  el[3] = elinks[3];
		}
	      
		for(unsigned iel(0);iel<neTx;iel++){
		  std::cout << "\t\t el " << std::setw(3) << iel << " = 0x"
			    << std::hex << std::setfill('0')
			    << std::setw(8) << el[iel]
			    << std::dec << std::setfill(' ')
			    << std::endl;	      		
		}
	      
		TPGFEDataformat::TcRawDataPacket rdp;
		TPGStage1Emulation::Stage1IO::convertElinksToTcRawData(TPGFEDataformat::BestC, 9, el, rdp);
		//TPGStage1Emulation::Stage1IO::convertElinksToTcRawData(TPGFEDataformat::BestC, 6, el, rdp);
		//rdp.print();

		//=====================================
		tcpktarr.push_back(rdp);	      
		//=====================================
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
	      
	    }//bx lopp
	    tpgarray[eventId] = std::make_pair(modId,tcpktarr);
	  }//module condn
	}
	
	tsh=tsh->nextSubpacketHeader();
      }
      noffecafe++;
    }

    // for(unsigned j(0);j<n64;j++) {
    //   std::cout << "Word " << std::setw(6) << j << " = 0x"
    // 		<< std::hex << std::setfill('0')
    // 		<< std::setw(16) << p[j]
    // 		<< std::dec << std::setfill(' ')
    // 	      << std::endl;
    // }
    // std::cout << std::endl;
 
    std::cout << "Processing event : " <<  iEvent++ << ", noffecafe: " << noffecafe << std::endl;
    return true;
  }

  bool runStop() {
    return true;
  }
private:
  uint64_t iEvent;
};

typedef TPGEventReaderTB25 TPGReader;
typedef DAQEventReaderTB25 DAQReader;

int main(int argc, char** argv){
  
  std::cout << "Nof arguments : " << argc << std::endl;
  if(argc < 3){
    std::cerr << argv[0] << ": no run numbers specified" << std::endl;
    return 1;
  }
  
  //Command line arg assignment
  //Assign relay and run numbers
  uint32_t runNumber(12600113);
  uint32_t sourceIdTPG(1600);
  uint32_t sourceIdDAQ(1601);
  uint32_t firstLs(1);
  
  unsigned dumpEvent(0);
  
  std::istringstream issRun(argv[1]);
  issRun >> runNumber;
  std::istringstream issLinkTPG(argv[2]);
  issLinkTPG >> sourceIdTPG;
  std::istringstream issLinkDAQ(argv[3]);
  issLinkDAQ >> sourceIdDAQ;
  std::istringstream isfirstLs(argv[4]);
  isfirstLs >> firstLs;
  
  std::map<uint64_t,TPGFEDataformat::HRocarray> hrocarray; //event,rocdata_array
  std::map<uint64_t,TPGFEDataformat::TcModuleBxPackets> tpgarray; //event,tpgdata_array
  
  uint32_t zside = 0, sector = 0, link = 1, det = 0;
  uint32_t econd = 0, econt = 0, selTC4 = 1, module = 0;
  TPGFEConfiguration::TPGFEIdPacking pck;
  
  uint32_t testmodid = pck.packModId(zside, sector, link, det, econt, selTC4, module); //we assume same ECONT and ECOND number for a given module
  /////============== Reading DAQ block ========================================================
  std::vector<Hgcal10gLinkReceiver::OrbitReaderEvent> vEventsDAQ;
  
  DAQReader daqReader;
  assert(daqReader.runStart(runNumber,sourceIdDAQ));
  
  std::string oDir("dat/");
  
  Hgcal10gLinkReceiver::OrbitReader oReader;
  oReader.setPrint(false);
  
  uint64_t nEventsDAQ(0);
  
  const Hgcal10gLinkReceiver::OrbitHeader *oh;
  //Hgcal10gLinkReceiver::FragmentTrailer *ft;
  
  bool done(false);
  for(unsigned i(firstLs);!done;i++) {
    std::ostringstream oss;
    oss << oDir << std::setfill('0') << "run" << std::setw(6) << runNumber
	<< "/run" << std::setw(6) << runNumber << "_ls" << std::setw(4)
	<< i << "_index000000_source" << std::setw(4) << sourceIdDAQ << ".raw";
    
    if(!oReader.open(oss.str())) {
      if(i==firstLs) std::cout << "Failed " << oss.str() << std::endl << std::endl;
      else           std::cout << "No more files to open" << std::endl << std::endl;
      done=true;
      continue;
    }
    std::cout << "Opened " << oss.str() << std::endl << std::endl;
    
    while((oh=oReader.readOrbit(vEventsDAQ))!=nullptr and nEventsDAQ<4) {
      assert(daqReader.orbit(*oh));      
      for(unsigned j(0);j<vEventsDAQ.size();j++) {
	assert(daqReader.event(vEventsDAQ[j],nEventsDAQ,testmodid,hrocarray));
	nEventsDAQ++;
      }      
    }
    
    oReader.close();
  }

  std::cout << "Total number DAQ of events in run = " << nEventsDAQ << std::endl;

  assert(daqReader.runStop());
  /////============== Reading DAQ block ========================================================


  /////============== Reading TPG block ========================================================
  std::vector<Hgcal10gLinkReceiver::OrbitReaderEvent> vEventsTPG;
  
  TPGReader tpgReader;
  assert(tpgReader.runStart(runNumber,sourceIdTPG));

  //These are alreader defined above in DAQ block
  // std::string oDir("dat/");  
  // Hgcal10gLinkReceiver::OrbitReader oReader;
  // oReader.setPrint(false);
  
  uint64_t nEventsTPG(0);

  //These are alreader defined above in DAQ block
  // const Hgcal10gLinkReceiver::OrbitHeader *oh;
  // //Hgcal10gLinkReceiver::FragmentTrailer *ft;
  
  done = false;
  for(unsigned i(firstLs);!done;i++) {
    std::ostringstream oss;
    oss << oDir << std::setfill('0') << "run" << std::setw(6) << runNumber
	<< "/run" << std::setw(6) << runNumber << "_ls" << std::setw(4)
	<< i << "_index000000_source" << std::setw(4) << sourceIdTPG << ".raw";
    
    if(!oReader.open(oss.str())) {
      if(i==firstLs) std::cout << "Failed " << oss.str() << std::endl << std::endl;
      else           std::cout << "No more files to open" << std::endl << std::endl;
      done=true;
      continue;
    }
    std::cout << "Opened " << oss.str() << std::endl << std::endl;
    
    while((oh=oReader.readOrbit(vEventsTPG))!=nullptr and nEventsTPG<4) {
      assert(tpgReader.orbit(*oh));
      
      for(unsigned j(0);j<vEventsTPG.size();j++) {
	assert(tpgReader.event(vEventsTPG[j],nEventsTPG,testmodid,tpgarray));
	nEventsTPG++;
      }
      
    }
    
    oReader.close();
  }

  std::cout << "Total number TPG of events in run = " << nEventsTPG << std::endl;

  assert(tpgReader.runStop());
  /////============== Reading TPG block ========================================================

  ////================= Dummy configuration ====================================================
  TPGFEConfiguration::Configuration cfgs;
  cfgs.setSiChMapFile("cfgmap/WaferCellMapTraces.txt");
  cfgs.setSciChMapFile("cfgmap/channels_sipmontile_TB2024.tsv");
  cfgs.initId();
  cfgs.readSiChMapping();
  cfgs.readSciChMapping();
  cfgs.loadModIdxToNameMapping();
  cfgs.loadMuxMapping();
  
  uint32_t idx = pck.packModId(zside, sector, link, det, econd, selTC4, module); //we assume same ECONT and ECOND number for a given module      
  cfgs.setModulePath(zside, sector, link, det, econd, selTC4, module);

  cfgs.setEconTFile("cfgmap/init_econt.yaml"); 
  cfgs.readEconTConfigYaml(idx);
  
  cfgs.setEconDFile("cfgmap/init_econd.yaml");
  cfgs.readEconDConfigYaml(idx);      
  
  for(uint32_t iroc = 0; iroc < 3 ; iroc++){
    uint32_t rocid_0 = pck.packRocId(zside, sector, link, det, econd, selTC4, module, iroc, 0);
    uint32_t rocid_1 = pck.packRocId(zside, sector, link, det, econd, selTC4, module, iroc, 1);
    cfgs.setRocFile("cfgmap/configs_v3b_full/320MLF3WXIH0014_roc0_e0.yaml"); 
    cfgs.readRocConfigYaml(rocid_0, rocid_1); //only roc0 is active corresponding to chip 3	
  }

  ///=========================== Modifying the parameters for emulation ======================
  cfgs.setPedThZero();
  cfgs.setPedZero();
  cfgs.printCfgPedTh(testmodid);
  std::map<uint32_t,TPGFEConfiguration::ConfigEconD>& econDPar =  cfgs.getEconDPara();
  for(const auto& it : econDPar){
    std::cout << "Econ-D:: it.first: "<< it.first << std::endl;
    econDPar.at(it.first).setPassThrough(true);
    econDPar.at(it.first).setNeRx(6);
    econDPar.at(it.first).print();
  }
  
  uint32_t density(0);
  uint32_t droplsb(1);
  uint32_t nElinks(4);
  uint32_t select(2); //2 for BC
  std::map<uint32_t,TPGFEConfiguration::ConfigEconT>& econTPar =  cfgs.getEconTPara();
  for(const auto& it : econTPar){
    std::cout << "Econ-T:: it.first: "<< it.first << std::endl;
    econTPar[it.first].setDensity(density);
    econTPar[it.first].setDropLSB(droplsb);
    econTPar[it.first].setNElinks(nElinks);
    econTPar[it.first].setSelect(select);
    bool usesum = (econd==0 and econt==0 and link==0)?false:true;
    econTPar[it.first].setMSSumType(usesum);
    econTPar[it.first].print();
  }
  ///=========================== Modifying the parameters for emulation ======================
  ////================= Dummy configuration ====================================================

  uint64_t refEvent = 0 ;
  uint32_t refModId = testmodid;
  /////============== Emulation block ========================================================
  TPGFEModuleEmulation::HGCROCTPGEmulation rocTPGEmul(cfgs);
  TPGFEModuleEmulation::ECONTEmulation econtEmul(cfgs);
  std::map<uint64_t,std::vector<std::pair<uint32_t,TPGFEDataformat::ModuleTcData>>> modarray; //event,moduleId
  modarray.clear();
  std::cout << "hrocarray.size(): " << hrocarray.size() << std::endl;
  std::cout << "tpgarray.size(): " << tpgarray.size() << std::endl;  
  for(const auto& hrocevent : hrocarray){
    uint64_t event =  hrocevent.first;
    std::vector<std::pair<uint32_t,TPGFEDataformat::HalfHgcrocData>> hrocvec = hrocevent.second;
    std::map<uint32_t,TPGFEDataformat::HalfHgcrocData> rocdata;
    for(const auto& data : hrocvec){
      rocdata[data.first] = data.second ;
      const TPGFEDataformat::HalfHgcrocData& hrocdata = data.second ;
      std::cout << "iEvent: " << event << ", id: " << data.first << std::endl;
      hrocdata.print();
    }

    uint32_t moduleId = testmodid;
    std::pair<uint32_t,TPGFEDataformat::ModuleTcData> modTcdata;   
    std::map<uint32_t,TPGFEDataformat::ModuleTcData> moddata;

    //================================================
    //HGCROC emulation for a given module
    //================================================
    bool isSim = false; //true for CMSSW simulation and false for beam-test analysis
    //if(event==refEvent)
    //rocTPGEmul.Emulate(isSim, event, moduleId, rocdata, modTcdata, event);
    // else
    rocTPGEmul.Emulate(isSim, event, moduleId, rocdata, modTcdata);
    //================================================

    bool hasModTCshowed = false;
    //refEvent = event;
    modarray[event].push_back(modTcdata);
    if(event==refEvent){
      uint32_t modTcid = modTcdata.first;
      TPGFEDataformat::ModuleTcData& mtcdata = modTcdata.second;
      if(modTcid==refModId) {
      //if(modTcid==refModId and !mtcdata.isTcTp1() and !mtcdata.isTcTp2()) {
	for(int i=0;i<70;i++) std::cout << "=";
	std::cout << "Event: " << event ;
	for(int i=0;i<70;i++) std::cout << "=";
	std::cout << std::endl;
	std::cout << "Event: " << event << ", modTcid: " << modTcid << ", slink BxCounter: " << rocdata[modTcid].getSlinkBx() << std::endl;
	mtcdata.print();
	hasModTCshowed = true;
      }
    }

    moddata.clear();
    for(const auto& data : modarray.at(event))
      if(data.first==moduleId) moddata[data.first] = data.second ;
    
    //================================================
    //ECONT emulation for a given module
    //================================================
    econtEmul.Emulate(isSim, event, moduleId, moddata);
    TPGFEDataformat::TcModulePacket& TcRawdata = econtEmul.accessTcRawDataPacket();
    //================================================
    if(TcRawdata.second.type()==TPGFEDataformat::BestC) TcRawdata.second.sortCh();

    if(event==refEvent) {
      std::cout << "============+++++++++Emulation Begin: Event: " << event <<"++++++++================================================" << std::endl;
      TcRawdata.second.print();
      std::cout << "============+++++++++Emulation End: Event: " << event <<"++++++================================================" << std::endl;
    }

    std::cout << "============+++++++++ECON-T Begin: Event: " << event <<"++++++++================================================" << std::endl;
    for(const auto& data : tpgarray[event].second){
      data.print();
      //if(data.bx()==TcRawdata.second.bx()) data.print();
    }
    std::cout << "============+++++++++ECON-T End: Event: " << event <<"++++++++================================================" << std::endl;
    
  }//loop over roc array
  /////============== Emulation block ========================================================
  
  return 0;
}

