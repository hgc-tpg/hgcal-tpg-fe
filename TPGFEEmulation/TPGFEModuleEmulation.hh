#ifndef TPGFEModuleEmulation_h
#define TPGFEModuleEmulation_h

#include <iostream>
#include <cstring>
#include <vector>
#include <cassert>
#include <cstddef>

#include "TPGFEConfiguration.hh"

namespace TPGFEModuleEmulation{
  
  //The following emulation work per event per module
  class HGCROCTPGEmulation{
  public:
    HGCROCTPGEmulation(TPGFEConfiguration::Configuration& cfgs) : configs(cfgs) {}
    //The following emulation function performs 1) pedestal subtraction, 2) linearization, 3) compression
    void Emulate(bool isSim, uint64_t ievent, uint32_t& moduleId, std::map<uint32_t,TPGFEDataformat::HalfHgcrocData>&, std::pair<uint32_t,TPGFEDataformat::ModuleTcData>&, uint64_t);
    
    uint16_t CompressHgroc(uint32_t val, bool isldm){
      
      //////The configuration of half of ROC based on HGCROC3a [doc. no. v2.0] (Read subsection 1.3.3 of Page-39)
      //////EDMS ROCv3a: https://edms.cern.ch/ui/#!master/navigator/document?D:100570166:100570166:subDocs
      //4E+3M
      //Ref:https://graphics.stanford.edu/%7Eseander/bithacks.html

      //controlled by the SelTC4 parameter of HGCROC (1/0 for LD/HD)
      val = (isldm)?val>>1:val>>3;
            
      uint32_t r = 0; // r will be lg(v)
      uint32_t sub ;
      uint32_t shift ;
      uint32_t mant ;
      
      if(val>7){
	uint32_t v = val; 
	r = 0; 
	while (v >>= 1) r++;
	sub = r - 2;
	shift = r - 3;
	if(sub<=0xF){
	  mant = (val>>shift) & 0x7;
	}else{
	  sub = 0xF;
	  mant = 0x7;
	}
      }else{
	r = 0;
	sub = 0;
	shift = 0;
	mant = val & 0x7;
      }
  
      uint16_t cdata = (sub<<3) | mant;
  
      return cdata;
    }
    
  private:    
    TPGFEConfiguration::Configuration& configs;
    TPGFEConfiguration::TPGFEIdPacking pck;
  };
  
  void HGCROCTPGEmulation::Emulate(bool isSim, uint64_t ievent, uint32_t& moduleId, std::map<uint32_t,TPGFEDataformat::HalfHgcrocData>& rocdata, std::pair<uint32_t,TPGFEDataformat::ModuleTcData>& moddata, uint64_t refEvent = 0xFFFFFFFFFFFFFFFF){

    pck.setModId(moduleId);    
    const std::map<std::tuple<uint32_t,uint32_t,uint32_t>,std::string>& modNameMap = configs.getModIdxToName();
    const std::string& modName = modNameMap.at(std::make_tuple(pck.getDetType(),pck.getSelTC4(),pck.getModule()));
    if(ievent==refEvent) std::cout << "TPGFEModuleEmulation::HGCROCTPGEmulation::Emulate: modName: " << modName << std::endl;
    if(ievent==refEvent) std::cout << "TPGFEModuleEmulation::HGCROCTPGEmulation::Emulate: pck.getDetType(): " << pck.getDetType() << std::endl;
    const std::map<std::string,std::vector<uint32_t>>& modTClist = (pck.getDetType()==0)?configs.getSiModTClist():configs.getSciModTClist();
    const std::vector<uint32_t>& tclist = modTClist.at(modName) ;
    const std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>& tcPinMap = (pck.getDetType()==0)?configs.getSiTCToROCpin():configs.getSciTCToROCpin();

    uint16_t bx = 0xFFFF;
    const uint32_t nTCs = tclist.size();
    if(ievent==refEvent) std::cout << "TPGFEModuleEmulation::HGCROCTPGEmulation::Emulate: nTCs : " << nTCs << std::endl;
    if(ievent==refEvent) for(const auto& itc : tclist) std::cout << "TPGFEModuleEmulation::HGCROCTPGEmulation::Emulate: tclist["<<itc<<"] " << std::endl;
    TPGFEDataformat::ModuleTcData mdata;
    mdata.setNofTCs(nTCs);
    TPGFEDataformat::HgcrocTcData hrtcdata[nTCs];
    
    for(const auto& itc : tclist){
      const std::vector<uint32_t>& pinlist = tcPinMap.at(std::make_pair(modName,itc)) ;
      uint32_t totadc = 0;
      bool isTot = false;
      bool isTcTp1 = false;
      bool isTcTp2 = false;
      bool isTcTp3 = false;
      for(const auto& tcch : pinlist){
	uint32_t rocpin = tcch%36 ;
	uint32_t rocn = TMath::Floor(tcch/72);
	uint32_t half = (int(TMath::Floor(tcch/36))%2==0)?0:1;
	uint32_t rocid = pck.getRocIdFromModId(moduleId,rocn,half);
	if(rocdata.find(rocid)==rocdata.end()){
	  std::cerr << "TPGFEModuleEmulation::HGCROCTPGEmulation::Emulate: HalfRoc not found in data for Event "<< ievent << ", rocid: " << rocid <<", Tcch: "<< tcch << ", rocn: " << rocn << ", half: " << half << ", rocpin: " << rocpin << std::endl;
	  continue ; 
	}
	if(ievent==refEvent) std::cout<<"TPGFEModuleEmulation::HGCROCTPGEmulation::Emulate: TC : " << itc <<", rocid: "<<rocid<< ", tcch: " << tcch <<", rocpin : "<<rocpin<<", rocn: "<<rocn<<", half: "<<half<<std::endl;
	TPGFEDataformat::HalfHgcrocChannelData& chdata = rocdata.at(rocid).getChannelData(rocpin);
	bx = rocdata.at(rocid).getBx();
	if(chdata.getTcTp()==1) isTcTp1 = true;
	if(chdata.getTcTp()==2) isTcTp2 = true;
	if(chdata.getTcTp()==3) isTcTp3 = true;
	if(ievent==refEvent) chdata.print();
	if(!isSim){ //if not simulation, assumed beamtest data
	  const TPGFEConfiguration::ConfigHfROC& rocpara = configs.getRocPara().at(rocid);
	  uint32_t ped = configs.getChPara().at(pck.packChId(rocid,rocpin)).getAdcpedestal();
	  //if(ievent==refEvent) rocpara.print();
	  if(!chdata.isTot()){
	    unsigned thr = rocpara.getAdcTH();
	    uint32_t adc = chdata.getAdc();
	    // if(chdata.getTcTp()==1)
	    //   if(itc!=20 or rocid!=770) adc = 0;
	    adc = (adc>(ped+thr) and !(rocpara.isChMasked(rocpin)) and (ped<0xFF)) ? adc-ped : 0 ;
	    totadc += adc;
	    if(ievent==refEvent) std::cout<<"\t TPGFEModuleEmulation::HGCROCTPGEmulation::Emulate: ped: " << ped << ", thr: " << thr <<", adcs0 : "<<chdata.getAdc()<<", adcs1 : "<<adc<<", rocpara.isChMasked(rocpin): "<< rocpara.isChMasked(rocpin) <<std::endl;
	  }else{
	    uint32_t tot = chdata.getTot(); //+7 for totup
	    uint32_t tot1 = (tot>=rocpara.getTotTH(rocpin)) ? (tot-rocpara.getTotP(rocpin)) : (rocpara.getTotTH(rocpin)-rocpara.getTotP(rocpin)) ;
	    //ideally ped<255 condition should only be restricted for ADC case, however TOT signal triggers in emulation for event 146245 of relay 1695829026 and link 1, which is not that seen by ECONT data
	    uint32_t totlin = (!(rocpara.isChMasked(rocpin)) and (ped<0xFF))?tot1*rocpara.getMultFactor():0; 
	    totadc += totlin;
	    isTot = true;
	    if(ievent==refEvent) std::cout<<"\t TPGFEModuleEmulation::HGCROCTPGEmulation::Emulate: tot1: " << tot1 << ", tot: " << chdata.getTot() <<", thr : "<< rocpara.getTotTH(rocpin) <<", ped: "<< rocpara.getTotP(rocpin) << ", totlin : "<< totlin <<std::endl;
	  }//istot or adc
	}else{
	  //check CMSSW for details should be 
	  if(!chdata.isTot())
	    totadc += chdata.getAdc(); 
	  else
	    totadc += chdata.getTot(); 
	}//isCMSSW simulation of beam-test analysis
      }//loop of sensor channels
      hrtcdata[itc].setCharge(totadc);
      hrtcdata[itc].setTot(isTot); //stand for tctp==3
      hrtcdata[itc].setTcTp1(isTcTp1);
      hrtcdata[itc].setTcTp2(isTcTp2);
      hrtcdata[itc].setTcTp3(isTcTp3);
      if(!isSim){
	pck.setModId(moduleId);    
	hrtcdata[itc].setCdata(CompressHgroc(totadc, pck.getSelTC4()));
      }else
	hrtcdata[itc].setCdata(CompressHgroc(totadc,1));
    }//loop over TCs
    mdata.setTCs(hrtcdata);
    mdata.setBx(bx);
    moddata = std::make_pair(moduleId,mdata);
  }

  class ECONTEmulation{
  public:
    ECONTEmulation(TPGFEConfiguration::Configuration& cfgs) : configs(cfgs), isVerbose(false) {}
    //~ECONTEmulation() {tcrawlist.clear();} Not needed
    ~ECONTEmulation() {}
    
    //The following emulation function performs 1) decompression, 2) calibration, 3) compression
    void Emulate(bool isSim, uint64_t ievent, uint32_t& moduleId, const std::map<uint32_t,TPGFEDataformat::ModuleTcData>&, TPGFEDataformat::TcModulePacket& );
    void EmulateJulTB(bool isSim, uint64_t ievent, uint32_t& moduleId, const std::map<uint32_t,TPGFEDataformat::ModuleTcData>&, TPGFEDataformat::TcModulePacket& );
    void EmulateSTC(bool isSim, uint64_t ievent, uint32_t& moduleId, const std::map<uint32_t,TPGFEDataformat::ModuleTcData>&);
    void EmulateJulSTC(bool isSim, uint64_t ievent, uint32_t& moduleId, const std::map<uint32_t,TPGFEDataformat::ModuleTcData>&);
    void EmulateBC(bool isSim, uint64_t ievent, uint32_t& moduleId, const std::map<uint32_t,TPGFEDataformat::ModuleTcData>&);
    
    void Emulate(bool isSim, uint64_t ievent, uint32_t& moduleId, const std::map<uint32_t,TPGFEDataformat::ModuleTcData>&);
    void EmulateJulTB(bool isSim, uint64_t ievent, uint32_t& moduleId, const std::map<uint32_t,TPGFEDataformat::ModuleTcData>&);
    const TPGFEDataformat::TcModulePacket& getTcRawDataPacket() const {return emulOut;}
    TPGFEDataformat::TcModulePacket& accessTcRawDataPacket() {return emulOut;}
    void setVerbose(bool verbose = true) {isVerbose = verbose;}
    
    // BX here may be absorbed into TcRawDataPacket?
    static void convertToElinkData(unsigned bx, const TPGFEDataformat::TcRawDataPacket &tcrdp, uint32_t *ve);
    void getElinkData(unsigned bx, uint32_t *ptr) const;

    void fillRandomTcRawData(unsigned bx);
    void fillZeroEnergyTcRawData(unsigned bx);
    
    static void generateTcRawData(bool zero, unsigned bx,
				  TPGFEDataformat::Type type,
				  unsigned nTc,
				  TPGFEDataformat::TcRawDataPacket &vtcrp);
    static void generateRandomTcRawData(unsigned bx,
					TPGFEDataformat::Type type,
					unsigned nTc,
					TPGFEDataformat::TcRawDataPacket &vtcrp);
    static void generateZeroEnergyTcRawData(unsigned bx,
					    TPGFEDataformat::Type type,
					    unsigned nTc,
					    TPGFEDataformat::TcRawDataPacket &vtcrp);
    
  private:
    uint32_t DecompressEcont(uint16_t compressed, bool density){
      //4E+3M with midpoint correction
      //controlled by the density parameter of ECON-T (0/1 stands for LD/HD and represents 1/3 bit shifts)
      uint32_t mant = compressed & 0x7;
      uint32_t expo = (compressed>>3) & 0xF;
      
      if(expo==0) 
	return (density) ? (mant<<3)+4 : (mant<<1)+1 ; //The the +4/+1 are midpoint corrections for expo==0
      
      uint32_t shift = expo+2; 
               shift += (density) ? 3 : 1;
      uint32_t decomp = 1<<shift; 
      uint32_t mpdeco = 1<<(shift-4);
      decomp = decomp | (mant<<(shift-3));
      decomp = decomp | mpdeco;

      return decomp;
    }

    uint16_t CompressEcontStc4E3M(uint64_t val, uint32_t dropLSB){
      //4E+3M
      //controlled by the dropLSB parameter of ECON-T 
      val = val>>dropLSB;  
      
      uint32_t r = 0; // r will be lg(v)
      uint32_t sub ;
      uint32_t shift ;
      uint32_t mant ;
      
      if(val>7){
	uint64_t v = val; 
	r = 0; 
	while (v >>= 1) r++;
	sub = r - 2;
	shift = r - 3;
	if(sub<=0xF){
	  mant = (val>>shift) & 0x7;
	}else{
	  sub = 0xF;
	  mant = 0x7;
	}
      }else{
	r = 0;
	sub = 0;
	shift = 0;
	mant = val & 0x7;
      }
      
      sub = sub & 0xF;
      mant = mant & 0x7;
      
      uint16_t packed = (sub<<3) | mant;
      
      return packed;
    }
    
    uint16_t CompressEcontStc5E4M(uint64_t val){
      //5E+4M
      //dropLSB is not applicable for 5E+4M

      uint32_t r = 0; // r will be lg(v)
      uint32_t sub ;
      uint32_t shift ;
      uint32_t mant ;
      
      if(val>0xF){
	uint64_t v = val; 
	r = 0; 
	while (v >>= 1) r++;
	sub = r - 3;         
	shift = r - 4;       
	if(sub<=0x1F){
	  mant = (val>>shift) & 0xF;
	}else{
	  sub = 0x1F;
	  mant = 0xF;
	}
      }else{
	r = 0;
	sub = 0;
	shift = 0;
	mant = val & 0xF;
      }

      sub = sub & 0x1F;
      mant = mant & 0xF;

      uint16_t packed = (sub<<4) | mant;
      
      return packed;
    }
    
    uint16_t CompressEcontBc(uint64_t val, uint32_t dropLSB){
      //4E+3M
      //controlled by the dropLSB parameter of ECON-T 
      val = val>>dropLSB;  
      
      uint32_t r = 0; // r will be lg(v)
      uint32_t sub ;
      uint32_t shift ;
      uint32_t mant ;
  
      if(val>7){
	uint64_t v = val; 
	r = 0; 
	while (v >>= 1) r++;
	sub = r - 2;
	shift = r - 3;
	if(sub<=0xF){
	  mant = (val>>shift) & 0x7;
	}else{
	  sub = 0xF;
	  mant = 0x7;
	}
      }else{
	r = 0;
	sub = 0;
	shift = 0;
	mant = val & 0x7;
      }
      
      sub = sub & 0xF;
      mant = mant & 0x7;

      uint16_t packed = (sub<<3) | mant;
      
      return packed;
    }
    
    uint16_t CompressEcontModsum(uint64_t val, uint32_t dropLSB){
      //5E+3M
      //controlled by the dropLSB parameter of ECON-T 
      //val = val>>dropLSB;  

      uint32_t r = 0; // r will be lg(v)
      uint32_t sub ;
      uint32_t shift ;
      uint32_t mant ;
      
      if(val>7){
	uint64_t v = val; 
	r = 0; 
	while (v >>= 1) r++;
	sub = r - 2;
	shift = r - 3;
	if(sub<=0x1F){
	  mant = (val>>shift) & 0x7;
	}else{
	  sub = 0x1F;
	  mant = 0x7;
	}
      }else{
	r = 0;
	sub = 0;
	shift = 0;
	mant = val & 0x7;
      }
      
      sub = sub & 0x1F;
      mant = mant & 0x7;
      
      uint16_t packed = (sub<<3) | mant;
 
      return packed;
    }
    

    void batcherOEMSort(std::vector<TPGFEDataformat::TcRawData>& tc) {
      //C++ adaptation of batcher odd-even sorting of https://github.com/dnoonan08/ECONT_Emulator/ASICBlocks/bestchoice.py
      //Requires 49 Tcs for sorting
      TPGFEDataformat::TcRawData dummy(TPGFEDataformat::Unknown,0x3f, 0, 0, false, false, false);
      tc.push_back(dummy);
      
      uint32_t N = uint32_t(tc.size());
      uint32_t t = uint32_t(ceil(log(N)/log(2)));
      uint32_t p = uint32_t(pow(2,(t-1)));
      while(p>=1){
	uint32_t q = uint32_t(pow(2,(t-1)));
	uint32_t r = 0;
	uint32_t d = p;
	while (q>=p){
	  for(uint32_t i=0; i<(N-d) ; i++){
	    if ((i & p) != r) continue;
	    if (tc[i] < tc[i+d]) std::swap(tc[i], tc[i+d]);
	  }
	  d = q - p;
	  q = floor(q/2);
	  r = p;
	}
	p = floor(p/2);
      }
    }

    uint32_t findLastMax(std::vector<TPGFEDataformat::TcRawData>& tc) {
      for(uint32_t itc = 0 ; itc< uint32_t(tc.size()-1); itc++) if(tc[itc+1] < tc[0]) return itc;
      return uint32_t(tc.size()-1);
    }

    bool isTie(std::vector<TPGFEDataformat::TcRawData>& tc) {
      for(uint32_t itc = 0 ; itc< tc.size(); itc++) if(tc[itc].energy() !=  tc[0].energy()) return false;
      return true;
    }

    bool isAllZero(std::vector<TPGFEDataformat::TcRawData>& tc) {
      for(uint32_t itc = 0 ; itc< tc.size(); itc++) if(tc[itc].energy()>0) return false;
      return true;
    }


    TPGFEConfiguration::Configuration& configs;
    TPGFEConfiguration::TPGFEIdPacking pck;
    TPGFEDataformat::TcModulePacket emulOut;
    bool isVerbose ; 
  };

  void ECONTEmulation::Emulate(bool isSim, uint64_t ievent, uint32_t& moduleId, const std::map<uint32_t,TPGFEDataformat::ModuleTcData>& moddata, TPGFEDataformat::TcModulePacket& econtOut){
    Emulate(isSim, ievent, moduleId, moddata);
    econtOut=emulOut;
  }

  void ECONTEmulation::EmulateJulTB(bool isSim, uint64_t ievent, uint32_t& moduleId, const std::map<uint32_t,TPGFEDataformat::ModuleTcData>& moddata, TPGFEDataformat::TcModulePacket& econtOut){
    EmulateJulTB(isSim, ievent, moduleId, moddata);
    econtOut=emulOut;
  }

  
  void ECONTEmulation::Emulate(bool isSim, uint64_t ievent, uint32_t& moduleId, const std::map<uint32_t,TPGFEDataformat::ModuleTcData>& moddata){
    emulOut.first=moduleId;
    emulOut.second.reset();
    
    const TPGFEDataformat::Type& outputType = configs.getEconTPara().at(moduleId).getOutType();
    if(outputType==TPGFEDataformat::BestC)
      EmulateBC(isSim, ievent, moduleId, moddata);
    else if(outputType==TPGFEDataformat::STC4A or outputType==TPGFEDataformat::STC4B or outputType==TPGFEDataformat::STC16 or TPGFEDataformat::CTC4A or TPGFEDataformat::CTC4B)
      EmulateSTC(isSim, ievent, moduleId, moddata);
  }

  void ECONTEmulation::EmulateJulTB(bool isSim, uint64_t ievent, uint32_t& moduleId, const std::map<uint32_t,TPGFEDataformat::ModuleTcData>& moddata){
    emulOut.first=moduleId;
    emulOut.second.reset();
    
    const TPGFEDataformat::Type& outputType = configs.getEconTPara().at(moduleId).getOutType();
    if(outputType==TPGFEDataformat::BestC)
      EmulateBC(isSim, ievent, moduleId, moddata);
    else if(outputType==TPGFEDataformat::STC4A or outputType==TPGFEDataformat::STC4B or outputType==TPGFEDataformat::STC16 or TPGFEDataformat::CTC4A or TPGFEDataformat::CTC4B)
      EmulateJulSTC(isSim, ievent, moduleId, moddata);
  }
  
  void ECONTEmulation::EmulateSTC(bool isSim, uint64_t ievent, uint32_t& moduleId, const std::map<uint32_t,TPGFEDataformat::ModuleTcData>& moddata) {

    pck.setModId(moduleId);    
    const std::map<std::tuple<uint32_t,uint32_t,uint32_t>,std::string>& modNameMap = configs.getModIdxToName();
    const std::map<uint32_t,uint32_t>& refMuxMap = configs.getMuxMapping() ;
    const std::string& modName = modNameMap.at(std::make_tuple(pck.getDetType(),pck.getSelTC4(),pck.getModule()));
    uint32_t dropLSB = configs.getEconTPara().at(moduleId).getDropLSB() ;
    uint32_t nofSTCs = configs.getEconTPara().at(moduleId).getNofSTCs();
    TPGFEDataformat::TcRawData dummy(TPGFEDataformat::Unknown,0x3f, 0, 0, false, false, false);
    
    const TPGFEDataformat::Type& outputType = configs.getEconTPara().at(moduleId).getOutType();
   
    if(outputType==TPGFEDataformat::STC4A or outputType==TPGFEDataformat::STC4B or outputType==TPGFEDataformat::CTC4A or outputType==TPGFEDataformat::CTC4B){
      
      const std::map<std::string,std::vector<uint32_t>>& modSTClist = (pck.getDetType()==0)?configs.getSiModSTClist():configs.getSciModSTClist();
      const std::vector<uint32_t>& stclist = modSTClist.at(modName) ;
      const std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>& stcTcMap = (pck.getDetType()==0)?configs.getSiSTCToTC():configs.getSciSTCToTC();
      uint16_t bx = 0xffff;
      emulOut.second.reset();
      for(const auto& istc : stclist){
	bool isTcTp1 = false;
	bool isTcTp2 = false;
	bool isTcTp3 = false;
	if(istc>=nofSTCs) continue;
	const std::vector<uint32_t>& tclist = stcTcMap.at(std::make_pair(modName,istc));
	const TPGFEDataformat::ModuleTcData& mdata = moddata.at(moduleId);
	bx = (mdata.getBx()==3564) ? 0xF : mdata.getBx() & 0x7 ; 
	uint64_t decompressedSTC = 0;
	std::vector<TPGFEDataformat::TcRawData> tcrawdatalist;
	TPGFEDataformat::TcRawData tcdata;
	for(const auto& econtc : tclist){
	  uint32_t hgctc = configs.getEconTPara().at(moduleId).getInputMux(econtc) ;
	  bool hasFound = false;
	  uint32_t emultc = 0xffffffff;
	  for (const auto& it : refMuxMap){
	    if (it.second == hgctc)  {
	      hasFound = true;
	      emultc = it.first;
	    }
	  }
	  if(!hasFound){
	    std::cerr << "TPGFEModuleEmulation::ECONTEmulation::EmulateSTC (moduleid="<<moduleId<<") : Mux not set for TC " << econtc << std::endl;
	    continue;
	  }
	  if(mdata.getTC(emultc).isTcTp1()) isTcTp1 = true;
	  if(mdata.getTC(emultc).isTcTp2()) isTcTp2 = true;
	  if(mdata.getTC(emultc).isTcTp3()) isTcTp3 = true;
	  uint64_t decompressed = DecompressEcont(mdata.getTC(emultc).getCdata(),configs.getEconTPara().at(moduleId).getDensity());
	  uint64_t decomp64bit =  decompressed * configs.getEconTPara().at(moduleId).getCalibration(econtc) ;
	  if(isVerbose) std::cout << "TPGFEModuleEmulation::ECONTEmulation::EmulateSTC4 (moduleid="<<moduleId<<", TC="<<econtc<<", hgctc="<<hgctc<<") decompressed: " << decompressed << ", decompressed*calib: " << decomp64bit << std::endl; 
	  decompressed =  decomp64bit >> 11;
	  decompressedSTC += decompressed  ;
	  if(isVerbose) std::cout << "TPGFEModuleEmulation::ECONTEmulation::EmulateSTC4 (moduleid="<<moduleId<<", TC="<<econtc<<") decompressed*calib>>11: " << decompressed << ", decompressedSTC: " << decompressedSTC << std::endl; 
	  uint16_t compressed_bc = CompressEcontBc(decompressed, dropLSB);
	  tcdata.setTriggerCell(TPGFEDataformat::BestC, econtc, compressed_bc, decompressed, mdata.getTC(emultc).isTcTp1(), mdata.getTC(emultc).isTcTp2(), mdata.getTC(emultc).isTcTp3()) ;
	  tcrawdatalist.push_back(tcdata);
	}
	if(isVerbose) std::cout << "TPGFEModuleEmulation::ECONTEmulation::EmulateSTC4 (moduleid="<<moduleId<<", STC="<< istc <<") tcrawdatalist.size: " << tcrawdatalist.size() << std::endl;
	std::sort(tcrawdatalist.begin(),tcrawdatalist.end(),TPGFEDataformat::TcRawDataPacket::customGTE);
	uint32_t lMaxId = isTie(tcrawdatalist) ? 3 : findLastMax(tcrawdatalist);
	tcrawdatalist.resize(lMaxId+1);
	if(isVerbose) std::cout << "TPGFEModuleEmulation::ECONTEmulation::EmulateSTC4 Before sorting....(moduleid="<<moduleId<<", STC="<< istc <<"), modified size: " << tcrawdatalist.size() << std::endl;
	if(isVerbose) for(int i=0;i<tcrawdatalist.size();i++) tcrawdatalist[i].print();
	std::sort(tcrawdatalist.begin(),tcrawdatalist.end(),TPGFEDataformat::TcRawDataPacket::customLTA);
	if(isVerbose) std::cout << "TPGFEModuleEmulation::ECONTEmulation::EmulateSTC4 After sorting....(moduleid="<<moduleId<<", STC="<< istc <<"), modified size: " << tcrawdatalist.size() << std::endl;
	if(isVerbose) for(int i=0;i<tcrawdatalist.size();i++) tcrawdatalist[i].print();
	TPGFEDataformat::TcRawData firstMinA = tcrawdatalist[0] ;
	if(isVerbose) firstMinA.print();
	uint16_t compressed_energy = (outputType==TPGFEDataformat::STC4A or outputType==TPGFEDataformat::CTC4A) ? CompressEcontStc4E3M(decompressedSTC, dropLSB) : CompressEcontStc5E4M(decompressedSTC);
	uint64_t raw_E = (outputType==TPGFEDataformat::STC4A or outputType==TPGFEDataformat::CTC4A) ? decompressedSTC>>dropLSB : decompressedSTC ;
	emulOut.second.setTBM(outputType, bx, 0); 
	if(outputType==TPGFEDataformat::STC4A or outputType==TPGFEDataformat::STC4B)
	  emulOut.second.setTcData(outputType, firstMinA.address()%4, compressed_energy, raw_E, isTcTp1, isTcTp2, isTcTp3);
	else
	  emulOut.second.setTcData(outputType, istc, compressed_energy, raw_E, isTcTp1, isTcTp2, isTcTp3);
	tcrawdatalist.clear();
      }//stc loop
      
    }
    
    if(outputType==TPGFEDataformat::STC16){
      
      const std::map<std::string,std::vector<uint32_t>>& modSTC16list = (pck.getDetType()==0)?configs.getSiModSTC16list():configs.getSciModSTC16list();
      const std::vector<uint32_t>& stc16list = modSTC16list.at(modName) ;
      const std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>& stc16TcMap = (pck.getDetType()==0)?configs.getSiSTC16ToTC():configs.getSciSTC16ToTC();
      uint16_t bx = 0xffff;
      emulOut.second.reset();
      for(const auto& istc16 : stc16list){
	bool isTcTp1 = false;
	bool isTcTp2 = false;
	bool isTcTp3 = false;
	//std::cout<<"iSTC: " << istc16 << std::endl;
	if(istc16>=nofSTCs) continue;
	const std::vector<uint32_t>& tclist = stc16TcMap.at(std::make_pair(modName,istc16));
	const TPGFEDataformat::ModuleTcData& mdata = moddata.at(moduleId);
	bx = (mdata.getBx()==3564) ? 0xF : mdata.getBx() & 0x7 ; 
	uint64_t decompressedSTC16 = 0;
	std::vector<TPGFEDataformat::TcRawData> tcrawdatalist;
	TPGFEDataformat::TcRawData tcdata;
	for(const auto& econtc : tclist){
	  uint32_t hgctc = configs.getEconTPara().at(moduleId).getInputMux(econtc) ;
	  bool hasFound = false;
	  uint32_t emultc = 0xffffffff;
	  for (const auto& it : refMuxMap){
	    //std::cout<<"iSTC: " << istc16 << ", first : " << it.first << ", it.second: " << it.second << ", econtc: " << econtc << ", hgctc: " << hgctc << std::endl;
	    if (it.second == hgctc)  {
	      hasFound = true;
	      emultc = it.first;
	    }
	  }
	  if(!hasFound){
	    std::cerr << "TPGFEModuleEmulation::ECONTEmulation::EmulateBC (moduleid="<<moduleId<<") : Mux not set for TC " << econtc << std::endl;
	    continue;
	  }// else{
	  //   std::cout << "TPGFEModuleEmulation::ECONTEmulation::EmulateSTC (moduleid="<<moduleId<<") : Mux set for TC " << econtc << " --> " << emultc << std::endl;
	  // }
	  if(mdata.getTC(emultc).isTcTp1()) isTcTp1 = true;
	  if(mdata.getTC(emultc).isTcTp2()) isTcTp2 = true;
	  if(mdata.getTC(emultc).isTcTp3()) isTcTp3 = true;
	  uint64_t decompressed = DecompressEcont(mdata.getTC(emultc).getCdata(), configs.getEconTPara().at(moduleId).getDensity());
	  uint64_t decomp64bit = decompressed * configs.getEconTPara().at(moduleId).getCalibration(econtc) ;
	  if(isVerbose) std::cout << "TPGFEModuleEmulation::ECONTEmulation::EmulateSTC16 (moduleid="<<moduleId<<", TC="<<econtc<<", hgctc="<<hgctc<<") decompressed: " << decompressed << ", decompressed*calib: " << decomp64bit << std::endl; 
	  decompressed =  decomp64bit >> 11;
	  decompressedSTC16 += decompressed ;
	  if(isVerbose) std::cout << "TPGFEModuleEmulation::ECONTEmulation::EmulateSTC16 (moduleid="<<moduleId<<", TC="<<econtc<<") decompressed*calib>>11: " << decompressed << ", decompressedSTC16: " << decompressedSTC16 << std::endl; 
	  uint16_t compressed_bc = CompressEcontBc(decompressed,dropLSB);
	  tcdata.setTriggerCell(TPGFEDataformat::BestC, econtc, compressed_bc, decompressed, mdata.getTC(emultc).isTcTp1(), mdata.getTC(emultc).isTcTp2(), mdata.getTC(emultc).isTcTp3()) ;
	  tcrawdatalist.push_back(tcdata);
	}
	// batcherOEMSort(tcrawdatalist);
	// uint32_t lMaxId = findLastMax(tcrawdatalist);
	// tcrawdatalist.resize(lMaxId+1);
	std::sort(tcrawdatalist.begin(),tcrawdatalist.end(),TPGFEDataformat::TcRawDataPacket::customGTE);
	uint32_t lMaxId = isTie(tcrawdatalist) ? 3 : findLastMax(tcrawdatalist);
	tcrawdatalist.resize(lMaxId+1);
	std::sort(tcrawdatalist.begin(),tcrawdatalist.end(),TPGFEDataformat::TcRawDataPacket::customLTA);
	TPGFEDataformat::TcRawData firstMinA = tcrawdatalist[0] ;
	uint16_t compressed_energy = CompressEcontStc5E4M(decompressedSTC16);
	emulOut.second.setTBM(outputType, bx, 0); 
	emulOut.second.setTcData(outputType, (firstMinA.address())%16, compressed_energy, decompressedSTC16, isTcTp1, isTcTp2, isTcTp3);
	tcrawdatalist.clear();
      }//stc16 loop
    }//STC16 select condition
    
  }//Emulate STC

  void ECONTEmulation::EmulateJulSTC(bool isSim, uint64_t ievent, uint32_t& moduleId, const std::map<uint32_t,TPGFEDataformat::ModuleTcData>& moddata) {
    
    pck.setModId(moduleId);    
    const std::map<std::tuple<uint32_t,uint32_t,uint32_t>,std::string>& modNameMap = configs.getModIdxToName();
    const std::map<uint32_t,uint32_t>& refMuxMap = configs.getMuxMapping() ;
    const std::string& modName = modNameMap.at(std::make_tuple(pck.getDetType(),pck.getSelTC4(),pck.getModule()));
    uint32_t dropLSB = configs.getEconTPara().at(moduleId).getDropLSB() ;
    uint32_t nofSTCs = configs.getEconTPara().at(moduleId).getNofSTCs();
    
    const TPGFEDataformat::Type& outputType = configs.getEconTPara().at(moduleId).getOutType();
   
    if(outputType==TPGFEDataformat::STC4A or outputType==TPGFEDataformat::STC4B or outputType==TPGFEDataformat::CTC4A or outputType==TPGFEDataformat::CTC4B){
      
      const std::map<std::string,std::vector<uint32_t>>& modSTClist = (pck.getDetType()==0)?configs.getSiModSTClist():configs.getSciModSTClist();
      const std::vector<uint32_t>& stclist = modSTClist.at(modName) ;
      const std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>& stcTcMap = (pck.getDetType()==0)?configs.getSiSTCToTC():configs.getSciSTCToTC();
      uint16_t bx = 0xffff;
      emulOut.second.reset();
      for(const auto& istc : stclist){
	bool isTcTp1 = false;
	bool isTcTp2 = false;
	bool isTcTp3 = false;
	if(istc>=nofSTCs) continue;
	const std::vector<uint32_t>& tclist = stcTcMap.at(std::make_pair(modName,istc));
	const TPGFEDataformat::ModuleTcData& mdata = moddata.at(moduleId);
	bx = (mdata.getBx()==3564) ? 0xF : mdata.getBx() & 0x7 ; 
	uint64_t decompressedSTC = 0;
	std::vector<TPGFEDataformat::TcRawData> tcrawdatalist;
	TPGFEDataformat::TcRawData tcdata;
	for(const auto& econtc : tclist){
	  uint32_t hgctc = configs.getEconTPara().at(moduleId).getInputMux(econtc) ;
	  bool hasFound = false;
	  uint32_t emultc = 0xffffffff;
	  for (const auto& it : refMuxMap){
	    if (it.second == hgctc)  {
	      hasFound = true;
	      emultc = it.first;
	    }
	  }
	  if(!hasFound){
	    std::cerr << "TPGFEModuleEmulation::ECONTEmulation::EmulateJulSTC (moduleid="<<moduleId<<") : Mux not set for TC " << econtc << std::endl;
	    continue;
	  }
	  if(mdata.getTC(emultc).isTcTp1()) isTcTp1 = true;
	  if(mdata.getTC(emultc).isTcTp2()) isTcTp3 = true;
	  if(mdata.getTC(emultc).isTcTp3()) isTcTp3 = true;
	  uint64_t decompressed = DecompressEcont(mdata.getTC(emultc).getCdata(),configs.getEconTPara().at(moduleId).getDensity());
	  uint64_t decomp64bit =  decompressed * configs.getEconTPara().at(moduleId).getCalibration(econtc) ;
	  if(isVerbose) std::cout << "TPGFEModuleEmulation::ECONTEmulation::EmulateJulSTC4 (moduleid="<<moduleId<<", TC="<<econtc<<") decompressed: " << decompressed << ", decompressed*calib: " << decomp64bit << std::endl; 
	  decompressed =  decomp64bit >> 11;
	  decompressedSTC += decompressed ;
	  if(isVerbose) std::cout << "TPGFEModuleEmulation::ECONTEmulation::EmulateJulSTC4 (moduleid="<<moduleId<<", TC="<<econtc<<") decompressed*calib>>11: " << decompressed << ", decompressedSTC: " << decompressedSTC << std::endl; 
	  uint16_t compressed_bc = CompressEcontBc(decompressed, dropLSB);
	  tcdata.setTriggerCell(TPGFEDataformat::BestC, econtc, compressed_bc, decompressed, mdata.getTC(emultc).isTcTp1(), mdata.getTC(emultc).isTcTp2(), mdata.getTC(emultc).isTcTp3()) ;
	  tcrawdatalist.push_back(tcdata);
	}
	batcherOEMSort(tcrawdatalist);
	uint32_t lMaxId = findLastMax(tcrawdatalist);
	tcrawdatalist.resize(lMaxId+1);
	if(isVerbose) std::cout << "TPGFEModuleEmulation::ECONTEmulation::EmulateJulSTC4 (moduleid="<<moduleId<<", STC="<< istc <<") size: " << tcrawdatalist.size() << std::endl;
	if(isVerbose) for(int i=0;i<tcrawdatalist.size();i++) tcrawdatalist[i].print();
	std::sort(tcrawdatalist.begin(),tcrawdatalist.end(),TPGFEDataformat::TcRawDataPacket::customGTA);
	TPGFEDataformat::TcRawData lastMax = isAllZero(tcrawdatalist) ? tcrawdatalist[tcrawdatalist.size()-1] : tcrawdatalist[0] ;
	if(isVerbose) lastMax.print();
	uint16_t compressed_energy = (outputType==TPGFEDataformat::STC4A or outputType==TPGFEDataformat::CTC4A) ? CompressEcontStc4E3M(decompressedSTC, dropLSB) : CompressEcontStc5E4M(decompressedSTC);
	uint64_t raw_E = (outputType==TPGFEDataformat::STC4A or outputType==TPGFEDataformat::CTC4A) ? decompressedSTC>>dropLSB : decompressedSTC ;
	emulOut.second.setTBM(outputType, bx, 0); 
	if(outputType==TPGFEDataformat::STC4A or outputType==TPGFEDataformat::STC4B)
	  emulOut.second.setTcData(outputType, lastMax.address()%4, compressed_energy, raw_E, isTcTp1, isTcTp2, isTcTp3);
	else
	  emulOut.second.setTcData(outputType, istc, compressed_energy, raw_E, isTcTp1, isTcTp2, isTcTp3);
	tcrawdatalist.clear();
      }//stc loop
      
    }
    
    if(outputType==TPGFEDataformat::STC16){
      
      const std::map<std::string,std::vector<uint32_t>>& modSTC16list = (pck.getDetType()==0)?configs.getSiModSTC16list():configs.getSciModSTC16list();
      const std::vector<uint32_t>& stc16list = modSTC16list.at(modName) ;
      const std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>& stc16TcMap = (pck.getDetType()==0)?configs.getSiSTC16ToTC():configs.getSciSTC16ToTC();
      uint16_t bx = 0xffff;
      emulOut.second.reset();
      for(const auto& istc16 : stc16list){
	bool isTcTp1 = false;
	bool isTcTp2 = false;
	bool isTcTp3 = false;
	if(istc16>=nofSTCs) continue;
	const std::vector<uint32_t>& tclist = stc16TcMap.at(std::make_pair(modName,istc16));
	const TPGFEDataformat::ModuleTcData& mdata = moddata.at(moduleId);
	bx = (mdata.getBx()==3564) ? 0xF : mdata.getBx() & 0x7 ; 
	uint64_t decompressedSTC16 = 0;
	std::vector<TPGFEDataformat::TcRawData> tcrawdatalist;
	TPGFEDataformat::TcRawData tcdata;
	for(const auto& econtc : tclist){
	  uint32_t hgctc = configs.getEconTPara().at(moduleId).getInputMux(econtc) ;
	  bool hasFound = false;
	  uint32_t emultc = 0xffffffff;
	  for (const auto& it : refMuxMap){
	    if (it.second == hgctc)  {
	      hasFound = true;
	      emultc = it.first;
	    }
	  }
	  if(!hasFound){
	    std::cerr << "TPGFEModuleEmulation::ECONTEmulation::EmulateJulSTC16 (moduleid="<<moduleId<<") : Mux not set for TC " << econtc << std::endl;
	    continue;
	  }
	  if(mdata.getTC(emultc).isTcTp1()) isTcTp1 = true;
	  if(mdata.getTC(emultc).isTcTp2()) isTcTp2 = true;
	  if(mdata.getTC(emultc).isTcTp3()) isTcTp3 = true;
	  uint64_t decompressed = DecompressEcont(mdata.getTC(emultc).getCdata(), configs.getEconTPara().at(moduleId).getDensity());
	  uint64_t decomp64bit = decompressed * configs.getEconTPara().at(moduleId).getCalibration(econtc) ;
	  if(isVerbose) std::cout << "TPGFEModuleEmulation::ECONTEmulation::EmulateJulSTC16 (moduleid="<<moduleId<<", TC="<<econtc<<") decompressed: " << decompressed << ", decompressed*calib: " << decomp64bit << std::endl; 
	  decompressed =  decomp64bit >> 11;
	  decompressedSTC16 += decompressed ;
	  if(isVerbose) std::cout << "TPGFEModuleEmulation::ECONTEmulation::EmulateJulSTC16 (moduleid="<<moduleId<<", TC="<<econtc<<") decompressed*calib>>11: " << decompressed << ", decompressedSTC16: " << decompressedSTC16 << std::endl; 
	  uint16_t compressed_bc = CompressEcontBc(decompressed,dropLSB);
	  tcdata.setTriggerCell(TPGFEDataformat::BestC, econtc, compressed_bc, decompressed, mdata.getTC(emultc).isTcTp1(), mdata.getTC(emultc).isTcTp2(), mdata.getTC(emultc).isTcTp3()) ;
	  tcrawdatalist.push_back(tcdata);
	}
	batcherOEMSort(tcrawdatalist);
	uint32_t lMaxId = findLastMax(tcrawdatalist);
	tcrawdatalist.resize(lMaxId+1);
	std::sort(tcrawdatalist.begin(),tcrawdatalist.end(),TPGFEDataformat::TcRawDataPacket::customGTA);
	TPGFEDataformat::TcRawData lastMax = isAllZero(tcrawdatalist) ? tcrawdatalist[tcrawdatalist.size()-1] : tcrawdatalist[0] ;
	uint16_t compressed_energy = CompressEcontStc5E4M(decompressedSTC16);
	emulOut.second.setTBM(outputType, bx, 0); 
	emulOut.second.setTcData(outputType, (lastMax.address())%16, compressed_energy, decompressedSTC16, isTcTp1, isTcTp2, isTcTp3);
	tcrawdatalist.clear();
      }//stc16 loop
    }//STC16 select condition
    
  }//Emulate Jul STC

  
  void ECONTEmulation::EmulateBC(bool isSim, uint64_t ievent, uint32_t& moduleId, const std::map<uint32_t,TPGFEDataformat::ModuleTcData>& moddata) {
    pck.setModId(moduleId);    
    const std::map<std::tuple<uint32_t,uint32_t,uint32_t>,std::string>& modNameMap = configs.getModIdxToName();
    const std::map<uint32_t,uint32_t>& refMuxMap = configs.getMuxMapping() ;
    const std::string& modName = modNameMap.at(std::make_tuple(pck.getDetType(),pck.getSelTC4(),pck.getModule()));
    uint32_t dropLSB = configs.getEconTPara().at(moduleId).getDropLSB() ;
    const std::map<std::string,std::vector<uint32_t>>& modTClist = (pck.getDetType()==0)?configs.getSiModTClist():configs.getSciModTClist();
    const std::vector<uint32_t>& tclist = modTClist.at(modName) ;
    
    const TPGFEDataformat::Type& outputType = configs.getEconTPara().at(moduleId).getOutType();
    const TPGFEDataformat::ModuleTcData& mdata = moddata.at(moduleId);

    uint16_t bx = (mdata.getBx()==3564) ? 0xF : mdata.getBx() & 0x7 ; //make 8 modulo
    uint64_t decompressedMS = 0;
    std::vector<TPGFEDataformat::TcRawData> tcrawdatalist;
    TPGFEDataformat::TcRawData tcdata;
    for(const auto& econtc : tclist){
      uint32_t hgctc = configs.getEconTPara().at(moduleId).getInputMux(econtc) ;
      bool hasFound = false;
      uint32_t emultc = 0xffff;
      for (const auto& it : refMuxMap){
	if (it.second == hgctc)  {
	  hasFound = true;
	  emultc = it.first;
	}
      }
      if(!hasFound){
	std::cerr << "TPGFEModuleEmulation::ECONTEmulation::EmulateBC (moduleid="<<moduleId<<") : Mux not set for TC " << econtc << std::endl;
	continue;
      }
      uint64_t decompressed = DecompressEcont(mdata.getTC(emultc).getCdata(), configs.getEconTPara().at(moduleId).getDensity()); //funny that return needs to defined as 64-bit
      uint64_t decomp64bit = decompressed * configs.getEconTPara().at(moduleId).getCalibration(econtc); //overflows for 12bit TOT
      if(isVerbose) std::cout << "TPGFEModuleEmulation::ECONTEmulation::EmulateBC (moduleid="<<moduleId<<", TC="<<econtc<<") decompressed: " << decompressed << ", decompressed*calib: " << decomp64bit << std::endl; 
      decompressed =  decomp64bit >> 11;
      decompressedMS += (decompressed >> dropLSB) ;
      if(isVerbose) std::cout << "TPGFEModuleEmulation::ECONTEmulation::EmulateBC (moduleid="<<moduleId<<", TC="<<econtc<<") decompressed*calib>>11: " << decompressed << ", decompressedMS: " << decompressedMS << std::endl; 
      uint16_t compressed_bc = CompressEcontBc(decompressed,dropLSB);
      tcdata.setTriggerCell(outputType, econtc, compressed_bc, decompressed>>dropLSB, mdata.getTC(emultc).isTcTp1(), mdata.getTC(emultc).isTcTp2(), mdata.getTC(emultc).isTcTp3()) ;
      tcrawdatalist.push_back(tcdata);
    }
    batcherOEMSort(tcrawdatalist);
    if(isVerbose) for(int itc=0;itc<48;itc++) tcrawdatalist[itc].print();
    
    uint16_t compressed_modsum = CompressEcontModsum(decompressedMS,dropLSB);
    emulOut.second.reset();
    emulOut.second.setTBM(outputType, bx, compressed_modsum, decompressedMS>>dropLSB); 
    uint32_t nofBCTcs = configs.getEconTPara().at(moduleId).getBCType();
    for(uint32_t itc = 0 ; itc<nofBCTcs ; itc++)
      emulOut.second.setTcData(outputType, tcrawdatalist[itc].address(), tcrawdatalist[itc].energy(), tcrawdatalist[itc].rawE(),tcrawdatalist[itc].isTcTp1(),tcrawdatalist[itc].isTcTp2(),tcrawdatalist[itc].isTcTp3());
    
  }//Emulate BC
  
  void ECONTEmulation::convertToElinkData(unsigned bx, const TPGFEDataformat::TcRawDataPacket &tcrdp, uint32_t *ve) {
    
    const std::vector<TPGFEDataformat::TcRawData>& vtc = tcrdp.getTcData();
    
    unsigned nVe(0);
    
    assert(vtc.size()>0);
    
    uint64_t data(uint64_t(bx)<<60);
    unsigned last(60);
    unsigned w(0);
    
    TPGFEDataformat::Type tcType(tcrdp.type());
    bool bcMap(false);
    uint64_t bcBits(0);
    
    if(tcType==TPGFEDataformat::BestC) {
      //assert(vtc[0].isModuleSum());
      data|=(uint64_t(tcrdp.moduleSum())<<52);
      last-=8;
      w++;
      
      bcMap=(vtc.size()>7); //High occupancy if nTC>=8
      
      ////////////// construct map for high occupancy ///////////////
      if(bcMap) {
        for(unsigned i(0);i<vtc.size();i++) {
          data|=uint64_t(1)<<(51-vtc[i].address());
        }
        last-=48;
      }
      ////////////////////////////////////////////
    }
    
    if(!bcMap) {
      for(w = 0;w<vtc.size() && !bcMap;w++) {
        if(last<32) {
          ve[nVe++]=data>>32;
          data=(data<<32);
          last+=32;
        }

        if     (tcType==TPGFEDataformat::BestC) last-=6;
        else if(tcType==TPGFEDataformat::STC4A) last-=2;
        else if(tcType==TPGFEDataformat::STC4B) last-=2;
        else if(tcType==TPGFEDataformat::CTC4A) last-=0;
        else if(tcType==TPGFEDataformat::CTC4B) last-=0;
        else if(tcType==TPGFEDataformat::STC16) last-=4;
        else assert(false);

	if(tcType!=TPGFEDataformat::CTC4A and tcType!=TPGFEDataformat::CTC4B)
	  data|=(uint64_t(vtc[w].address())<<last);
      }
    }
    
    //unsigned wLo(tcType==TPGFEDataformat::BestC?1:0);
    for(unsigned w(0);w<vtc.size();w++) {
      if(last<32) {
        ve[nVe++]=data>>32;
        data=(data<<32);
        last+=32;
      }

      if     (tcType==TPGFEDataformat::BestC) last-=7;
      else if(tcType==TPGFEDataformat::STC4A) last-=7;
      else if(tcType==TPGFEDataformat::STC4B) last-=9;
      else if(tcType==TPGFEDataformat::CTC4A) last-=7;
      else if(tcType==TPGFEDataformat::CTC4B) last-=9;
      else if(tcType==TPGFEDataformat::STC16) last-=9;
      else assert(false);

      data|=(uint64_t(vtc[w].energy())<<last);
    }
    
    // Catch last words
    if(last<32) {
      ve[nVe++]=data>>32;
      data=(data<<32);
      last+=32;
    }
    
    ve[nVe++]=data>>32;
  }
  
  void ECONTEmulation::getElinkData(unsigned bx, uint32_t *ptr) const {
    convertToElinkData(bx,emulOut.second,ptr);
  }

  void ECONTEmulation::fillRandomTcRawData(unsigned bx) {
    // FIXME - HOW TO GET TYPE AND NUMBER OF TCS?
    //generateRandomTcRawData(bx,configs.getEconTPara().at(moduleId).getOutType(),9,emulOut);
  }

  void ECONTEmulation::fillZeroEnergyTcRawData(unsigned bx) {
    // FIXME - HOW TO GET TYPE AND NUMBER OF TCS?
    //generateZeroEnergyTcRawData(bx,configs.getEconTPara().at(moduleId).getOutType(),9,emulOut);
  }

  void ECONTEmulation::generateRandomTcRawData(unsigned bx,
					       TPGFEDataformat::Type type,
					       unsigned nTc,
					       TPGFEDataformat::TcRawDataPacket &vtcrp) {
    generateTcRawData(false,bx,type,nTc,vtcrp);
  }
  
  void ECONTEmulation::generateZeroEnergyTcRawData(unsigned bx,
						   TPGFEDataformat::Type type,
						   unsigned nTc,
						   TPGFEDataformat::TcRawDataPacket &vtcrp) {

    generateTcRawData(true,bx,type,nTc,vtcrp);
  }
  
  void ECONTEmulation::generateTcRawData(bool zero, unsigned bx,
					 TPGFEDataformat::Type type,
					 unsigned nTc,
					 TPGFEDataformat::TcRawDataPacket &vtcrp) {
    
    //std::vector<TPGFEDataformat::TcRawData> &vtc(vtcrp.second);
    std::vector<TPGFEDataformat::TcRawData> &vtc(vtcrp.setTcData());
    vtcrp.setType(type);
  
    if(type==TPGFEDataformat::BestC) {
      vtc.resize(nTc);
      vtcrp.setModuleSum(zero?0:rand()&0xff);
    
      unsigned step(48/nTc);
      //for(unsigned i(1);i<vtc.size();i++) {
      for(unsigned i(0);i<vtc.size();i++) {
	      vtc[i].setTriggerCell(type,step*(i)+(rand()%step),zero?0:rand()&0x7f);
      }
    
    } else if(type==TPGFEDataformat::STC4A) {
      vtc.resize(nTc);
      for(unsigned i(0);i<vtc.size();i++) {
	vtc[i].setTriggerCell(type,rand()&0x3,zero?0:rand()&0x7f);
      }
    
    } else if(type==TPGFEDataformat::STC4B) {
      vtc.resize(nTc);
      for(unsigned i(0);i<vtc.size();i++) {
	vtc[i].setTriggerCell(type,rand()&0x3,zero?0:rand()&0x1ff);
      }
    
    } else if(type==TPGFEDataformat::STC16) {
      vtc.resize(nTc);
      for(unsigned i(0);i<vtc.size();i++) {
	vtc[i].setTriggerCell(type,rand()&0xf,zero?0:rand()&0x1ff);
      }
    
    } else {
      assert(false);
    }
  }

  
}//end of namespace


#endif
