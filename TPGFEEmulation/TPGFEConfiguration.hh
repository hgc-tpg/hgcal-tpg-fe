#ifndef TPGFEConfiguration_h
#define TPGFEConfiguration_h

#include <iostream>
#include <iomanip>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <map>

#include <yaml-cpp/yaml.h>

#include "TMath.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"

#include "TPGFEDataformat.hh"

namespace TPGFEConfiguration{
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////The configuration of half of ROC based on HGCROC3a [doc. no. v2.0] (See Table@Page-43)
  //////EDMS ROCv3a: https://edms.cern.ch/ui/#!master/navigator/document?D:100570166:100570166:subDocs
  //////EDMS ROCv3b(recent): https://edms.cern.ch/ui/#!master/navigator/document?D:101362066:101362066:subDocs
  class ConfigHfROC {    
  public:
    ConfigHfROC() {}
    uint32_t getAdcTH() const { return uint32_t(Adc_TH);}
    uint64_t getClrAdcTottrig() const { return ClrAdcTot_trig;}
    bool isChMasked(uint32_t ich) const {
      int chnl = ich%36;
      return (getClrAdcTottrig()>>chnl) & 0x1 ;
    }
    uint32_t getTotTH(uint32_t ich) const {
      uint32_t chnl = ich%36;
      uint32_t  tot_idx = TMath::FloorNint(chnl/9);
      return uint32_t(Tot_TH[tot_idx]);
    }
    uint32_t getTotP(uint32_t ich) const {
      uint32_t chnl = ich%36;
      uint32_t tot_idx = TMath::FloorNint(chnl/9);
      return uint32_t(Tot_P[tot_idx]);
    }
    uint32_t getMultFactor() const { return uint32_t(MultFactor);}    
    void setAdcTH(uint32_t  adcth) { Adc_TH = adcth & 0x1F;}
    void setClrAdcTottrig(uint64_t clradctottrig) { ClrAdcTot_trig = clradctottrig & 0xFFFFFFFF;}
    void setMultFactor(uint32_t multfactor) { MultFactor = multfactor & 0x1F;}
    void setTotTH(uint32_t tot_idx, uint32_t tot_th) { Tot_TH[tot_idx] = tot_th & 0xFF;}
    void setTotP(uint32_t tot_idx, uint32_t tot_p) { Tot_P[tot_idx] = tot_p & 0x7F;}
    void print() const {
      std::cout << std::dec << ::std::setfill(' ')
		<< "ConfigHfROC(" << this << ")::print(): "
		<<"Adc_TH = "<< std::setw(4) << getAdcTH()
		<<", MultFactor = "<< std::setw(3) << getMultFactor()
		<< std::endl;
      
      std::cout << std::dec << ::std::setfill(' ')
		<< "ConfigHfROC(" << this << ")::print(): "
		<<"ClrAdcTot_trig = ";
      for(uint32_t ich=0;ich<36;ich++)
	std::cout << std::setw(2) << "("<< ich <<": " << isChMasked(ich) <<") ";
      std::cout << std::endl;
      
      std::cout << std::dec << ::std::setfill(' ')
		<< "ConfigHfROC(" << this << ")::print(): "
		<<"Tot_P = ";
      for(uint32_t itotch=0;itotch<4;itotch++)
	std::cout << std::setw(4) << "("<< itotch <<": " << uint32_t(Tot_P[itotch]) <<") ";
      std::cout << std::endl;
      
      std::cout << std::dec << ::std::setfill(' ')
		<< "ConfigHfROC(" << this << ")::print(): "
		<<"Tot_TH = ";
      for(uint32_t itotch=0;itotch<4;itotch++)
	std::cout << std::setw(4) << "("<< itotch <<": " << uint32_t(Tot_TH[itotch]) <<") ";
      std::cout << std::endl;

    }

  private:    
    //Digital Info
    uint8_t Adc_TH; //5-bits
    uint64_t ClrAdcTot_trig;  //36-bits
    uint8_t MultFactor; //5-bits
    uint8_t Tot_P[4];  //one per 9 channel (each with 7-bits):not present in 2023 beam test
    uint8_t Tot_TH[4]; //one per 9 channel (each with 8 bits):not present in 2023 beam test    
  };
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //The configuration of channel corresponding to ADC per module
  class ConfigCh {
  public:
    ConfigCh() {}
    uint32_t getAdcpedestal() const { return uint32_t(Adc_pedestal);}
    void setAdcpedestal(uint32_t ped) { Adc_pedestal = ped & 0xFF;}
    void print() {
      std::cout << std::dec << ::std::setfill(' ')
		<< "ConfigCh(" << this << ")::print(): "
		<<"Adc_pedestal = "<< std::setw(4)<< getAdcpedestal()
		<< std::endl;
    }
    void print(uint32_t ich) {
      std::cout << std::dec << ::std::setfill(' ')
		<< "ConfigCh(" << this << ")::print(): "
		<<"ich: "<< ich <<", Adc_pedestal = "<< std::setw(4)<< getAdcpedestal()
		<< std::endl;
    }
    
  private:
    uint8_t Adc_pedestal; //8-bits 
  };
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////ECON-D [doc. no: v1.1]
  //////https://edms.cern.ch/ui/#!master/navigator/document?P:100053490:100904542:subDocs
  class ConfigEconD {
  public:
    ConfigEconD() : isPassThrough(false), neRx(0) {}
    bool passThrough() const { return isPassThrough;}
    uint32_t getNeRx() const { assert(neRx!=0); return uint32_t(neRx);}
    void setPassThrough(bool isPT) { isPassThrough = isPT;}
    void setNeRx(uint32_t nofeRx) { assert(nofeRx!=0); neRx = nofeRx;}
    void print() {
      std::cout << std::dec << ::std::setfill(' ')
		<< "ConfigEconD(" << this << ")::print(): "
		<<"isPassThrough mode = "
		<< std::setw(2) << passThrough()
		<<", \tnof eRx = "
		<< std::setw(2) << getNeRx()
		<< std::endl;
    }
    
  private:    
    bool isPassThrough; //1-bit
    uint8_t neRx; 
  };
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////ECON-T [doc. no: v10]
  //////https://edms.cern.ch/ui/#!master/navigator/document?P:100053490:100430098:subDocs
  class ConfigEconT {
  public:
    ConfigEconT() : density(0), dropLSB(0), select(0), stc_type(0) {
      for(uint32_t itc=0;itc<48;itc++) {
	calv[itc] = 0;
	mux[itc] = 0;
      }
    }
    bool getDensity() const { assert(density==0 or density==1); return (density==1) ? true : false;}
    uint32_t getDropLSB() const { assert(dropLSB>=0 and dropLSB<=4); return uint32_t(dropLSB);}
    uint32_t getSelect() const { assert(select==1 or select==2); return uint32_t(select);}
    uint32_t getSTCType() const { assert(stc_type>=0 and stc_type<=4); return uint32_t(stc_type);}
    uint32_t getNElinks() const { assert(eporttx_numen!=0);  return uint32_t(eporttx_numen);}
    uint32_t getNofTCs() const { return getBCType(); }
    uint32_t getBCType() const {
      uint32_t maxTcs = 0;
      if(getOutType()==TPGFEDataformat::BestC){	
	switch(getNElinks()){
	case 1:
	  maxTcs = 1;
	  break;
	case 2:
	  maxTcs = 4;
	  break;
	case 3:
	  maxTcs = 6;
	  break;
	case 4:
	  maxTcs = 9;
	  break;
	case 5:
	  maxTcs = 14;
	  break;
	case 6:
	  maxTcs = 18;
	  break;
	case 7:
	  maxTcs = 23;
	  break;
	case 8:
	  maxTcs = 28;
	  break;
	case 9:
	  maxTcs = 32;
	  break;
	case 10:
	  maxTcs = 37;
	  break;
	case 11:
	  maxTcs = 41;
	  break;
	case 12:
	  maxTcs = 46;
	  break;
	default: //same as case13 or above
	  maxTcs = 48;
	  break;
	}	  
      }
      return maxTcs;
    }
    uint32_t getNofSTCs() const {
      uint32_t maxSTCs = 0;
      switch(getSTCType()){
      case 0: //0 = STC4B(5E+4M)
	switch(getNElinks()){
	case 1:
	  maxSTCs = 2;
	  break;
	case 2:
	  maxSTCs = 5;
	  break;
	case 3:
	  maxSTCs = 8;
	  break;
	case 4:
	  maxSTCs = 11;
	  break;
	default: //case 5 or higher
	  maxSTCs = 12;
	  break;
	}
	break;
      case 1: //1 = STC16(5E+4M)
	switch(getNElinks()){
	case 1:
	  maxSTCs = 2;
	  break;
	default: //case 2 or higher
	  maxSTCs = 3;
	  break;
	}
	break;
      case 2: //2 = CTC4A(4E+3M)
	switch(getNElinks()){
	case 1:
	  maxSTCs = 4;
	  break;
	case 2:
	  maxSTCs = 8;
	  break;
	default: //case 3 or higher
	  maxSTCs = 12;
	  break;
	}
	break;
      default: //3 = STC4A(4E+3M) //4 = CTC4B(5E+3M)
	switch(getNElinks()){
	case 1:
	  maxSTCs = 3;
	  break;
	case 2:
	  maxSTCs = 6;
	  break;
	case 3:
	  maxSTCs = 10;
	  break;
	default: //case 4 or higher
	  maxSTCs = 12;
	  break; 
	}	
      }//stctype
      return maxSTCs;
    }
    uint32_t getCalibration(uint32_t itc) const {
      assert(itc<=47) ;
      return uint32_t(calv[itc]);
    }
    uint32_t getInputMux(uint32_t itc) const {
      assert(itc<=47) ;
      return (isConnectedMux(itc)) ? uint32_t(mux[itc]) : 0x80 ;
    }
    bool isConnectedMux(uint32_t itc) const {
      assert(itc<=47) ;
      return (mux[itc]>>7) ? false : true ;
    }
    TPGFEDataformat::Type getOutType() const {
      TPGFEDataformat::Type type;
      switch(select){
      case 0:
	type = TPGFEDataformat::TS;
	break;
      case 1:
	switch(stc_type){
	case 0:
	  type = TPGFEDataformat::STC4B;
	  break;
	case 1:
	  type = TPGFEDataformat::STC16;
	  break;
	case 2:
	  type = TPGFEDataformat::CTC4A;
	  break;
	case 3:
	  type = TPGFEDataformat::STC4A;
	  break;
	case 4:
	  type = TPGFEDataformat::CTC4B;
	  break;
	default:
	  type = TPGFEDataformat::Unknown;
	  break;
	}
	break;
      case 2:
	type = TPGFEDataformat::BestC;
	break;
      case 3:
	type = TPGFEDataformat::RA;
	break;
      case 4:
	type = TPGFEDataformat::AE;
	break;
      default:
	type = TPGFEDataformat::Unknown;
	break;
      }
      return type;
    }
    void setDensity(uint32_t den) { assert(density==0 or density==1); density = den;}
    void setDropLSB(uint32_t dLSB) { assert(dLSB>=0 and dLSB<=4); dropLSB = dLSB;}
    void setSelect(uint32_t sel) { assert(sel==1 or sel==2); select = sel;}
    void setSTCType(uint32_t stctype) { assert(stctype>=0 and stctype<=4); stc_type = stctype;}
    void setNElinks(uint32_t nlinks) { assert(nlinks!=0); eporttx_numen = nlinks;}
    void setCalibration(uint32_t itc, uint32_t calib) {
      assert(itc<=47) ;
      calv[itc] = (calib & 0xFFF);
    }
    void setInputMux(uint32_t itc, uint32_t muxval) {
      assert((muxval<=47 or muxval==0x80) and itc<=47) ; //since default value is not known, set to 0x80 for unconnected TC
      mux[itc] = muxval & 0xff;
    }
    void print() {
      std::cout << std::dec << ::std::setfill(' ')
		<<"ConfigEconT(" << this << ")::print(): "
		// <<", \tECON-T mode = "
		// << std::setw(8) << getOutType()
		<<", \tLSB used for TC input = "
		<< std::setw(2) << getDensity()
		<<", \tLSB in ECON-T output = "
		<< std::setw(2) << getDropLSB()
		<<", \tECON-T select = "
		<< std::setw(2) << getSelect()
		<<", \tECON-T stc_type = "
		<< std::setw(2) << getSTCType()
		<<", \tECON-T maxTCs = "
		<< std::setw(2) << ((getSelect()==2) ? getBCType() : getNofSTCs())
		<< std::endl;
      for(uint32_t itc=0;itc<48;itc++)
	std::cout <<"ConfigEconT(" << this << ")::print(): "
		  <<" itc: " << itc
		  <<", \tCALV = "
		  << std::setw(5) << getCalibration(itc)
		  <<", \tMUX = "
		  << std::setw(2) << getInputMux(itc)
		  << std::endl;
    }
    
  private:    
    uint8_t density; //lsb at the input TC from ROC
    uint8_t dropLSB; //lsb at the output during the packing
    uint8_t select; //0 = Threshold Sum (TS), 1 = Super Trigger Cell (STC), 2 = Best Choice (BC), 3 = Repeater, 4=Autoencoder (AE).
    uint8_t stc_type; //0 = STC4B(5E+4M), 1 = STC16(5E+4M), 2 = CTC4A(4E+3M), 3 = STC4A(4E+3M), 4 = CTC4B(5E+3M)
    uint8_t eporttx_numen;//number of elinks
    uint16_t calv[48]; //12-bit calibration for 48 TCs
    uint8_t mux[48];   //multiplexer between HGCROC and TC to ECONT
  };
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////Packing several IDs for listing objects
  class TPGFEIdPacking {
  public:
    TPGFEIdPacking() : packedVal(0) {}
    
    //June/2024 : Proposal1 : 1+2+11+1+3+1+3+3+1 = 1(+/-z) + 2(120deg sector) + 11(link max ~1848) + 1(Si/Sci) + 3(ECONT number max 7) + 1(LD/HD) + 4[(Si type)/(Sci type)] +  3(ROC number) + 1 (Half) = 27 bits
    //June/2024 : Proposal2 : 1+2+11+1+3+1+7+3+1 = 1(+/-z) + 2(120deg sector) + 11(link max ~1848) + 1(Si/Sci) + 3(ECONT number max 7) + 1(LD/HD) + [3(Si)/8(Sci)] +  3(ROC number) + 1 (Half) = 31 bits
    //June/2024 : Proposal3 : 1+2+6+6+1+3+1+3+3+1 = 1(+/-z) + 2(120deg sector) + 6(layers) + 6(max 40 per LD/HD/layer) + 1(Si/Sci) + 3(ECONT number max 7) + 1(LD/HD) + 4[(Si type)/(Sci type)] +  3(ROC number) + 1 (Half) = 28 bits
    //June/2024 : Proposal4 : 1+2+6+6+1+3+1+7+3+1 = 1(+/-z) + 2(120deg sector) + 6(layers) + 6(max 40 per LD/HD/layer) + 1(Si/Sci) + 3(ECONT number max 7) + 1(LD/HD) + [3(Si)/8(Sci)] +  3(ROC number) + 1 (Half) = 32 bits
    
    //Below we choose for proposal 1
    uint32_t getModId() const { return uint32_t(packedVal & 0x3FFFFF);}
    uint32_t getRocId() const { return uint32_t(packedVal & 0x3FFFFFF);}
    uint32_t getZside() const { return uint32_t((packedVal>>26) & 0x1);}
    uint32_t getSector() const { return uint32_t((packedVal>>24) & 0x3);}
    uint32_t getLink() const { return uint32_t((packedVal>>13) & 0x7FF);}
    uint32_t getDetType() const { return uint32_t((packedVal>>12) & 0x1);}
    uint32_t getEconN() const { return uint32_t((packedVal>>9) & 0x7);}
    uint32_t getSelTC4() const { return uint32_t((packedVal>>8) & 0x1);} 
    uint32_t getModule() const { return uint32_t((packedVal>>4) & 0xF);} 
    uint32_t getRocN() const { return uint32_t(packedVal>>1 & 0x7);}
    uint32_t getHalf() const { return uint32_t(packedVal & 0x1);}
    
    uint64_t getChId() const { return (packedVal & 0xFFFFFFFF);}
    uint32_t getRocChId(uint64_t rocchId) const { return uint32_t(rocchId & 0x3F);}
    
    uint32_t packModId(uint32_t zside, uint32_t sector, uint32_t link, uint32_t det, uint32_t econt, uint32_t selTC4, uint32_t module) {
      assert(zside<=0x1 and sector<=0x3 and link<=0x7FF and det<=0x1 and econt<=0x7 and selTC4<=0x1 and module<=0xF);
      //packedVal = (zside<<21) | (sector<<19) | (link<<8) | (det<<7) | (econt<<4) | (selTC4<<3) | module ;
      //This way modId is same as rocId of first ROC and channel 0, however it makes easier to use above methods to access the zside, link, sector,..... 
      uint32_t rocn = 0, half = 0;
      packedVal = (zside<<26) | (sector<<24) | (link<<13) | (det<<12) | (econt<<9) | (selTC4<<8) | (module<<4) |  (rocn<<1) | half ; 
      return packedVal;
    }
    uint32_t getRocIdFromModId(uint32_t moduleId, uint32_t rocn, uint32_t half){
      packedVal = moduleId | (rocn<<1 | half);
      return packedVal;
    }
    uint32_t getModIdFromRocId(uint32_t rocId){
      uint32_t rocn = 0, half = 0;
      uint32_t id = (rocId>>4) ;
      packedVal = (id<<4) | (rocn<<1 | half);
      return packedVal;
    }
    uint32_t getRocIdFromChId(uint64_t rocchId){
      packedVal = rocchId >> 6 ;
      return packedVal;
    }
    uint32_t getModIdFromChId(uint64_t rocchId){
      uint32_t rocn = 0, half = 0;
      uint32_t id = (rocchId>>10) ;
      packedVal = (id<<4) | (rocn<<1 | half);
      return packedVal;
    }
    uint32_t packRocId(uint32_t zside, uint32_t sector, uint32_t link, uint32_t det, uint32_t econt, uint32_t selTC4, uint32_t module, uint32_t rocn, uint32_t half) {
      assert(zside<=0x1 and sector<=0x3 and link<=0x7FF and det<=0x1 and econt<=0x7 and selTC4<=0x1 and module<=0xF and rocn<=0x6 and half<=0x1);
      packedVal = (zside<<26) | (sector<<24) | (link<<13) | (det<<12) | (econt<<9) | (selTC4<<8) | (module<<4) |  (rocn<<1) | half ; 
      return packedVal;
    }
    void printRocId(){
      std::cout << std::dec << ::std::setfill(' ') << std::setw(2)
		<<"  RocId : "<< getRocId()
		<< "\t (zside: " << getZside() <<", sector: "<< getSector() <<", link: "<< getLink() <<", isSi: "<< getDetType()
		<< ", EconN: "<< getEconN() << ", isLD : "<< getSelTC4() << ", module type : "<< getModule() << ", roc : " << getRocN() << ", half: "<< getHalf() << ")"
		<< std::endl;
    }
    uint64_t packChId(uint32_t rocid, uint32_t ch) {
      assert(ch<=0x3F);
      packedVal = (rocid<<6) | ch ;
      return packedVal;
    }

    void setModId(uint32_t id) { packedVal = id;}    
    void setRocId(uint32_t id) { packedVal = id;}    
    void setChId(uint64_t id) { packedVal = id ;}
    void setZero() {packedVal = 0;}
    uint32_t packij(uint32_t i, uint32_t j) {
      return (i<<6)|j  ;
    }
    std::pair<uint32_t,uint32_t> unpackij(uint32_t ij) {
      uint32_t j = ij & 0x3F ;
      uint32_t i = (ij>>6) & 0x3F ;
      return std::make_pair(i,j); 
    }
    
  private:
    uint64_t packedVal; //could be ROCid/chid/compressed TC energy....any generic packed value
    ////////////////////
    //The pin configuration could be different depending on the hexaboard/tileboard type, we need to account that.
    //////// Si-- ////////////
    //There are 5/4 hexaboard types for LD/HD as per "Si Modules" section of https://edms.cern.ch/ui/#!master/navigator/document?D:101059405:101059405:subDocs
    //1 bit for LD/HD + 3 bits hexboard types
    /////////////////////////
    //////// Sci-- ////////////
    //https://edms.cern.ch/ui/#!master/navigator/document?D:101205340:101205340:subDocs
    //https://edms.cern.ch/ui/#!master/navigator/document?P:100053490:101447883:subDocs
    //1 bit for LD/HD + 3 (8 categories A,B,C,D,E,G,J,K) +4(0-15 nof rings)) = 8 bits for scintillators
    
    //Note that the link number and HexBoard/TileBoard types are not in the pedestal yaml
  };
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  class Configuration{
  public:
    Configuration() {initId(); setTrainEWIndices(0,'w',0);}
    
    //setters
    void setSiChMapFile(std::string mapfile) {SiMapfname = mapfile;}
    void setSciChMapFile(std::string mapfile) {SciMapfname = mapfile;}
    void setRocFile(std::string pedThfile) {PedThfname = pedThfile;}
    void setEconDFile(std::string econdfile) {EconDfname = econdfile;}
    void setEconTFile(std::string econtfile) {EconTfname = econtfile;}
    void setRocPara(const std::map<uint32_t,TPGFEConfiguration::ConfigHfROC>& cfghroc){
      for(auto const& hroc : cfghroc)
	hroccfg[hroc.first] = cfghroc.at(hroc.first);
    }
    void setChPara(const std::map<uint64_t,TPGFEConfiguration::ConfigCh>& cfghrocch) {
      for(auto const& hrocch : cfghrocch)
	hrocchcfg[hrocch.first] = cfghrocch.at(hrocch.first);
    }
    void setEconDPara(const std::map<uint32_t,TPGFEConfiguration::ConfigEconD>& cfgeconD) {
      for(auto const& econD : cfgeconD)
	econDcfg[econD.first] = cfgeconD.at(econD.first);
    }
    void setEconTPara(const std::map<uint32_t,TPGFEConfiguration::ConfigEconT>& cfgeconT) {
      for(auto const& econT : cfgeconT)
	econTcfg[econT.first] = cfgeconT.at(econT.first);
    }
    
    void setTrainEWIndices(uint32_t tr_index, char ew_c, uint32_t ew_index) { train_idx = tr_index ; ew = ew_c ; ew_idx = ew_index ;}
    //set the module related parameters
    void setModulePath(uint32_t zs, uint32_t sect, uint32_t lnk, uint32_t SiorSci, uint32_t econt_index, uint32_t LDorHD, uint32_t mod_index){
      zside = zs; sector = sect; link = lnk; det = SiorSci;
      econt = econt_index; selTC4 = LDorHD; module = mod_index;
    }
    void initId(){
      zside = 0, sector = 0, link = 0, det = 0;
      econt = 0, selTC4 = 1, module = 0, rocn = 0, half = 0;
    }
    
    //Read the channel <--> pin mapping from text file
    void readSiChMapping();
    void readSciChMapping();
    void loadModIdxToNameMapping();
    void loadMuxMapping();
    void loadMuxMappingJul24(); //Hack, but needed to not break Jul24 testbeam script
    
    //Read the ROC/ECOND/ECONT configs from yaml file
    void readRocConfigYaml(const std::string& moduletype);
    void readRocConfigYaml(const uint32_t& rocid_0, const uint32_t& rocid_1);
    void readEconDConfigYaml();
    void readEconTConfigYaml();
    void readEconDConfigYaml(uint32_t idx);
    void readEconTConfigYaml(uint32_t idx);
    
    void setPedThZero(); //set the pedestal and thresholds to zero
    void setPedZero(); //set only the ped values to zero keep the other roc config values as loaded from config
    
    //getters
    const uint32_t getSiModNhroc(std::string typecode) { return  SiModNhroc[typecode];}
    const uint32_t getSciModNhroc(std::string typecode) { return  SciModNhroc[typecode];}
    const std::map<std::string,std::vector<uint32_t>>& getSiModTClist() { return  SiModTClst;}
    const std::map<std::string,std::vector<uint32_t>>& getSiModSTClist() { return  SiModSTClst;}
    const std::map<std::string,std::vector<uint32_t>>& getSiModSTC16list() { return  SiModSTC16lst;}
    const std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>& getSiTCToROCpin() {return  SiTCToROCpin;}
    const std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>& getSiSTCToTC() {return  SiSTCToTC;}
    const std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>& getSiSTC16ToTC() {return  SiSTC16ToTC;}
    const std::map<std::pair<std::string,std::tuple<uint32_t,uint32_t,uint32_t>>,uint32_t>& getSiSeqToROCpin() {return SiSeqToRocpin;}
    const std::map<std::pair<std::string,uint32_t>,uint32_t>& getSiRocpinToAbsSeq() {return SiRocpinToAbsSeq;}
    
    const std::map<std::string,std::vector<uint32_t>>& getSciModTClist() { return  SciModTClst;}
    const std::map<std::string,std::vector<uint32_t>>& getSciModSTClist() { return  SciModSTClst;}
    const std::map<std::string,std::vector<uint32_t>>& getSciModSTC16list() { return  SciModSTC16lst;}
    const std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>& getSciTCToROCpin() {return  SciTCToROCpin;}
    const std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>& getSciSTCToTC() {return  SciSTCToTC;}
    const std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>& getSciSTC16ToTC() {return  SciSTC16ToTC;}
    const std::map<std::pair<std::string,std::tuple<uint32_t,uint32_t,uint32_t>>,uint32_t>& getSciSeqToROCpin() {return SciSeqToRocpin;}
    const std::map<std::pair<std::string,uint32_t>,uint32_t>& getSciRocpinToAbsSeq() {return SciRocpinToAbsSeq;}
    
    const std::map<std::tuple<uint32_t,uint32_t,uint32_t>,std::string>& getModIdxToName() {return modIdxToName;}
    const std::map<uint32_t,uint32_t>& getMuxMapping() {return refMuxMap;}
    
    std::map<uint32_t,TPGFEConfiguration::ConfigHfROC>& getRocPara() { return hroccfg;}
    std::map<uint64_t,TPGFEConfiguration::ConfigCh>& getChPara() { return hrocchcfg;}
    std::map<uint32_t,TPGFEConfiguration::ConfigEconD>& getEconDPara() { return econDcfg;}
    std::map<uint32_t,TPGFEConfiguration::ConfigEconT>& getEconTPara() { return econTcfg;}
    ////////////////////////////////////////
    void printCfgPedTh(uint32_t moduleId);
  private:
    ////////////////////////////////////////
    //channel <--> pin mapping and related variables
    ////////////////////////////////////////
    std::string SiMapfname ;
    std::map<std::string,std::vector<uint32_t>>  SiModTClst ; std::map<std::string,std::vector<uint32_t>>  SiModSTClst; std::map<std::string,std::vector<uint32_t>>  SiModSTC16lst;
    std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>  SiTCToROCpin; std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>  SiTCToIJ;
    //std::map<std::pair<std::string,uint32_t>,uint32_t>  SiROCpinToTC; std::map<std::pair<std::string,uint32_t>,uint32_t>  SiIJToTC;
    //std::map<std::pair<std::string,uint32_t>,uint32_t>  SiTCToSTC;
    std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>  SiSTCToTC; std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>  SiSTC16ToTC;
    std::map<std::string,uint32_t> SiModNhroc;
    std::map<std::pair<std::string,std::tuple<uint32_t,uint32_t,uint32_t>>,uint32_t>  SiSeqToRocpin;
    std::map<std::pair<std::string,uint32_t>,uint32_t>  SiRocpinToAbsSeq;
    
    std::string SciMapfname ;
    std::map<std::string,std::vector<uint32_t>>  SciModTClst ; std::map<std::string,std::vector<uint32_t>>  SciModSTClst; std::map<std::string,std::vector<uint32_t>>  SciModSTC16lst;
    std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>  SciTCToROCpin; std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>  SciTCToIJ;
    //std::map<std::pair<std::string,uint32_t>,uint32_t>  SciROCpinToTC; std::map<std::pair<std::string,uint32_t>,uint32_t>  SciIJToTC;
    //std::map<std::pair<std::string,uint32_t>,uint32_t>  SciTCToSTC;
    std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>  SciSTCToTC; std::map<std::pair<std::string,uint32_t>,std::vector<uint32_t>>  SciSTC16ToTC;
    std::map<std::string,uint32_t> SciModNhroc;
    std::map<std::pair<std::string,std::tuple<uint32_t,uint32_t,uint32_t>>,uint32_t>  SciSeqToRocpin;
    std::map<std::pair<std::string,uint32_t>,uint32_t>  SciRocpinToAbsSeq;
    
    std::map<std::tuple<uint32_t,uint32_t,uint32_t>,std::string>  modIdxToName; //detType, LD/HD, modindex
    ////////////////////////////////////////
    
    ////////////////////////////////////////
    //id definition
    ////////////////////////////////////////
    uint32_t zside, sector, link, det;
    uint32_t econt, selTC4, module, rocn, half;
    TPGFEConfiguration::TPGFEIdPacking pck;
    ////////////////////////////////////////

    ////////////////////////////////////////
    //roc configuration
    ////////////////////////////////////////
    std::string PedThfname;
    uint32_t train_idx;
    char ew;
    uint32_t ew_idx;    
    std::map<uint32_t,TPGFEConfiguration::ConfigHfROC> hroccfg; std::map<uint64_t,TPGFEConfiguration::ConfigCh> hrocchcfg;
    ////////////////////////////////////////

    ////////////////////////////////////////
    //econd configuration
    ////////////////////////////////////////
    std::string EconDfname;
    std::map<uint32_t,TPGFEConfiguration::ConfigEconD> econDcfg;
    ////////////////////////////////////////

    ////////////////////////////////////////
    //econt configuration
    ////////////////////////////////////////
    std::string EconTfname;
    std::map<uint32_t,TPGFEConfiguration::ConfigEconT> econTcfg;
    std::map<uint32_t,uint32_t> refMuxMap;
    ////////////////////////////////////////
    
  };
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Configuration::readSiChMapping(){
  
    std::string Typecode;
    uint32_t ROC, HalfROC, Seq;
    std::string ROCpin;
    uint32_t ROCCH;
    int SiCell, TrLink, TrCell, iu, iv;
    float trace;
    int t;
    //Typecode ROC HalfROC Seq ROCpin SiCell TrLink TrCell iu iv trace t
    std::ifstream inwafermap(SiMapfname);
    std::stringstream ss;
    std::string s;
    
    int isLD = -1;
    TPGFEConfiguration::TPGFEIdPacking pck;  

    uint32_t prevHalfROC = 999;
    std::string prevTypecode = "";
    while(getline(inwafermap,s)){
      //std::cout << s.size() << std::endl;
      if(s.find("ML")!=std::string::npos or s.find("MH")!=std::string::npos){
	//std::cout << s << std::endl;
	ss << s.data() << std::endl;
	ss >> Typecode >> ROC >> HalfROC >> Seq >> ROCpin >> SiCell  >> TrLink  >> TrCell  >> iu  >> iv  >> trace  >> t ;
	if(ROCpin.find("CALIB")==std::string::npos and TrLink!=-1 and TrCell!=-1){
	  //std::cout << s ;//<< std::endl;
	  ROCCH = stoi(ROCpin);
	  isLD = (Typecode.find("MH-")==std::string::npos)?1:0;
	  uint32_t absTC = (isLD==1) ? (ROC*16 + TrLink*4 + TrCell) : (ROC*8 + TrLink*2 + TrCell) ;
	  uint32_t absSTC = uint32_t(TMath::FloorNint(absTC/4));
	  uint32_t absSTC16 = uint32_t(TMath::FloorNint(absTC/16));
	  uint32_t rocpin = ROC*72 + ROCCH;
	  uint32_t iUiV = pck.packij(uint32_t(iu),uint32_t(iv)) ;
	  std::tuple<uint32_t,uint32_t,uint32_t> seqch = std::make_tuple( ROC, HalfROC, Seq);
	  uint32_t absSeq = (2*ROC+HalfROC)*37 + Seq ;
	  
	  SiTCToROCpin[std::make_pair(Typecode,absTC)].push_back( rocpin );
	  SiTCToIJ[std::make_pair(Typecode,absTC)].push_back( iUiV );
	  // SiROCpinToTC[std::make_pair(Typecode,rocpin)] = absTC;
	  // SiIJToTC[std::make_pair(Typecode,iUiV)] = absTC;
	  SiSeqToRocpin[std::make_pair(Typecode,seqch)] = rocpin;
	  SiRocpinToAbsSeq[std::make_pair(Typecode,rocpin)] = absSeq;
	  
	  if(prevHalfROC!=HalfROC or prevTypecode.compare(Typecode)!=0){
	    if(prevTypecode.compare(Typecode)!=0){
	      SiModNhroc[Typecode] = 1;
	    }else{
	      SiModNhroc[Typecode]++;
	    }
	    prevHalfROC = HalfROC;
	    prevTypecode = Typecode ; 
	  }
	  if (std::find(SiModTClst[Typecode].begin(), SiModTClst[Typecode].end(), absTC) == SiModTClst[Typecode].end()) {
	    SiModTClst[Typecode].push_back(absTC);
	  }
	  if (std::find(SiModSTClst[Typecode].begin(), SiModSTClst[Typecode].end(), absSTC) == SiModSTClst[Typecode].end()) {
	    SiModSTClst[Typecode].push_back(absSTC);
	  }
	  if (std::find(SiModSTC16lst[Typecode].begin(), SiModSTC16lst[Typecode].end(), absSTC16) == SiModSTC16lst[Typecode].end()) {
	    SiModSTC16lst[Typecode].push_back(absSTC16);
	  }
	  // SiTCToSTC[std::make_pair(Typecode,absTC)] = absSTC;
	  if (std::find(SiSTCToTC[std::make_pair(Typecode,absSTC)].begin(), SiSTCToTC[std::make_pair(Typecode,absSTC)].end(), absTC) == SiSTCToTC[std::make_pair(Typecode,absSTC)].end()) {
	    SiSTCToTC[std::make_pair(Typecode,absSTC)].push_back( absTC );
	  }
	  if (std::find(SiSTC16ToTC[std::make_pair(Typecode,absSTC16)].begin(), SiSTC16ToTC[std::make_pair(Typecode,absSTC16)].end(), absTC) == SiSTC16ToTC[std::make_pair(Typecode,absSTC16)].end()) {
	    SiSTC16ToTC[std::make_pair(Typecode,absSTC16)].push_back( absTC );
	  }

	  //std::cout <<"\t"<< absSTC << "\t" << absTC  <<"\t"<< rocpin <<"\t"<< iu <<"\t"<< iv << std::endl;
	}//skip unconnected or calibration cells
      }//check type 
    }//loop over files
  
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Configuration::readSciChMapping(){
    
    std::string Typecode;
    uint32_t ROC, HalfROC;
    int Seq, ROCpin;
    int TrLink, TrCell;
    std::string iring;
    int iphi;
    float t, trace;

    std::string debugmodule = "TM-B11B12" ;
    bool isDebug = false;
    
    //Typecode ROC HalfROC Seq ROCpin TrLink TrCell iring iphi trace t
    std::ifstream inwafermap(SciMapfname);
    std::stringstream ss;
    std::string s;
    
    int isLD = -1;
    uint32_t ring = 0;
    TPGFEConfiguration::TPGFEIdPacking pck;
    
    uint32_t prevHalfROC = 999;
    std::string prevTypecode = "";
    while(getline(inwafermap,s)){
      //std::cout << s.size() << std::endl;
      if(s.find("TM-")!=std::string::npos){
	//std::cout << s << std::endl;
	ss << s.data() << std::endl;
	ss >> Typecode >> ROC >> HalfROC >> Seq >> ROCpin >> TrLink  >> TrCell  >> iring  >> iphi  >> t >> trace ;
	if(TrLink!=-1 and TrCell!=-1 and ROCpin>=0 and Seq>=0 and iring.find("-1")==std::string::npos and iphi!=-1 and t!=-1){
	  
	  if(s.find(debugmodule)!=std::string::npos and isDebug) std::cout << s << std::endl;
	  
	  isLD = (iring.find("h")==std::string::npos)?1:0;
	  if(iring.find("h")==std::string::npos)
	    ring = stoi(iring);
	  else
	    ring = stoi(iring.substr(1,iring.size()));
	  TrLink = (isLD==1)? ((ROC==1)?(TrLink-5):(TrLink-1)) : ((ROC==1)?(TrLink-3):(TrLink-1));
	  //uint32_t absTC = (isLD==1 ) ? (ROC*16 + TrLink*4 + TrCell) : (ROC*8 + TrLink*2 + TrCell) ;
	  uint32_t absTC = (isLD==1 ) ? (ROC*16 + TrLink*4 + TrCell) : (ROC*8 + HalfROC*4 + TrCell) ;
	  uint32_t absSTC = uint32_t(TMath::FloorNint(absTC/4));
	  uint32_t absSTC16 = uint32_t(TMath::FloorNint(absTC/16));
	  if(s.find(debugmodule)!=std::string::npos and isDebug) std::cout << "TPGFEConfiguration::Configuration::readSciChMapping: isLD: " << isLD << ", ROC: " << ROC <<", TrLink: " << TrLink <<", TrCell: " << TrCell << ", absTC: " << absTC  << ", absSTC: " << absSTC << ", absSTC16: " << absSTC16 << std::endl;
	  uint32_t rocpin = ROC*72 + uint32_t(ROCpin);
	  uint32_t iRiP = pck.packij(ring,uint32_t(iphi));
	  std::tuple<uint32_t,uint32_t,uint32_t> seqch = std::make_tuple( ROC, HalfROC, uint32_t(Seq));
	  uint32_t absSeq = (2*ROC+HalfROC)*37 + Seq ;
	  
	  if(prevHalfROC!=HalfROC or prevTypecode.compare(Typecode)!=0){
	    if(prevTypecode.compare(Typecode)!=0){
	      SciModNhroc[Typecode] = 1;
	    }else{
	      SciModNhroc[Typecode]++;
	    }
	    prevHalfROC = HalfROC;
	    prevTypecode = Typecode ; 
	  }
	  SciTCToROCpin[std::make_pair(Typecode,absTC)].push_back( rocpin );
	  SciTCToIJ[std::make_pair(Typecode,absTC)].push_back( iRiP );
	  // SciROCpinToTC[std::make_pair(Typecode,rocpin)] = absTC;
	  // SciIJToTC[std::make_pair(Typecode,rocpin)] = absTC;
	  SiSeqToRocpin[std::make_pair(Typecode,seqch)] = rocpin;
	  SciRocpinToAbsSeq[std::make_pair(Typecode,rocpin)] = absSeq;
	  
	  if (std::find(SciModTClst[Typecode].begin(), SciModTClst[Typecode].end(), absTC) == SciModTClst[Typecode].end()) {
	    SciModTClst[Typecode].push_back(absTC);
	  }
	  if (std::find(SciModSTClst[Typecode].begin(), SciModSTClst[Typecode].end(), absSTC) == SciModSTClst[Typecode].end()) {
	    SciModSTClst[Typecode].push_back(absSTC);
	  }
	  if (std::find(SciModSTC16lst[Typecode].begin(), SciModSTC16lst[Typecode].end(), absSTC16) == SciModSTC16lst[Typecode].end()) {
	    SciModSTC16lst[Typecode].push_back(absSTC16);
	  }

	  if (std::find(SciSTCToTC[std::make_pair(Typecode,absSTC)].begin(), SciSTCToTC[std::make_pair(Typecode,absSTC)].end(), absTC) == SciSTCToTC[std::make_pair(Typecode,absSTC)].end()) {
	    SciSTCToTC[std::make_pair(Typecode,absSTC)].push_back( absTC );
	  }

	  if (std::find(SciSTC16ToTC[std::make_pair(Typecode,absSTC16)].begin(), SciSTC16ToTC[std::make_pair(Typecode,absSTC16)].end(), absTC) == SciSTC16ToTC[std::make_pair(Typecode,absSTC16)].end()) {
	    SciSTC16ToTC[std::make_pair(Typecode,absSTC16)].push_back( absTC );
	  }
	  //std::cout <<"\t"<< absSTC << "\t" << absTC  <<"\t"<< rocpin <<"\t"<< ring <<"\t"<< iphi << std::endl;
	}//skip the unconnected channels
      }//Read Scintillator modules
    }//loop over full file 
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Configuration::loadModIdxToNameMapping(){
    //June/2024 : Proposal1 : 1+2+11+1+3+1+3+3+1 = 1(+/-z) + 2(120deg sector) + 11(link max ~1848) + 1(Si/Sci) + 3(ECONT number max 7) + 1(LD/HD) + 4[(Si type)/(Sci type)] +  3(ROC number) + 1 (Half) = 26 bits
    //detType (Si=0,Sci=1), LD/HD (LD=1,HD=0), modindex
    //(LD: Full (F/0), Top (T/1), Bottom (B/2), Left (L/3), Right (R/4), Five (5/5)=======HD: Full (F/0), Top (T/1) aka Half Minus, Bottom (B/2) aka Chop Two, Left (L/3), Right(R/4))
    
    /////////    //To be confirmed from geometry group modnum<-->modtype
    
    modIdxToName[std::make_tuple(0,1,0)] = "ML-F";
    modIdxToName[std::make_tuple(0,1,1)] = "ML-T";
    modIdxToName[std::make_tuple(0,1,2)] = "ML-B";
    modIdxToName[std::make_tuple(0,1,3)] = "ML-L";
    modIdxToName[std::make_tuple(0,1,4)] = "ML-R";
    modIdxToName[std::make_tuple(0,1,5)] = "ML-5";
    
    modIdxToName[std::make_tuple(0,0,0)] = "MH-F";
    modIdxToName[std::make_tuple(0,0,1)] = "MH-T";
    modIdxToName[std::make_tuple(0,0,2)] = "MH-B";
    modIdxToName[std::make_tuple(0,0,3)] = "MH-L";
    modIdxToName[std::make_tuple(0,0,4)] = "MH-R";

    modIdxToName[std::make_tuple(1,1,0)] = "TM-A5A6";
    modIdxToName[std::make_tuple(1,1,1)] = "TM-B11B12";
    modIdxToName[std::make_tuple(1,1,2)] = "TM-C5";
    modIdxToName[std::make_tuple(1,1,3)] = "TM-D8";
    modIdxToName[std::make_tuple(1,1,4)] = "TM-E8";
    modIdxToName[std::make_tuple(1,1,5)] = "TM-G3";
    modIdxToName[std::make_tuple(1,1,6)] = "TM-G5";
    modIdxToName[std::make_tuple(1,1,7)] = "TM-G7";
    modIdxToName[std::make_tuple(1,1,8)] = "TM-G8";

    modIdxToName[std::make_tuple(1,0,0)] = "TM-K6";
    modIdxToName[std::make_tuple(1,0,1)] = "TM-K8";
    modIdxToName[std::make_tuple(1,0,2)] = "TM-K11";
    modIdxToName[std::make_tuple(1,0,3)] = "TM-K12";
    modIdxToName[std::make_tuple(1,0,4)] = "TM-J12";
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Configuration::loadMuxMapping(){
    //uint32_t type_F_muxmap[48] = {3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8, 15, 14, 13, 12, 19, 18, 17, 16, 23, 22, 21, 20, 27, 26, 25, 24, 31, 30, 29, 28, 35, 34, 33, 32, 39, 38, 37, 36, 43, 42, 41, 40, 47, 46, 45, 44};
    uint32_t type_F_muxmap[48] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47};
    for(uint32_t itc=0;itc<48;itc++) refMuxMap[itc] = type_F_muxmap[itc];
      
  }
  void Configuration::loadMuxMappingJul24(){
    uint32_t type_F_muxmap[48] = {3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8, 15, 14, 13, 12, 19, 18, 17, 16, 23, 22, 21, 20, 27, 26, 25, 24, 31, 30, 29, 28, 35, 34, 33, 32, 39, 38, 37, 36, 43, 42, 41, 40, 47, 46, 45, 44};
    for(uint32_t itc=0;itc<48;itc++) refMuxMap[itc] = type_F_muxmap[itc];
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Configuration::readRocConfigYaml(const std::string& modName)
  {
    std::cout<<"Configuration::readRocConfigYaml:: Module config file : "<<PedThfname<<std::endl;
    YAML::Node node(YAML::LoadFile(PedThfname));

    const uint32_t nhfrocs = (pck.getDetType()==0)?getSiModNhroc(modName):getSciModNhroc(modName);
    const uint32_t nrocs = TMath::CeilNint(nhfrocs/2);
    const uint32_t nchs = 2*TPGFEDataformat::HalfHgcrocData::NumberOfChannels;
    for(uint32_t iroc=0;iroc<nrocs;iroc++){
      std::string EW(1,ew);
      std::string rocname = "train_" + std::to_string(train_idx) + ".roc" + std::to_string(iroc) + "_" + EW + std::to_string(ew_idx);
      uint32_t th_0 = node["Configuration"]["roc"][rocname]["cfg"]["DigitalHalf"]["0"]["Adc_TH"].as<uint32_t>();
      uint32_t th_1 = node["Configuration"]["roc"][rocname]["cfg"]["DigitalHalf"]["1"]["Adc_TH"].as<uint32_t>();
      uint64_t chmask_0 = node["Configuration"]["roc"][rocname]["cfg"]["DigitalHalf"]["0"]["ClrAdcTot_trig"].as<uint64_t>();
      uint64_t chmask_1 = node["Configuration"]["roc"][rocname]["cfg"]["DigitalHalf"]["1"]["ClrAdcTot_trig"].as<uint64_t>();
      uint32_t seltc4_0 = node["Configuration"]["roc"][rocname]["cfg"]["DigitalHalf"]["0"]["SelTC4"].as<uint32_t>();
      uint32_t seltc4_1 = node["Configuration"]["roc"][rocname]["cfg"]["DigitalHalf"]["1"]["SelTC4"].as<uint32_t>();
      uint32_t multfactor_0 = node["Configuration"]["roc"][rocname]["cfg"]["DigitalHalf"]["0"]["MultFactor"].as<uint32_t>();
      uint32_t multfactor_1 = node["Configuration"]["roc"][rocname]["cfg"]["DigitalHalf"]["1"]["MultFactor"].as<uint32_t>();
      //std::cout<<"rocname : "<<rocname<<", th_0 : "<<th_0<<", th_1 : "<<th_1<<", chmask_0 : "<<chmask_0<<", chmask_1 : "<<chmask_1<<std::endl;
      rocn = iroc;
      selTC4 = seltc4_0; half = 0;
      uint32_t rocid_0 = pck.packRocId(zside, sector, link, det, econt, selTC4, module, rocn, half);
      selTC4 = seltc4_1; half = 1;
      uint32_t rocid_1 = pck.packRocId(zside, sector, link, det, econt, selTC4, module, rocn, half);

      uint32_t dummy_tot_th[4], dummy_tot_p[4];
      for(int itot=0;itot<4;itot++){
	dummy_tot_th[itot] = 0 ;
	dummy_tot_p[itot] = 0 ;
      }
      
      TPGFEConfiguration::ConfigHfROC hroc_0, hroc_1;
      
      hroc_0.setAdcTH(th_0);
      hroc_0.setClrAdcTottrig(chmask_0);
      hroc_0.setMultFactor(multfactor_0);
      for(int itot=0;itot<4;itot++){
	hroc_0.setTotTH(itot, dummy_tot_th[itot]);
	hroc_0.setTotP(itot, dummy_tot_p[itot]);
      }

      hroc_1.setAdcTH(th_1);
      hroc_1.setClrAdcTottrig(chmask_1);
      hroc_1.setMultFactor(multfactor_1);
      for(int itot=0;itot<4;itot++){
	hroc_1.setTotTH(itot, dummy_tot_th[itot]);
	hroc_1.setTotP(itot, dummy_tot_p[itot]);
      }
    
      hroccfg[rocid_0] = hroc_0;
      hroccfg[rocid_1] = hroc_1;
    
      for(uint32_t ich=0;ich<nchs;ich++){
	uint32_t ihalf = (ich<TPGFEDataformat::HalfHgcrocData::NumberOfChannels)?0:1;
	uint32_t chnl = ich%TPGFEDataformat::HalfHgcrocData::NumberOfChannels;
	uint32_t ped = node["Configuration"]["roc"][rocname]["cfg"]["ch"][std::to_string(ich)]["Adc_pedestal"].as<uint32_t>() ;
	if(ihalf==0){
	  TPGFEConfiguration::ConfigCh ch_0;
	  ch_0.setAdcpedestal(ped);
	  hrocchcfg[pck.packChId(rocid_0,chnl)] = ch_0;
	}else{
	  TPGFEConfiguration::ConfigCh ch_1;
	  ch_1.setAdcpedestal(ped);
	  hrocchcfg[pck.packChId(rocid_1,chnl)] = ch_1;
	}
      }//channel loop
    }//roc loop
    //std::cout<<"============"<<std::endl;  
  }//end of read ped class
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Configuration::readRocConfigYaml(const uint32_t& rocid_0, const uint32_t& rocid_1)
  {
    std::cout<<"Configuration::readRocConfigYaml:: Module config file : "<<PedThfname<<std::endl;
    YAML::Node node(YAML::LoadFile(PedThfname));

    std::string totname;
    uint32_t tot_th_0[4], tot_p_0[4], tot_th_1[4], tot_p_1[4];
    const uint32_t nchs = 2*TPGFEDataformat::HalfHgcrocData::NumberOfChannels;
    uint32_t th_0 = node["DigitalHalf"]["0"]["Adc_TH"].as<uint32_t>();
    uint32_t th_1 = node["DigitalHalf"]["1"]["Adc_TH"].as<uint32_t>();
    // uint64_t chmask_0 = node["DigitalHalf"]["0"]["ClrAdcTot_trig"].as<uint64_t>();
    // uint64_t chmask_1 = node["DigitalHalf"]["1"]["ClrAdcTot_trig"].as<uint64_t>();
    // uint32_t seltc4_0 = node["DigitalHalf"]["0"]["SelTC4"].as<uint32_t>();
    // uint32_t seltc4_1 = node["DigitalHalf"]["1"]["SelTC4"].as<uint32_t>();
    uint64_t chmask_0 = 0,  chmask_1 = 0;
    uint32_t multfactor_0 = node["DigitalHalf"]["0"]["MultFactor"].as<uint32_t>();
    uint32_t multfactor_1 = node["DigitalHalf"]["1"]["MultFactor"].as<uint32_t>();
    for(int itot=0;itot<4;itot++) {
      totname = "Tot_P" + std::to_string(itot);
      tot_p_0[itot] = node["DigitalHalf"]["0"][totname].as<uint32_t>();
      tot_p_1[itot] = node["DigitalHalf"]["1"][totname].as<uint32_t>();
      totname = "Tot_TH" + std::to_string(itot);
      tot_th_0[itot] = node["DigitalHalf"]["0"][totname].as<uint32_t>();
      tot_th_1[itot] = node["DigitalHalf"]["1"][totname].as<uint32_t>();
      // tot_p_0[itot] = 0;
      // tot_p_1[itot] = 0;
      // tot_th_0[itot] = 0;
      // tot_th_1[itot] = 0;
    }    
    //std::cout<<"rocname : "<<rocname<<", th_0 : "<<th_0<<", th_1 : "<<th_1<<", chmask_0 : "<<chmask_0<<", chmask_1 : "<<chmask_1<<std::endl;

    // The following are set as input or modulepath
    // rocn = iroc;
    // selTC4 = seltc4_0; half = 0;
    // uint32_t rocid_0 = pck.packRocId(zside, sector, link, det, econt, selTC4, module, rocn, half);
    // selTC4 = seltc4_1; half = 1;
    // uint32_t rocid_1 = pck.packRocId(zside, sector, link, det, econt, selTC4, module, rocn, half);

    
    TPGFEConfiguration::ConfigHfROC hroc_0, hroc_1;      
    hroc_0.setAdcTH(th_0);
    hroc_0.setClrAdcTottrig(chmask_0);
    hroc_0.setMultFactor(multfactor_0);
    for(int itot=0;itot<4;itot++){
      hroc_0.setTotTH(itot, tot_th_0[itot]);
      hroc_0.setTotP(itot, tot_p_0[itot]);
    }

    hroc_1.setAdcTH(th_1);
    hroc_1.setClrAdcTottrig(chmask_1);
    hroc_1.setMultFactor(multfactor_1);
    for(int itot=0;itot<4;itot++){
      hroc_1.setTotTH(itot, tot_th_1[itot]);
      hroc_1.setTotP(itot, tot_p_1[itot]);
    }
    
    hroccfg[rocid_0] = hroc_0;
    hroccfg[rocid_1] = hroc_1;
    
    for(uint32_t ich=0;ich<nchs;ich++){
      uint32_t ihalf = (ich<TPGFEDataformat::HalfHgcrocData::NumberOfChannels)?0:1;
      uint32_t chnl = ich%TPGFEDataformat::HalfHgcrocData::NumberOfChannels;
      uint32_t ped = node["ch"][std::to_string(ich)]["Adc_pedestal"].as<uint32_t>() ;
      if(ihalf==0){
	TPGFEConfiguration::ConfigCh ch_0;
	ch_0.setAdcpedestal(ped);
	hrocchcfg[pck.packChId(rocid_0,chnl)] = ch_0;
      }else{
	TPGFEConfiguration::ConfigCh ch_1;
	ch_1.setAdcpedestal(ped);
	hrocchcfg[pck.packChId(rocid_1,chnl)] = ch_1;
      }
    }//channel loop
  }//end of read ped class
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Configuration::readEconDConfigYaml()
  {
    
    //something like below should loop over ECON-D config(s) for different modules
    uint32_t idx = pck.packModId(zside, sector, link, det, econt, selTC4, module); //we assume same ECONT and ECOND number for a given module
    readEconDConfigYaml(idx);
    
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Configuration::readEconDConfigYaml(uint32_t idx)
  {
    std::cout<<"Configuration::readEconDConfigYaml:: ECOND init config file : "<<EconDfname<<std::endl;
    YAML::Node node(YAML::LoadFile(EconDfname));
    
    //something like below should loop over ECON-D config(s) for different modules
    //uint32_t idx = pck.packModId(zside, sector, link, det, econt, selTC4, module); //we assume same ECONT and ECOND number for a given module
    uint32_t pass_thru_mode = node["RocDaqCtrl"]["Global"]["pass_thru_mode"].as<uint32_t>();
    uint32_t active_erxs = node["RocDaqCtrl"]["Global"]["active_erxs"].as<uint32_t>();
    bool passTM = (pass_thru_mode==1) ? true : false ; 
    TPGFEConfiguration::ConfigEconD econd;
    econd.setPassThrough(passTM);
    econd.setNeRx(active_erxs);
    econDcfg[idx] = econd;
    
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Configuration::readEconTConfigYaml()
  {
    //something like below should loop over ECON-T config(s) for different modules
    uint32_t idx = pck.packModId(zside, sector, link, det, econt, selTC4, module);  
    readEconTConfigYaml(idx);
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Configuration::readEconTConfigYaml(uint32_t idx)
  {
    std::cout<<"Configuration::readEconTConfigYaml:: ECONT init config file : "<<EconTfname<<std::endl;
    YAML::Node node(YAML::LoadFile(EconTfname));
    
    ///////////////////////////
    uint32_t dummyCalib = 0x800; //this corresponds to the default value of unit gain
    uint32_t dummyMux = 0x80;    //default setting not know so set this value to disconnect the corresponding TC
    ///////////////////////////

    bool isFirst = true;
    TPGFEConfiguration::ConfigEconT econt;
    if(auto search_it = econTcfg.find(idx) ;  search_it != econTcfg.end()) isFirst = false;
    
    
    bool hasFound = false;
    uint32_t select = 0;
    if(node["Algorithm"]["Global"]["select"]){
      select = node["Algorithm"]["Global"]["select"].as<uint32_t>();
      hasFound = true;
    }
    if(isFirst) econt.setSelect(select);
    if(hasFound and !isFirst) econTcfg[idx].setSelect(select);
      
    hasFound = false;
    uint32_t density = 0;
    if(node["Algorithm"]["Global"]["density"]){
      density = node["Algorithm"]["Global"]["density"].as<uint32_t>();
      hasFound = true;
    }
    if(isFirst) econt.setDensity(density);
    if(hasFound and !isFirst) econTcfg[idx].setDensity(density);
    
    hasFound = false;
    uint32_t dropLSB = 0;
    if(node["AlgoDroplsb"]["Global"]){
      dropLSB = node["AlgoDroplsb"]["Global"].as<uint32_t>();
      hasFound = true;
    }
    if(isFirst) econt.setDropLSB(dropLSB);
    if(hasFound and !isFirst) econTcfg[idx].setDropLSB(dropLSB);
    
    hasFound = false;
    uint32_t stctype = 0;
    if(node["FmtBuf"]["Global"]["stc_type"]){
      stctype = node["FmtBuf"]["Global"]["stc_type"].as<uint32_t>();
      hasFound = true;
    }
    if(isFirst) econt.setSTCType(stctype);
    if(hasFound and !isFirst) econTcfg[idx].setSTCType(stctype);
				
    hasFound = false;
    uint32_t eporttx_numen = 0;
    if(node["FmtBuf"]["Global"]["eporttx_numen"]){
      eporttx_numen = node["FmtBuf"]["Global"]["eporttx_numen"].as<uint32_t>();
      hasFound = true;
    }
    if(isFirst) econt.setNElinks(eporttx_numen);
    if(hasFound and !isFirst) econTcfg[idx].setNElinks(eporttx_numen);
    
    uint32_t value = 0; 
    if (YAML::Node mux = node["Mux"]) {
      for(uint32_t itc=0;itc<48;itc++){
	std::stringstream ss ;
	//ss << std::setw(2) << std::setfill('0') << itc; //0 padding for TCs between 0-9
	ss << itc;
	if(mux[ss.str().c_str()]){
	  value = mux[ss.str().c_str()].as<uint32_t>();
	}else{
	  value = dummyMux; //set this 
	  std::cerr << "TPGFEconfiguration::Configuration::readEconTConfigYaml(moduleid="<<idx<<") configfile="<<EconTfname<<" : Mux not set for TC " << ss.str().c_str() << std::endl;
	}
	//std::cout << "itc : " << itc <<", ss: " <<ss.str().c_str()<< ", value : " << value << std::endl;
	if(isFirst)
	  econt.setInputMux(itc,value);
	else
	  econTcfg[idx].setInputMux(itc,value);
      }
    }else{
      if(isFirst)
	for(uint32_t itc=0;itc<48;itc++)
	  econt.setInputMux(itc,dummyMux);
    }
    
    if (YAML::Node cal = node["Cal"]) {
      for(uint32_t itc=0;itc<48;itc++){
	std::stringstream ss ;
	ss << std::setw(2) << std::setfill('0') << itc;
	value = (cal[ss.str().c_str()]) ? cal[ss.str().c_str()].as<uint32_t>() : dummyCalib;
	if(isFirst)
	  econt.setCalibration(itc,value);
	else
	  econTcfg[idx].setCalibration(itc,value);
      }
    }else{
      if(isFirst)
	for(uint32_t itc=0;itc<48;itc++)
	  econt.setCalibration(itc,dummyCalib);
    }
    
    if(isFirst) econTcfg[idx] = econt;
    
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Configuration::setPedThZero(){
    uint32_t default_mulfactor = 25;
    
    for(auto const& hrocch : hrocchcfg)
      hrocchcfg[hrocch.first].setAdcpedestal(0);
    
    for(auto const& hroc : hroccfg){
      hroccfg[hroc.first].setAdcTH(0);
      hroccfg[hroc.first].setMultFactor(default_mulfactor);
      for(int itot=0;itot<4;itot++){
	hroccfg[hroc.first].setTotTH(itot, 0);
	hroccfg[hroc.first].setTotP(itot, 0);
      }
    }
    
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Configuration::setPedZero(){
    
    for(auto const& hrocch : hrocchcfg)
      hrocchcfg[hrocch.first].setAdcpedestal(0);
    
    // for(auto const& hroc : hroccfg){
    //   for(int itot=0;itot<4;itot++){
    // 	hroccfg[hroc.first].setTotTH(itot, 0);
    // 	hroccfg[hroc.first].setTotP(itot, 0);
    //   }
    // }    
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Configuration::printCfgPedTh(uint32_t moduleId){
    
    for(auto const& hroc : hroccfg){
      std::cout << "Configuration::printCfgPedTh moduleId: "<< moduleId
		<<", RocId: " <<  hroc.first << std::endl;
      if(moduleId==pck.getModIdFromRocId(uint32_t(hroc.first))){
	hroccfg[hroc.first].print();
	for(auto const& hrocch : hrocchcfg){
	  if(hroc.first==pck.getRocIdFromChId(uint64_t(hrocch.first))){
	    uint32_t ich = pck.getRocChId(hrocch.first) ;
	    hrocchcfg[hrocch.first].print(ich);
	  }
	}
      }
    }
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}


#endif
