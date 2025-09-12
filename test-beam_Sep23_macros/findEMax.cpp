/**********************************************************************
 Created on : 01/07/2024
 Purpose    : Find the maximum energy value for TC,STC,Modsum for a given module type (Silicon or Scintillator)
 Author     : Indranil Das, Visiting Fellow
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

#include "yaml-cpp/yaml.h"
#include <deque>

#include <string>
#include <fstream>
#include <cassert>

#include "TPGFEDataformat.hh"
#include "TPGFEConfiguration.hh"
#include "TPGFEModuleEmulation.hh"

#include "TPGBEDataformat.hh"
#include "Stage1IO.hh"

int main(int argc, char** argv)
{
  
  //===============================================================================================================================
  //input paramaters for testing
  //===============================================================================================================================
  uint32_t zside = 0, sector = 0, link = 0, econt = 0;
  
  //A combinatin of following three will uniquiely identfy the module type in HGCAL
  //(outlook : some modification expected for the ECONTs connected to different types of modules, which are not considered in current tests)
  //det: (0,1) = (Si,Sci)
  //selTC4: (0,1) = (HD,LD);
  //module:
  //Si:==
  //LD: Full (F/0), Top (T/1), Bottom (B/2), Left (L/3), Right (R/4), Five (5/5)
  //HD: Full (F/0), Top (T/1) aka Half Minus, Bottom (B/2) aka Chop Two, Left (L/3), Right(R/4)
  //Sci:==
  //LD: A5A6(0), B11B12(1), C5(2), D8(3), E8(4), G3(5), G5(6), G7(7), G8(8)
  //HD: K6(0), K8(1), K11(2), K12(3), J12(4)
  uint32_t det = 0, selTC4 = 1, module = 0;
  
  uint32_t multfactor = 31;//TOT mult factor other values could be 14 or 8
  uint32_t inputLSB = (selTC4==0) ? 1 : 0; //lsb at the input TC from ROC
  uint32_t dropLSB = 0;  //lsb at the output during the packing [0 to 4]
  uint32_t select = 1 ;  //1 = Super Trigger Cell (STC), 2 = Best Choice (BC)
  //Not in use: 0 = Threshold Sum (TS), 3 = Repeater, 4=Autoencoder (AE).
  uint32_t stc_type = 1; //0 = STC4B(5E+4M), 1 = STC16(5E+4M), 2 = CTC4A(4E+3M), 3 = STC4A(4E+3M), 4 = CTC4B(5E+3M)
  const uint32_t nelinks = 2; //1 = BC1, 2 = BC4, 3 = BC6, 4 = BC9, 5 = BC14..... (see details https://edms.cern.ch/ui/#!master/navigator/document?P:100053490:100430098:subDocs)
  uint32_t calibration = 0x800; //0x400 = 0.5, 0x800 = 1.0, 0xFFF = 1.99951 (max)
  
  // uint32_t maxADC = 0x3FF ; //10 bit input in TPG path
  // uint32_t maxTOT = 0xFFF ; //12 bit input in TPG path
  // uint32_t maxADC = 134 ; 
  // uint32_t maxTOT = 0x0 ; 
  // uint16_t bx = 4;
  // uint32_t ped = 144;
  // uint32_t adcTh = 5;
  
  uint32_t maxADC = 800 ; 
  uint32_t maxTOT = 0x0 ; 
  uint16_t bx = 0;
  uint32_t ped = 0;
  uint32_t adcTh = 0;
  bool hasset_singlech = false;
  
  //===============================================================================================================================
  
  
  //===============================================================================================================================
  //Read channel mapping
  //===============================================================================================================================
  TPGFEConfiguration::Configuration cfgs;
  cfgs.setSiChMapFile("cfgmap/WaferCellMapTraces.txt");
  cfgs.setSciChMapFile("cfgmap/channels_sipmontile_HDtypes.hgcal.txt");
  cfgs.initId();
  cfgs.readSiChMapping();
  cfgs.readSciChMapping();
  cfgs.loadModIdxToNameMapping();
  cfgs.loadMuxMapping();
  
  //June/2024 : Proposal1 : 1+2+11+1+3+1+3+3+1 = 1(+/-z) + 2(120deg sector) + 11(link max ~1848) + 1(Si/Sci) + 3(ECONT number max 7) + 1(LD/HD) + 4[(Si type)/(Sci type)] +  3(ROC number) + 1 (Half) = 27 bits
  //June/2024 : Proposal3 : 1+2+6+6+1+3+1+3+3+1 = 1(+/-z) + 2(120deg sector) + 6(layers) + 6(max 40 per LD/HD/layer for MB) + 1(Si/Sci) + 3(ECONT number max 7) + 1(LD/HD) + 4[(Si type)/(Sci type)] +  3(ROC number) + 1 (Half) = 28 bits  
  //Below we choose for proposal 1
  TPGFEConfiguration::TPGFEIdPacking pck;
  uint32_t moduleId = pck.packModId(zside, sector, link, det, econt, selTC4, module);
  const std::map<std::tuple<uint32_t,uint32_t,uint32_t>,std::string>& modNameMap = cfgs.getModIdxToName();
  const std::string& modName = modNameMap.at(std::make_tuple(pck.getDetType(),pck.getSelTC4(),pck.getModule()));
  const uint32_t nhfrocs = (pck.getDetType()==0)?cfgs.getSiModNhroc(modName):cfgs.getSciModNhroc(modName);
  //===============================================================================================================================
  
  //===============================================================================================================================
  //Set adc pedestal and threshold to zero the way it is read from yaml module file
  //===============================================================================================================================
  std::map<uint32_t,TPGFEConfiguration::ConfigHfROC> hroccfg;
  std::map<uint64_t,TPGFEConfiguration::ConfigCh> hrocchcfg;
  const uint32_t nrocs = ceil(nhfrocs/2.);
  const uint32_t nchs = 2*TPGFEDataformat::HalfHgcrocData::NumberOfChannels;  
  std::cout<<"Check1 : modName: "<<modName<<", nhfrocs : " << nhfrocs << ", nrocs: " << nrocs << std::endl;
  
  for(uint32_t iroc=0;iroc<nrocs;iroc++){
    uint32_t half0 = 0, half1 = 1;
    uint32_t rocid_0 = pck.packRocId(zside, sector, link, det, econt, selTC4, module, iroc, half0);
    uint32_t rocid_1 = pck.packRocId(zside, sector, link, det, econt, selTC4, module, iroc, half1);      
    TPGFEConfiguration::ConfigHfROC hroc_0, hroc_1;
      
    hroc_0.setAdcTH(adcTh);
    hroc_0.setClrAdcTottrig(0);
    hroc_0.setMultFactor(multfactor);
    for(int itot=0;itot<4;itot++){
      hroc_0.setTotTH(itot, 0);
      hroc_0.setTotP(itot, 0);
    }

    hroc_1.setAdcTH(adcTh);
    hroc_1.setClrAdcTottrig(0);
    hroc_1.setMultFactor(multfactor);
    for(int itot=0;itot<4;itot++){
      hroc_1.setTotTH(itot, 0);
      hroc_1.setTotP(itot, 0);
    }
    
    hroccfg[rocid_0] = hroc_0;
    hroccfg[rocid_1] = hroc_1;
    
    for(uint32_t ich=0;ich<nchs;ich++){
      uint32_t ihalf = (ich<TPGFEDataformat::HalfHgcrocData::NumberOfChannels)?0:1;
      uint32_t chnl = ich%TPGFEDataformat::HalfHgcrocData::NumberOfChannels;
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

  cfgs.setRocPara(hroccfg);
  cfgs.setChPara(hrocchcfg);
  //===============================================================================================================================

    
  //===============================================================================================================================
  //Set ECON-D and ECON-T parameters manually for September, 2023 beam-test
  //===============================================================================================================================
  std::map<uint32_t,TPGFEConfiguration::ConfigEconD> econDPar;
  econDPar[moduleId].setPassThrough(true);
  econDPar[moduleId].setNeRx(nhfrocs);
  cfgs.setEconDPara(econDPar);
  econDPar[moduleId].print();
  
  std::map<uint32_t,TPGFEConfiguration::ConfigEconT> econTPar ;
  econTPar[moduleId].setDensity(inputLSB);
  econTPar[moduleId].setDropLSB(dropLSB);
  econTPar[moduleId].setSelect(select);
  econTPar[moduleId].setSTCType(stc_type); //STC type
  econTPar[moduleId].setNElinks(nelinks); //BC type
  for(uint32_t itc=0;itc<48;itc++) econTPar[moduleId].setCalibration(itc,calibration);
  for(uint32_t itc=0;itc<48;itc++) econTPar[moduleId].setInputMux(itc,itc);
  cfgs.setEconTPara(econTPar);
  std::map<uint32_t,TPGFEConfiguration::ConfigEconT>& cfgTPar = cfgs.getEconTPara() ;
  for(const auto& it : cfgTPar) cfgTPar.at(it.first).print();
  //===============================================================================================================================

  
  //===============================================================================================================================
  //Set ADC/TOT for Emulation
  //===============================================================================================================================
  std::map<uint64_t,std::vector<std::pair<uint32_t,TPGFEDataformat::HalfHgcrocData>>> hrocarray; //event,rocId
  uint64_t eventId = 0;
  uint32_t rocn = 0;
  for(uint32_t ihroc = 0 ; ihroc < nhfrocs ; ihroc++){
    uint32_t half = (ihroc%2==0) ? 0 : 1 ;
    TPGFEDataformat::HalfHgcrocChannelData chdata[TPGFEDataformat::HalfHgcrocData::NumberOfChannels];
    for(uint32_t ich = 0 ; ich < TPGFEDataformat::HalfHgcrocData::NumberOfChannels ; ich++){
      //chdata[ich].setTot(maxTOT,0x3);
      //chdata[ich].setAdc(maxADC,0x0);
      if(!hasset_singlech){
	chdata[ich].setAdc(maxADC,0x0);
	hasset_singlech = true;
      }else{
	chdata[ich].setAdc(0,0x0);
      }
    }
    TPGFEDataformat::HalfHgcrocData hrocdata;
    hrocdata.setBx(bx);
    hrocdata.setChannels(chdata);
    pck.setZero();
    hrocarray[eventId].push_back(std::make_pair(pck.packRocId(zside, sector, link, det, econt, selTC4, module, rocn, half),hrocdata));
    if(ihroc%2==1)rocn++;
  }
  //===============================================================================================================================
  
  //===============================================================================================================================
  //Emulation of HGCROC
  //===============================================================================================================================
  //Rearrange the rocdata for emulation input
  std::map<uint32_t,TPGFEDataformat::HalfHgcrocData> rocdata; 
  rocdata.clear();
  for(const auto& data : hrocarray.at(eventId)){
    rocdata[data.first] = data.second ;
  }
  
  //emulation output format
  std::pair<uint32_t,TPGFEDataformat::ModuleTcData> modTcdata;
  
  bool isSim = false; //true for CMSSW simulation and false for beam-test analysis (Note: not yet set for CMSSW)

  ////////////////////// HRCROC emulation //////////////////////////
  TPGFEModuleEmulation::HGCROCTPGEmulation rocTPGEmul(cfgs);
  rocTPGEmul.Emulate(isSim, eventId, moduleId, rocdata, modTcdata);
  /////////////////////////////////////////////////////////////////

 
  //while needed for later inspection
  std::map<uint64_t,std::vector<std::pair<uint32_t,TPGFEDataformat::ModuleTcData>>> modarray; //event,moduleId
  modarray[eventId].push_back(modTcdata);
  //===============================================================================================================================

  //===============================================================================================================================
  //Emulation of ECONT
  //===============================================================================================================================
  //Rearrange the TC output of module for ECONT emulation input
  std::map<uint32_t,TPGFEDataformat::ModuleTcData> moddata;
  for(const auto& data : modarray.at(eventId))
    moddata[data.first] = data.second ;
  
  ////////////////////// ECONT emulation //////////////////////////
  TPGFEModuleEmulation::ECONTEmulation econtEmul(cfgs);
  //econtEmul.setVerbose(true);
  econtEmul.Emulate(isSim, eventId, moduleId, moddata);
  /////////////////////////////////////////////////////////////////
  
  //emulation output format
  TPGFEDataformat::TcModulePacket& TcRawdata = econtEmul.accessTcRawDataPacket();
  //===============================================================================================================================
  
  //===============================================================================================================================
  //Packing for elink
  //===============================================================================================================================
  uint32_t elinkemul[nelinks];
  TPGFEModuleEmulation::ECONTEmulation::convertToElinkData(bx, TcRawdata.second, elinkemul);
  //===============================================================================================================================
  
  //===============================================================================================================================
  // check I/O
  //===============================================================================================================================
  for(const auto& data : hrocarray.at(eventId)){
    const TPGFEDataformat::HalfHgcrocData& hrocdata = data.second ;
    std::cout << "ModuleId: "<<moduleId << ", RocId: "<< data.first << std::endl;
    if(moduleId==pck.getModIdFromRocId(uint32_t(data.first)))
      hrocdata.print();
  }
  
  const TPGFEDataformat::ModuleTcData& modtcdata = moddata.at(moduleId);
  modtcdata.print();

  std::cout <<"\t1: Module " << TcRawdata.first << ", size : " << TcRawdata.second.size() << std::endl;
  if(TcRawdata.second.type()==TPGFEDataformat::BestC) TcRawdata.second.sortCh();
  const TPGFEDataformat::TcRawDataPacket& tcpkt = TcRawdata.second ;
  tcpkt.print();

    //===============================================================================================================================

  
  //===============================================================================================================================
  //Unpacking for elink to TC rawdata
  //===============================================================================================================================
  TPGFEDataformat::TcRawDataPacket vTC1;
  TPGStage1Emulation::Stage1IO::convertElinksToTcRawData(TcRawdata.second.type(), TcRawdata.second.size(), elinkemul, vTC1);
  vTC1.print();
  //===============================================================================================================================


  //===============================================================================================================================
  //TC rawdata to stage1 input
  //===============================================================================================================================
  TPGBEDataformat::UnpackerOutputStreamPair up1;
  TPGStage1Emulation::Stage1IO::convertTcRawDataToUnpackerOutputStreamPair(bx, vTC1, up1);
  up1.print();
  //===============================================================================================================================
  
  hroccfg.clear();
  hrocchcfg.clear();
  econDPar.clear();
  econTPar.clear();
  hrocarray.clear();
  rocdata.clear();
  modarray.clear();
  moddata.clear();
  
  
  return true;
}
