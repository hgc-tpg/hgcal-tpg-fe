////Code from Paul


/*
g++ -std=c++11 \
-I ~/HgcalSoftware/hgcal10glinkreceiver \
-I ~/HgcalSoftware/hgcal10glinkreceiver/common/inc \
-I ~/HgcalSoftware/hgcal10glinkreceiver/offline/inc \
-I ~/HgcalSoftware/hgcal-tpg-fe/inc \
-I ~/HgcalSoftware/hgcal-tpg-fe/TPGStage1Emulation \
-I ~/HgcalSoftware/hgcal-tpg-fe/TPGStage2Emulation \
-I ~/HgcalSoftware/hgcal-tpg-fe/EMPTools \
-I ~/HgcalSoftware/hgcal-tpg-fe/EMPTools/HLS_arbitrary_Precision_Types/include \
src/TpgTimingCheck.cpp -o bin/TpgTimingCheck.exe \
-L /usr/lib64 -l yaml-cpp \
`root-config --libs --cflags`
*/

#ifndef Hgcal10gLinkReceiver_TpgTimingCheck_h
#define Hgcal10gLinkReceiver_TpgTimingCheck_h

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <thread>

#include <yaml-cpp/yaml.h>

#include "TH1D.h"
#include "TH2D.h"
#include "TFileHandler.h"

#include "SlinkBoe.h"
#include "SlinkEoe.h"
#include "TpgSubpacketHeader.h"
#include "OrbitReader.h"

#include "Stage1IO.hh"

void bestChoiceChannels(uint64_t word, std::vector<unsigned> &vLo, std::vector<unsigned> &vHi, bool doPrint=false) {
  vLo.resize(0);
  vHi.resize(0);

  bool done(false);
  for(unsigned i(0);i<8 && !done;i++) {
    unsigned ch((word>>(46-6*i))&0x3f);
    if(i==0) vLo.push_back(ch);
    else if(ch>vLo[i-1] && ch<48) vLo.push_back(ch);
    else done=true;
  }

  for(unsigned i(0);i<48;i++) {
    if((word&(uint64_t(1)<<(51-i)))!=0) vHi.push_back(i);
  }

  if(doPrint) {
    std::cout << "bestChoiceChannels 0x"
	      << std::hex << std::setfill('0')
	      << std::setw(16) << word
	      << std::dec << std::setfill(' ')
	      << std::endl
	      << " Low  occupancy: up to " << vLo.size() << " channels found =";
    
    for(unsigned i(0);i<vLo.size();i++)
    std::cout << " " << std::setw(2) << vLo[i];
    std::cout << std::endl;
  
    std::cout << " High occupancy: " << std::setw(2) << vHi.size() << " channels found =";

    for(unsigned i(0);i<vHi.size();i++)
      std::cout << " " << std::setw(2) << vHi[i];
    std::cout << std::endl;  
    std::cout << std::endl;  
  }
}

void bestChoiceChannels(uint32_t w0, uint32_t w1, std::vector<unsigned> &vLo, std::vector<unsigned> &vHi, bool doPrint=false) {
  uint64_t word(w0);
  word=(word<<32)|w1;
  bestChoiceChannels(word,vLo,vHi,doPrint);
}

void bestChoiceChannels(const uint32_t *w, std::vector<unsigned> &vLo, std::vector<unsigned> &vHi, bool doPrint=false) {
  uint64_t word(w[0]);
  word=(word<<32)|w[1];
  bestChoiceChannels(word,vLo,vHi,doPrint);
}

class TpgTimingCheck {
public:
  TpgTimingCheck() {
  }

  bool runStart(uint32_t nRun, uint32_t sid) {
    std::cout << "TpgTimingCheck: called runStart for run " << nRun
	      <<", source id " << sid << std::endl;
    
    _sourceId=sid;
    _nEvents=0;
    _initialOrbit=0;

    for(unsigned i(0);i<128;i++) {
      for(unsigned j(0);j<16;j++) {
	_hPlot[i][j]=nullptr;
      }
    }
    
    _tfh.initialise(std::string("TpgTimingCheckRun")+std::to_string(nRun)+"Sid"+std::to_string(sid));

    lNode=YAML::LoadFile("TpgTimingModules.yaml");
    std::cout << lNode << std::endl;
    
     return true;
  }

  bool event(unsigned n64, const uint64_t *p) {
    //std::cout << "TpgTimingCheck: called event" << std::endl;
    
    _nEvents++;
    
    Hgcal10gLinkReceiver::SlinkBoe *b(nullptr);
    if(n64>=2) {
      b=(Hgcal10gLinkReceiver::SlinkBoe*)p;
      //b->print();
      //std::cout << std::endl;
    }
    
    Hgcal10gLinkReceiver::SlinkEoe *e(nullptr);
    if(n64>=4) {
      e=(Hgcal10gLinkReceiver::SlinkEoe*)(p+n64-2);
      //e->print();
      //std::cout << std::endl;
    }

    uint32_t vElinks[14];
    
    const Hgcal10gLinkReceiver::TpgSubpacketHeader *tsh(reinterpret_cast<const Hgcal10gLinkReceiver::TpgSubpacketHeader*>(b+1));
    const Hgcal10gLinkReceiver::TpgSubpacketHeader *tshEnd(reinterpret_cast<const Hgcal10gLinkReceiver::TpgSubpacketHeader*>(e-1));

    if(!tsh->validPattern()) {
      tsh->print();
      return false;
    }
    tsh=tsh->nextSubpacketHeader();
    
    bool done(false);
    while(tsh<=tshEnd && !done) {
      if(!tsh->validPattern()) {
	done=true;
	
      } else {
	if((tsh->channelId()%2)==0) { // RX only
	  //tsh->print();

	  unsigned chan(tsh->channelId()/2);
	  if(_hPlot[chan][0]==nullptr) {
	    tsh->print();
	    for(unsigned m(0);m<lNode["Links"][chan]["Modules"].size();m++) {
	      _hChan[chan][m]=new TH2D((std::string("Chan")+std::to_string(tsh->channelId()/2)+"Module"+std::to_string(m)).c_str(),
					";BX;Channel",7,-3.5,3.5,48,0,48);
	      _hChanE[chan][m]=new TH2D((std::string("ChanE")+std::to_string(tsh->channelId()/2)+"Module"+std::to_string(m)).c_str(),
					";BX;Channel",7,-3.5,3.5,48,0,48);
	      for(unsigned j(0);j<48;j++) {
		if(j<10) _hEnergy[chan][m][j]=new TH2D((std::string("Energy")+std::to_string(tsh->channelId()/2)+"Module"+std::to_string(m)+"Channel0"+std::to_string(j)).c_str(),
						       ";BX;Energy",7,-3.5,3.5,128,0,128);
		else     _hEnergy[chan][m][j]=new TH2D((std::string("Energy")+std::to_string(tsh->channelId()/2)+"Module"+std::to_string(m)+"Channel" +std::to_string(j)).c_str(),
						       ";BX;Energy",7,-3.5,3.5,128,0,128);
	      }
	    }
	    
	    for(unsigned j(0);j<8;j++) {
	      _hBx[chan][j]=new TH2D((std::string("BX")+std::to_string(tsh->channelId()/2)+"Word"+std::to_string(j)).c_str(),
				     ";BX;TPG BX",3600,0,3600,16,0,16);
	    }
	    
	    for(unsigned j(0);j<16;j++) {
	      _hPlot[chan][j]=new TH2D((std::string("Plot")+std::to_string(tsh->channelId()/2)+"Word"+std::to_string(j)).c_str(),
				       ";BX;Value",7,-3.5,3.5,1024,0,65536);
	      _hDiff[chan][j]=new TH2D((std::string("Plot")+std::to_string(tsh->channelId()/2)+"Word"+std::to_string(j)+"Diff").c_str(),
				       ";BX;Value",7,-3.5,3.5,1024,-65536,65536);
	    }
	  }
	  
	  for(unsigned bx(0);bx<tsh->numberOfBxs();bx++) {
	    int diffBx(int(bx)-int(tsh->numberOfBxs()/2));
	    const uint64_t *d((const uint64_t*)(tsh+1+                    bx*tsh->numberOfWordsPerBx()));
	    const uint64_t *c((const uint64_t*)(tsh+1+(tsh->numberOfBxs()/2)*tsh->numberOfWordsPerBx()));
	    
	    const uint32_t *s((const uint32_t*)d);

	    //std::cout << "d,s = " << std::hex << std::setw(16) << d[0]
	    //	      << " " << std::setw(8) << s[0] << " " << std::setw(8) << s[1] << std::endl;

	    int realBx(int(e->bxId())+diffBx);
	    if(realBx<   1) realBx+=3564;
	    if(realBx>3564) realBx-=3564;
	    
	    for(unsigned j(0);j<tsh->numberOfWordsPerBx();j++) {
	      _hBx[chan][2*j  ]->Fill(realBx,(d[j]>>28)&0xf);
	      _hBx[chan][2*j+1]->Fill(realBx,(d[j]>>60)&0xf);
	      
	      _hPlot[chan][4*j+3]->Fill(diffBx,(d[j]    )&0xffff);
	      _hPlot[chan][4*j+2]->Fill(diffBx,(d[j]>>16)&0xffff);
	      _hPlot[chan][4*j+1]->Fill(diffBx,(d[j]>>32)&0xffff);
	      _hPlot[chan][4*j+0]->Fill(diffBx,(d[j]>>48)&0xffff);
	      
	      if(false && diffBx==0) {
		std::cout << "bx " << bx << ", n " << unsigned(tsh->numberOfBxs()) << " diffBx " << diffBx << std::endl;
		tsh->print();
		if((int((d[j]    )&0xffff)-int((c[j]    )&0xffff))!=0) std::cout << "DIFF!!! " << d[j] << ", " << c[j] << std::endl;
		if((int((d[j]>>16)&0xffff)-int((c[j]>>16)&0xffff))!=0) std::cout << "DIFF!!! " << d[j] << ", " << c[j] << std::endl;
		if((int((d[j]>>32)&0xffff)-int((c[j]>>32)&0xffff))!=0) std::cout << "DIFF!!! " << d[j] << ", " << c[j] << std::endl;
		if((int((d[j]>>48)&0xffff)-int((c[j]>>48)&0xffff))!=0) std::cout << "DIFF!!! " << d[j] << ", " << c[j] << std::endl;
	      }

	      _hDiff[chan][4*j+3]->Fill(diffBx,int((d[j]    )&0xffff)-int((c[j]    )&0xffff));
	      _hDiff[chan][4*j+2]->Fill(diffBx,int((d[j]>>16)&0xffff)-int((c[j]>>16)&0xffff));
	      _hDiff[chan][4*j+1]->Fill(diffBx,int((d[j]>>32)&0xffff)-int((c[j]>>32)&0xffff));
	      _hDiff[chan][4*j+0]->Fill(diffBx,int((d[j]>>48)&0xffff)-int((c[j]>>48)&0xffff));
	    }
	    
	    if(lNode["Links"][chan]["Modules"]) {
	      //std::cout << "LINK " << chan << std::endl << lNode["Links"][chan]["Modules"] << std::endl;
	      for(unsigned m(0);m<lNode["Links"][chan]["Modules"].size();m++) {
		if(lNode["Links"][chan]["Modules"][m]["Format"].as<std::string>()=="STC4A") {
		  for(unsigned j(0);j<lNode["Links"][chan]["Modules"][m]["Elinks"].size();j++) {
		    unsigned n(lNode["Links"][chan]["Modules"][m]["Elinks"][j].as<unsigned>());
		    //vElinks[j]=s[(n%2)==0?n+1:n-1];
		    vElinks[j]=s[n];
		  }
		  TPGFEDataformat::TcRawDataPacket rdp;
		  TPGStage1Emulation::Stage1IO::convertElinksToTcRawData(TPGFEDataformat::STC4A,lNode["Links"][chan]["Modules"][m]["NumberOfTcs"].as<unsigned>(),vElinks,rdp);
		  rdp.print();
		  const std::vector<TPGFEDataformat::TcRawData> &vTc(rdp.getTcData());
		  for(unsigned i(0);i<lNode["Links"][chan]["Modules"][m]["NumberOfTcs"].as<unsigned>();i++) {
		    _hChan[chan][m]->Fill(diffBx,vTc[i].address()+4*i);
		    _hChanE[chan][m]->Fill(diffBx,vTc[i].address()+4*i,vTc[i].energy());
		  }
		}
		
		else if(lNode["Links"][chan]["Modules"][m]["Format"].as<std::string>()=="STC16") {
		  for(unsigned j(0);j<lNode["Links"][chan]["Modules"][m]["Elinks"].size();j++) {
		    unsigned n(lNode["Links"][chan]["Modules"][m]["Elinks"][j].as<unsigned>());
		    //vElinks[j]=s[(n%2)==0?n+1:n-1];
		    vElinks[j]=s[n];
		  }
		  TPGFEDataformat::TcRawDataPacket rdp;
		  TPGStage1Emulation::Stage1IO::convertElinksToTcRawData(TPGFEDataformat::STC16,lNode["Links"][chan]["Modules"][m]["NumberOfTcs"].as<unsigned>(),vElinks,rdp);
		  rdp.print();
		}
		
		else if(lNode["Links"][chan]["Modules"][m]["Format"].as<std::string>()=="BestChoice") {
		  for(unsigned j(0);j<lNode["Links"][chan]["Modules"][m]["Elinks"].size();j++) {
		    unsigned n(lNode["Links"][chan]["Modules"][m]["Elinks"][j].as<unsigned>());
		    //vElinks[j]=s[(n%2)==0?n+1:n-1];
		    vElinks[j]=s[n];
		  }
		  
		  std::vector<unsigned> vLo,vHi;
		  bestChoiceChannels(vElinks[0],vElinks[1],vLo,vHi,_nEvents==1);

		  bool validCh(true);
		  //if(false && (chan==100 || chan==102 || chan==106 || chan==112 || chan==118)) {
		  if(lNode["Links"][chan]["Modules"][m]["NumberOfTcs"].as<unsigned>()>=9) {
		    //std::cout << chan << " Expected hi occ " << lNode["Links"][chan]["Modules"][m]["NumberOfTcs"].as<unsigned>() << std::endl;
		    if(vHi.size()!=lNode["Links"][chan]["Modules"][m]["NumberOfTcs"].as<unsigned>()) {
		      std::cout << std::endl << "===>> Event " << _nEvents << std::endl;
		      std::cout << std::endl << lNode["Links"][chan]["Modules"][m] << std::endl << std::endl;
		      bestChoiceChannels(vElinks[0],vElinks[1],vLo,vHi,true);
		      validCh=false;
		      //assert(false);
		    }
		  } else {
		    //std::cout << chan << " Expected lo occ " << lNode["Links"][chan]["Modules"][m]["NumberOfTcs"].as<unsigned>() << std::endl;
		    //bestChoiceChannels(vElinks[0],vElinks[1],vLo,vHi,true);
		    if(vLo.size()<lNode["Links"][chan]["Modules"][m]["NumberOfTcs"].as<unsigned>()) {
		      std::cout << std::endl << "===>> Event " << _nEvents << std::endl;
		      std::cout << std::endl << lNode["Links"][chan]["Modules"][m] << std::endl << std::endl;
		      bestChoiceChannels(vElinks[0],vElinks[1],vLo,vHi,true);
		      validCh=false;
		      //assert(false);
		    }
		  }
		  //}
		  /*		  
		  if(lNode["Links"][chan]["Modules"][m]["NumberOfTcs"].as<unsigned>()>=9) {
		    for(unsigned i(0);i<vHi.size();i++) {
		      _hChan[chan][m]->Fill(diffBx,vHi[i]);		      
		    }
		  } else {
		    for(unsigned i(0);i<vLo.size();i++) {
		      _hChan[chan][m]->Fill(diffBx,vLo[i]);		      
		    }
		  }
		  */
		  if(validCh) {
				      
		    if(lNode["Links"][chan]["Modules"][m]["NumberOfTcs"].as<unsigned>()>=9) {
		      TPGFEDataformat::TcRawDataPacket rdp;
		      TPGStage1Emulation::Stage1IO::convertElinksToTcRawData(TPGFEDataformat::BestC,lNode["Links"][chan]["Modules"][m]["NumberOfTcs"].as<unsigned>(),vElinks,rdp);
		      if(_nEvents==1) rdp.print();
		      const std::vector<TPGFEDataformat::TcRawData> &vTc(rdp.getTcData());
		      for(unsigned i(0);i<vHi.size();i++) {
			_hChan[chan][m]->Fill(diffBx,vTc[i].address());
			_hChanE[chan][m]->Fill(diffBx,vTc[i].address(),vTc[i].energy());
			_hEnergy[chan][m][vTc[i].address()]->Fill(diffBx,vTc[i].energy());
		      }
		    } else {
		      TPGFEDataformat::TcRawDataPacket rdp;
		      TPGStage1Emulation::Stage1IO::convertElinksToTcRawData(TPGFEDataformat::BestC,lNode["Links"][chan]["Modules"][m]["NumberOfTcs"].as<unsigned>(),vElinks,rdp);
		      if(_nEvents==1) rdp.print();
		      const std::vector<TPGFEDataformat::TcRawData> &vTc(rdp.getTcData());
		      for(unsigned i(0);i<vLo.size();i++) {
			_hChan[chan][m]->Fill(diffBx,vTc[i].address());
			_hChanE[chan][m]->Fill(diffBx,vTc[i].address(),vTc[i].energy());
			_hEnergy[chan][m][vTc[i].address()]->Fill(diffBx,vTc[i].energy());
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	
	tsh=tsh->nextSubpacketHeader();
      }
    }
    //exit(0);
    return true;
  }

  bool runStop() {
    std::cout << "TpgTimingCheck: called runStop" << std::endl;
    return true;
  }
  
private:
  TFileHandler _tfh;
  YAML::Node lNode;
  
  TH2D *_hBx[128][8];
  
  TH2D *_hChan[128][8];
  TH2D *_hChanE[128][8];
  TH2D *_hEnergy[128][8][48];

  TH2D *_hPlot[128][16];
  TH2D *_hDiff[128][16];
  
  uint32_t _sourceId;
  uint32_t _nEvents;
  uint32_t _initialOrbit;
};

#endif

//#include "TpgTimingCheck.h"

typedef TpgTimingCheck CheckTypedef;

#include "OrbitGenericCheck.h"

typedef OrbitGenericCheck OrbitCheckTypedef;

#include "OrbitCheck.hxx"
