#ifndef Hgcal10gLinkReceiver_EcontEnergies_h
#define Hgcal10gLinkReceiver_EcontEnergies_h

#include <iostream>
#include <iomanip>
#include <cstdint>
#include <vector>

namespace Hgcal10gLinkReceiver {

  void econtEnergies(const uint64_t *p, std::vector<uint16_t> &v) {
    v.resize(24);

    uint64_t e0,e1;
  
    e0=p[0]<<32|(p[1]&0xffffffff);
    e1=p[1]<<48|(p[2]&0xffffffff)<<16|(p[3]&0xffffffff)>>16;
    //std::cout << "e0,e1 = " << std::hex << e0 << "," << e1 << std::endl;
  
    v[ 4]=(e0>> 1)&0x7f;
    v[ 3]=(e0>> 8)&0x7f;
    v[ 2]=(e0>>15)&0x7f;
    v[ 1]=(e0>>22)&0x7f;
    v[ 0]=(e0>>29)&0x7f;
  
    v[11]=(e1    )&0x7f;
    v[10]=(e1>> 7)&0x7f;
    v[ 9]=(e1>>14)&0x7f;
    v[ 8]=(e1>>21)&0x7f;
    v[ 7]=(e1>>28)&0x7f;
    v[ 6]=(e1>>35)&0x7f;
    v[ 5]=(e1>>42)&0x7f;
  
    e0=(p[0]&0xffffffff00000000)    |(p[1]&0xffffffff00000000)>>32;
    e1=(p[1]&0xffffffff00000000)<<16|(p[2]&0xffffffff00000000)>>16|(p[3]&0xffffffff00000000)>>48;
    //std::cout << "e0,e1 = " << std::hex << e0 << "," << e1 << std::endl;
  
    v[16]=(e0>> 1)&0x7f;
    v[15]=(e0>> 8)&0x7f;
    v[14]=(e0>>15)&0x7f;
    v[13]=(e0>>22)&0x7f;
    v[12]=(e0>>29)&0x7f;
  
    v[23]=(e1    )&0x7f;
    v[22]=(e1>> 7)&0x7f;
    v[21]=(e1>>14)&0x7f;
    v[20]=(e1>>21)&0x7f;
    v[19]=(e1>>28)&0x7f;
    v[18]=(e1>>35)&0x7f;
    v[17]=(e1>>42)&0x7f;
  }

  unsigned unpackerEnergies(const uint64_t *p, std::vector<uint16_t> &v) {
    v.resize(24);

    for(unsigned i(0);i<6;i++) {
      v[i   ]=(p[i]    )&0x7f;
      v[i+ 6]=(p[i]>>13)&0x7f;
      v[i+12]=(p[i]>>32)&0x7f;
      v[i+18]=(p[i]>>45)&0x7f;
    }

    return (p[0]>>26)&0xf;                                                                                                                
  }

  class EcontTriggerCellData {
  public:
    EcontTriggerCellData() {
    }
    
    EcontTriggerCellData(uint8_t c, uint8_t e) {
      _channel=c;
      _energy=e;
    }

    unsigned channel() const {
      return _channel;
    }

    unsigned energy() const {
      return _energy;
    }

  private:
    uint8_t _channel;
    uint8_t _energy;
  };
  
  class ModuleEcontData {
  public:
    ModuleEcontData() {
    }

    void setFromUnpacker(uint64_t *p, bool lo=true) {
      _tcd.resize(12);
      
      unsigned offset(lo?0:32);
      for(unsigned i(0);i<6;i++) {
	_tcd[i  ]=EcontTriggerCellData((p[i]>>( 7+offset))&0x3f,(p[i]>>(   offset))&0x7f);
	_tcd[i+6]=EcontTriggerCellData((p[i]>>(20+offset))&0x3f,(p[i]>>(13+offset))&0x7f);
      }
    
      _valid=true;
      _bunchCounter=(p[0]>>26)&0xf;
    }
    
  private:
    bool _valid;
    uint8_t _bunchCounter;
    std::vector<EcontTriggerCellData> _tcd;    
  };
  
}

#endif
