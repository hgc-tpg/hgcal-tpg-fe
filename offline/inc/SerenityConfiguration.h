#ifndef Hgcal10gLinkReceiver_SerenityConfiguration_h
#define Hgcal10gLinkReceiver_SerenityConfiguration_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <vector>
#include <string>
#include <cstring>
#include <unordered_map>

#include "RecordConfigured.h"


namespace Hgcal10gLinkReceiver {

  class SerenityConfiguration {

  public:
  
    SerenityConfiguration() { 
    }
    
    virtual ~SerenityConfiguration() {
    }
    
    const std::unordered_map<std::string,uint32_t>& configuration() const {
      return _map;
    }
    
    void setConfiguration(const std::unordered_map<std::string,uint32_t> &m) {
      _map=m;
    }

    void setConfiguration(const RecordConfigured *r) {
      r->configuration(_map);
    }

    /*    
    bool uhalWrite(const std::string &s, uint32_t v) {
      std::cout << "uhalWrite() "
		<< s << ", value = " << v << std::endl;

      auto search=map_.find(s);
      if(search!=map_.end()) {
	search->second.write(_uhalData[search->second._dataWord],v);
	return true;
      }

      return false;
    }
    */
    uint32_t uhalRead(const std::string &s) {
      std::cout << "uhalRead called for "
		<< s << std::endl;

      auto search=accessInfoMap.find(s);
      if(search!=accessInfoMap.end()) {
	return search->second.read(_uhalData[search->second._dataWord]);
      }

      return 0xffffffff;
    }
    
    void print(std::ostream &o=std::cout) {
      o << "SerenityConfiguration::print()" << std::endl;
      o << " Current settings for " << _map.size()
	<< " values:" << std::endl;

      o << std::hex << std::setfill('0');

      for(unsigned i(0);i<_map.size();i++) {
	uint32_t v(uhalRead(_map[i]));
	std::cout << "  " << _map[i] << " = 0x"
		  << std::hex << std::setfill('0')
		  << std::setw(8) << v
		  << std::dec << std::setfill(' ')
		  << " = " << v
		  << std::endl;
      }
      o << std::dec << std::setfill(' ');
    }

  protected:
    std::unordered_map<std::string,uint32_t> map_;
  };

}

#endif
