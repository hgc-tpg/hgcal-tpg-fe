#include <string>
#include <chrono>
#include <thread>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <vector>
#include <unistd.h>

//using namespace Hgcal10gLinkReceiver;

int main(int argc, char *argv[]) {
  if(argc<2) return 1;

  uint64_t word(0);

  std::istringstream sin(argv[1]);
  sin >> std::hex >> word;

  unsigned bx((word>>60)&0xf);
  unsigned ms((word>>52)&0xff);
  
  std::vector<unsigned> cLo;
  for(unsigned i(0);i<8;i++) cLo.push_back((word>>(46-6*i))&0x3f);

  std::vector<unsigned> cHi;
  for(unsigned i(0);i<48;i++) {
    if((word&(uint64_t(1)<<(51-i)))!=0) cHi.push_back(i);
  }

  std::cout << "BX = " << std::setw(2) << bx
	    << ", module sum = " << std::setw(3) << ms << std::endl;
  
  std::cout << "Low  occ format " << std::setw(2) << cLo.size()
	    << " channels =";
  for(unsigned i(0);i<cLo.size();i++)
    std::cout << " " << std::setw(2) << cLo[i];
  std::cout << std::endl;
  
  std::cout << "High occ format " << std::setw(2) << cHi.size()
	    << " channels =";
  for(unsigned i(0);i<cHi.size();i++)
    std::cout << " " << std::setw(2) << cHi[i];
  std::cout << std::endl;  
  
  return 0;
}
