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
  
  std::vector<unsigned> cv;
  for(unsigned i(0);i<3;i++) cv.push_back((word>>(56-4*i))&0xf);

  std::vector<unsigned> ev;
  for(unsigned i(0);i<3;i++) ev.push_back((word>>(39-9*i))&0x1ff);

  unsigned empty(word&0x1fffff);
  
  std::cout << "BX = " << std::setw(2) << bx << std::endl;
  std::cout << "Channels =";
  for(unsigned i(0);i<cv.size();i++)
    std::cout << " " << std::setw(3) << cv[i]+16*i;
  std::cout << std::endl;
  
  std::cout << "Energies =";
  for(unsigned i(0);i<ev.size();i++)
    std::cout << " " << std::setw(3) << ev[i];
  std::cout << std::endl;  
  
  std::cout << "Empty = " << empty << std::endl;
  return 0;
}
