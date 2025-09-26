#ifndef Hgcal10gLinkReceiver_OrbitGenericCheck_h
#define Hgcal10gLinkReceiver_OrbitGenericCheck_h

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

#include "SlinkBoe.h"
#include "SlinkEoe.h"
#include "OrbitReader.h"

class OrbitGenericCheck {
public:
  OrbitGenericCheck() {
  }

  bool runStart(uint32_t nRun, uint32_t sid) {
    return _ct.runStart(nRun,sid);
  }

  bool orbit(const Hgcal10gLinkReceiver::OrbitHeader &oh) {
    return true;
  }

  bool event(const Hgcal10gLinkReceiver::OrbitReaderEvent &ore) {
    return _ct.event(2*ore._ft->fragmentSize(),ore._array);
  }

  bool runStop() {
    return _ct.runStop();
  }

private:
  CheckTypedef _ct;
};

#endif
