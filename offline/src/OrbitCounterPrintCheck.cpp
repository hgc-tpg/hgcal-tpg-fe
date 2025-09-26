/*
  g++ -std=c++11 -I ../hgcal10glinkreceiver -I ../hgcal10glinkreceiver/common/inc -I ../hgcal10glinkreceiver/offline/inc ../hgcal10glinkreceiver/offline/src/OrbitCounterCheck.cpp -o OrbitCounterCheck.exe
*/

#include "CounterPrintCheck.h"

typedef CounterPrintCheck CheckTypedef;

#include "OrbitGenericCheck.h"

typedef OrbitGenericCheck OrbitCheckTypedef;

#include "OrbitCheck.hxx"
