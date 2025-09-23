#all: emul_test-beam_Sep23 findEMax read_econt_Jul24 GenerateEmpRxFile
HOME=/home/hep/idas
CPPFLAGS=-I$(HOME)/Software/yaml-cpp/include -I common/inc -I inc -I TPGStage1Emulation/ -I TPGFEEmulation/ -I`root-config --incdir` 
BOOST=/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/boost/1.80.0-87b5de10acd2f2c8a325345ad058b814
LDFLAGS=-L$(HOME)/Software/yaml-cpp/lib64 -L$(BOOST)/lib
CPPFLAGSSTAGE2=-I inc/ -I TPGStage2Emulation/ -I EMPTools/ -I EMPTools/CMSSWCode -I EMPTools/HLS_arbitrary_Precision_Types/include/ -I`root-config --incdir` 

all: stage2SemiEmulator.exe stage2BkgRate.exe CalcResolution.exe #stage2CPCalc.exe #stage2BkgRate.exe #CalcResolution.exe stage2SemiEmulator.exe stage2CPCalc.exe #writeStage2MultEMPRx_TWCP.exe  compareStage2TWCP_FWvsEmul_allBxs.exe stage2TWCPEmulEMP.exe stage2SemiEmulationTCBitsAccumulatorFW.exe stage2SemiEmulationTCBits.exe stage2HtoTauTauEnergyCorrelation.exe testfp.exe  #compareCPPyCpp.exe stage2CPCalc.exe readCPPyEmulOut.exe writeStage2EMPRx_TWCP.exe writeStage2MultEMPRx_TWCP.exe compareStage2TWCP_FWvsEmul.exe compareStage2TWCP_FWvsEmul_allBxs.exe stage2TWCPEmulEMP.exe  #compareStage2TowerFWvsEmul.exe stage2TowerEmulEMP.exe readlpGBTpairEvents.exe  #checkStage2Config.exe testStage1Stage2Emulation.exe  writeStage2EMPRx.exe EmulTowerPRRTest.exe #readNTuple.exe  fillInputData.exe testStage2SemiClustering.exe stage2HtoTauTauEnergyCorrelation.exe ntupleMCInfo.exe vbfjet.exe  TowerPreEmulTest.exe  TestUnpackerTCProcInterface.exe #findFixedpattern.exe findEMax.exe tpgdata_3T_tcproc.exe  tpgdata_3T_fe.exe  scanadc_Sep24.exe dump_event.exe tpgdata_2T_fe.exe emul_Sep24.exe emul_3T_Sep24.exe  #loop_emul_Sep24.exe scanconfigval_Sep24.exe emul_Sep24.exe validateFixedADC.exe # findEMax.exe GenerateEmpRxFile.exe dump_event.exe emul_test-beam_Sep23.exe 

emul_test-beam_Sep23.exe:  test-beam_Sep23_macros/emul_test-beam_Sep23.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep23_macros/emul_test-beam_Sep23.cpp  -l yaml-cpp `root-config --libs --cflags` -o emul_test-beam_Sep23.exe

findEMax.exe: test-beam_Sep23_macros/findEMax.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh TPGStage2Emulation/*.hh common/inc/*.h
	g++  $(LDFLAGS) $(CPPFLAGS) $(CPPFLAGSSTAGE2) test-beam_Sep23_macros/findEMax.cpp  -l yaml-cpp `root-config --libs --cflags` -o findEMax.exe

read_econt_Jul24.exe: test-beam_Aug24_macros/read_econt_Jul24.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ -I common/inc -I inc -I . test-beam_Aug24_macros/read_econt_Jul24.cpp  -l yaml-cpp `root-config --libs --cflags` -o read_econt_Jul24.exe

validation_lpGBTs.exe: test-beam_Aug24_macros/validation_lpGBTs.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ -I common/inc -I inc test-beam_Aug24_macros/validation_lpGBTs.cpp  -l yaml-cpp `root-config --libs --cflags` -o validation_lpGBTs.exe

emul_Jul24.exe: test-beam_Aug24_macros/emul_Jul24.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) -I. TPGStage1Emulation/HGCalLayer1PhiOrderFwImpl.cc $(CPPFLAGS) test-beam_Aug24_macros/emul_Jul24.cpp  -l yaml-cpp `root-config --libs --cflags` -o emul_Jul24.exe -lm

tpgdata_3T_tcproc.exe: test-beam_Sep24_macros/tpgdata_3T_tcproc.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) -I. TPGStage1Emulation/HGCalLayer1PhiOrderFwImpl.cc $(CPPFLAGS) test-beam_Sep24_macros/tpgdata_3T_tcproc.cpp  -l yaml-cpp `root-config --libs --cflags` -o tpgdata_3T_tcproc.exe -lm

tpgdata_3T_fe.exe: test-beam_Sep24_macros/tpgdata_3T_fe.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep24_macros/tpgdata_3T_fe.cpp  -l yaml-cpp `root-config --libs --cflags` -o tpgdata_3T_fe.exe -lm

tpgdata_2T_fe.exe: test-beam_Sep24_macros/tpgdata_2T_fe.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep24_macros/tpgdata_2T_fe.cpp  -l yaml-cpp `root-config --libs --cflags` -o tpgdata_2T_fe.exe -lm

findFixedpattern.exe: test-beam_Sep24_macros/findFixedpattern.C inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep24_macros/findFixedpattern.C  -l yaml-cpp `root-config --libs --cflags` -o findFixedpattern.exe -lm

scanconfigval_Sep24.exe: test-beam_Sep24_macros/scanconfigval_Sep24.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep24_macros/scanconfigval_Sep24.cpp  -l yaml-cpp `root-config --libs --cflags` -o scanconfigval_Sep24.exe -lm

scanadc_Sep24.exe: test-beam_Sep24_macros/scanadc_Sep24.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep24_macros/scanadc_Sep24.cpp  -l yaml-cpp `root-config --libs --cflags` -o scanadc_Sep24.exe -lm

emul_Sep24.exe: test-beam_Sep24_macros/emul_Sep24.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep24_macros/emul_Sep24.cpp  -l yaml-cpp `root-config --libs --cflags` -o emul_Sep24.exe -lm

emul_3T_Sep24.exe: test-beam_Sep24_macros/emul_3T_Sep24.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep24_macros/emul_3T_Sep24.cpp  -l yaml-cpp `root-config --libs --cflags` -o emul_3T_Sep24.exe -lm

loop_emul_Sep24.exe: test-beam_Sep24_macros/loop_emul_Sep24.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep24_macros/loop_emul_Sep24.cpp  -l yaml-cpp `root-config --libs --cflags` -o loop_emul_Sep24.exe -lm

validateFixedADC.exe: test-beam_Sep24_macros/validateFixedADC.C inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep24_macros/validateFixedADC.C  -l yaml-cpp `root-config --libs --cflags` -o validateFixedADC.exe -lm

dump_event.exe: test-beam_Sep24_macros/dump_event.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep24_macros/dump_event.cpp  -l yaml-cpp `root-config --libs --cflags` -o dump_event.exe

GenerateEmpRxFile.exe: TPGStage1Emulation/GenerateEmpRxFile.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS)  TPGStage1Emulation/GenerateEmpRxFile.cpp  -l yaml-cpp `root-config --libs --cflags` -o GenerateEmpRxFile.exe

TowerPreEmulTest.exe: TPGStage1Emulation/TowerSums.h TPGStage1Emulation/Utilities.h
	g++ $(CPPFLAGS) -c TPGStage1Emulation/Utilities.cpp -o Utilities.o
	g++ $(CPPFLAGS) -c TPGStage1Emulation/TowerSums.cpp -o TowerSums.o	
	g++ $(CPPFLAGS) stage1-PRR/TowerPreEmulTest.cpp  Utilities.o TowerSums.o -o TowerPreEmulTest.exe
	rm *.o

TestUnpackerTCProcInterface.exe: stage1-PRR/TestUnpackerTCProcInterface.cpp TPGFEEmulation/*.hh  TPGStage1Emulation/*.hh common/inc/*.h inc/*.*
	g++ -I TPGStage1Emulation/ -I. TPGStage1Emulation/HGCalLayer1PhiOrderFwImpl.cc -I. -I common/inc -I inc -I TPGFEEmulation/ -I TPGStage2Emulation/ -I EMPTools/ stage1-PRR/TestUnpackerTCProcInterface.cpp -L /opt/local/lib  -l yaml-cpp `root-config --libs --cflags` -I`root-config --incdir` EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* -IEMPTools/CMSSWCode/ -IEMPTools/HLS_arbitrary_Precision_Types/include/ -I$(BOOST)/include -lboost_iostreams -lz -llzma $(LDFLAGS) $(CPPFLAGS) `root-config --libs --cflags` -o TestUnpackerTCProcInterface.exe

UnpackerTCProcInterfaceEMPAndDirect.exe: stage1-PRR/UnpackerTCProcInterfaceEMPAndDirect.cpp TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h inc/*.* TPGStage2Emulation/*.hh
	g++ -I TPGStage1Emulation/ -I. TPGStage1Emulation/HGCalLayer1PhiOrderFwImpl.cc -I. -I common/inc -I inc -I TPGFEEmulation/ -I TPGStage2Emulation/ -I EMPTools/ stage1-PRR/UnpackerTCProcInterfaceEMPAndDirect.cpp -L /opt/local/lib  -l yaml-cpp `root-config --libs --cflags` -I`root-config --incdir` EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* -IEMPTools/CMSSWCode/ -IEMPTools/HLS_arbitrary_Precision_Types/include/ -I$(BOOST)/include -lboost_iostreams -lz -llzma $(LDFLAGS) $(CPPFLAGS) `root-config --libs --cflags` -o UnpackerTCProcInterfaceEMPAndDirect.exe

# TestEMP: EMPTools/test.cpp EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* $(BOOST)/include
# 	g++  EMPTools/test.cpp EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* -IEMPTools/CMSSWCode/ -IEMPTools/HLS_arbitrary_Precision_Types/include/ -I$(BOOST)/include -lboost_iostreams -lz -llzma -o TestEMP.exe

TestEMP: EMPTools/test.cpp EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/*
	g++  EMPTools/test.cpp EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* -IEMPTools/CMSSWCode/ -IEMPTools/HLS_arbitrary_Precision_Types/include/  -lboost_iostreams -lz -llzma -o TestEMP.exe

EmulTowerPRRTest.exe: stage1-PRR/EmulTowerPRRTest.cpp EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/*
	g++ $(CPPFLAGS) -c TPGStage1Emulation/Utilities.cpp -o Utilities.o
	g++ $(CPPFLAGS) -c TPGStage1Emulation/TowerSums.cpp -o TowerSums.o	
	g++ $(LDFLAGS) $(CPPFLAGS) stage1-PRR/EmulTowerPRRTest.cpp EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* -I$(BOOST)/include -IEMPTools/CMSSWCode/ -IEMPTools/HLS_arbitrary_Precision_Types/include/ `root-config --libs --cflags` -lboost_iostreams -lz -llzma -l yaml-cpp -lm Utilities.o TowerSums.o -o EmulTowerPRRTest.exe
	rm *.o	

readNTuple.exe: test-stage2_Nov24/readNTuple.cpp inc/*.*  TPGStage2Emulation/*.hh
	g++ $(CPPFLAGSSTAGE2) test-stage2_Nov24/readNTuple.cpp `root-config --libs --cflags` -o readNTuple.exe

fillInputData.exe: test-stage2_Nov24/fillInputData.cpp inc/*.*  TPGStage2Emulation/*.hh
	g++ $(CPPFLAGSSTAGE2) test-stage2_Nov24/fillInputData.cpp `root-config --libs --cflags` -o fillInputData.exe

testStage2SemiClustering.exe: test-stage2_Nov24/testStage2SemiClustering.cpp inc/*.*  TPGStage2Emulation/*.hh
	g++ $(LDFLAGS) $(CPPFLAGSSTAGE2) test-stage2_Nov24/testStage2SemiClustering.cpp EMPTools/CMSSWCode/L1Trigger/L1THGCal/src/backend_emulator/*cc `root-config --libs --cflags` -l yaml-cpp -o testStage2SemiClustering.exe

ntupleMCInfo.exe: test-stage2_Nov24/ntupleMCInfo.cpp inc/*.*  TPGStage2Emulation/*.hh
	g++ $(CPPFLAGSSTAGE2) test-stage2_Nov24/ntupleMCInfo.cpp `root-config --libs --cflags` -o ntupleMCInfo.exe -lEG

stage2HtoTauTauEnergyCorrelation.exe: test-stage2_Nov24/stage2HtoTauTauEnergyCorrelation.cpp inc/*.*  TPGStage2Emulation/*.hh
	g++ $(CPPFLAGSSTAGE2) test-stage2_Nov24/stage2HtoTauTauEnergyCorrelation.cpp `root-config --libs --cflags` -o stage2HtoTauTauEnergyCorrelation.exe -lEG $(LDFLAGS) $(CPPFLAGS) -l yaml-cpp

vbfjet.exe: test-stage2_Nov24/vbfjet.cpp inc/*.*  TPGStage2Emulation/*.hh
	g++ $(CPPFLAGSSTAGE2) test-stage2_Nov24/vbfjet.cpp `root-config --libs --cflags` -o vbfjet.exe -lEG

writeStage2EMPRx.exe: test-stage2_Nov24/writeStage2EMPRx.cpp inc/*.*  TPGStage2Emulation/*.hh
	g++ $(CPPFLAGSSTAGE2) EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* test-stage2_Nov24/writeStage2EMPRx.cpp `root-config --libs --cflags` -L$(BOOST)/lib  -lboost_iostreams -lz -llzma -o writeStage2EMPRx.exe

stage2TowerEmulEMP.exe: test-stage2_Nov24/stage2TowerEmulEMP.cpp inc/*.*  TPGStage2Emulation/*.hh
	g++ $(CPPFLAGSSTAGE2) $(LDFLAGS) $(CPPFLAGS) EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* test-stage2_Nov24/stage2TowerEmulEMP.cpp `root-config --libs --cflags` -L$(BOOST)/lib  -lboost_iostreams -lz -llzma -l yaml-cpp -o stage2TowerEmulEMP.exe


EmulEMPCompTC.exe: stage1-PRR/EmulEMPCompTC.cpp inc/*.*  
	g++ $(CPPFLAGSSTAGE2) $(LDFLAGS) $(CPPFLAGS) EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* stage1-PRR/EmulEMPCompTC.cpp `root-config --libs --cflags` -L$(BOOST)/lib  -lboost_iostreams -lboost_program_options -lz -llzma -l yaml-cpp -o EmulEMPCompTC.exe

testStage1Stage2Emulation.exe: combined-emulation/testStage1Stage2Emulation.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh TPGStage2Emulation/*.hh TPGEmulation/*.hh
	g++ $(CPPFLAGSSTAGE2) $(CPPFLAGS) $(LDFLAGS) -I TPGEmulation/ -I TPGFEEmulation/  EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* combined-emulation/testStage1Stage2Emulation.cpp `root-config --libs --cflags` -L$(BOOST)/lib  -lboost_iostreams -lz -llzma -l yaml-cpp -o testStage1Stage2Emulation.exe

checkStage2Config.exe: test-stage2_Nov24/checkStage2Config.cpp inc/*.* TPGFEEmulation/*.hh TPGStage1Emulation/*.hh TPGStage2Emulation/*.hh TPGEmulation/*.hh
	g++ $(CPPFLAGSSTAGE2) $(CPPFLAGS) $(LDFLAGS) -I TPGEmulation/ -I TPGFEEmulation/  EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* test-stage2_Nov24/checkStage2Config.cpp `root-config --libs --cflags` -L$(BOOST)/lib  -lboost_iostreams -lz -llzma -l yaml-cpp -o checkStage2Config.exe

readlpGBTpairEvents.exe: stage1-PRR/readlpGBTpairEvents.cpp EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/*
	g++ $(LDFLAGS) $(CPPFLAGS) stage1-PRR/readlpGBTpairEvents.cpp EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* -I$(BOOST)/include -IEMPTools/CMSSWCode/ -IEMPTools/HLS_arbitrary_Precision_Types/include/ `root-config --libs --cflags` -lboost_iostreams -lz -llzma -l yaml-cpp -lm -o readlpGBTpairEvents.exe

compareStage2TowerFWvsEmul.exe: test-stage2_Nov24/compareStage2TowerFWvsEmul.cpp EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/
	g++  $(CPPFLAGSSTAGE2) EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* test-stage2_Nov24/compareStage2TowerFWvsEmul.cpp `root-config --libs --cflags` -L$(BOOST)/lib  -lboost_iostreams -lz -llzma -o compareStage2TowerFWvsEmul.exe

stage2CPCalc.exe: test-stage2_Nov24/stage2CPCalc.cpp inc/*.* EMPTools/CMSSWCode/DataFormats/L1THGCal/interface/HGCalCluster_HW.h  TPGStage2Emulation/*.hh
	g++ $(CPPFLAGSSTAGE2) $(LDFLAGS) $(CPPFLAGS) EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* test-stage2_Nov24/stage2CPCalc.cpp `root-config --libs --cflags` -L$(BOOST)/lib  -lboost_iostreams -lz -llzma -l yaml-cpp -o stage2CPCalc.exe

stage2TWCPEmulEMP.exe: test-stage2_Nov24/stage2TWCPEmulEMP.cpp inc/*.* EMPTools/CMSSWCode/DataFormats/L1THGCal/interface/HGCalCluster_HW.h  TPGStage2Emulation/*.hh
	g++ $(CPPFLAGSSTAGE2) $(LDFLAGS) $(CPPFLAGS) EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* test-stage2_Nov24/stage2TWCPEmulEMP.cpp `root-config --libs --cflags` -L$(BOOST)/lib  -lboost_iostreams -lz -llzma -l yaml-cpp -o stage2TWCPEmulEMP.exe

compareStage2TWCP_FWvsEmul.exe: test-stage2_Nov24/compareStage2TWCP_FWvsEmul.cpp EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/
	g++  $(CPPFLAGSSTAGE2) EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* test-stage2_Nov24/compareStage2TWCP_FWvsEmul.cpp `root-config --libs --cflags` -L$(BOOST)/lib  -lboost_iostreams -lz -llzma -o compareStage2TWCP_FWvsEmul.exe

writeStage2EMPRx_TWCP.exe: test-stage2_Nov24/writeStage2EMPRx_TWCP.cpp inc/*.*  TPGStage2Emulation/*.hh
	g++ $(CPPFLAGSSTAGE2) EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* test-stage2_Nov24/writeStage2EMPRx_TWCP.cpp `root-config --libs --cflags` -L$(BOOST)/lib  -lboost_iostreams -lz -llzma -o writeStage2EMPRx_TWCP.exe

writeStage2MultEMPRx_TWCP.exe: test-stage2_Nov24/writeStage2MultEMPRx_TWCP.cpp inc/*.*  TPGStage2Emulation/*.hh
	g++ $(CPPFLAGSSTAGE2) EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* test-stage2_Nov24/writeStage2MultEMPRx_TWCP.cpp `root-config --libs --cflags` -L$(BOOST)/lib  -lboost_iostreams -lz -llzma -o writeStage2MultEMPRx_TWCP.exe

readCPPyEmulOut.exe: test-stage2_Nov24/readCPPyEmulOut.cpp inc/*.*  TPGStage2Emulation/*.hh EMPTools/CMSSWCode/DataFormats/L1THGCal/interface/HGCalCluster_HW.h
	g++ $(CPPFLAGSSTAGE2) EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* test-stage2_Nov24/readCPPyEmulOut.cpp `root-config --libs --cflags` -L$(BOOST)/lib  -lboost_iostreams -lz -llzma -o readCPPyEmulOut.exe

compareCPPyCpp.exe: test-stage2_Nov24/compareCPPyCpp.cpp inc/*.*  TPGStage2Emulation/*.hh EMPTools/CMSSWCode/DataFormats/L1THGCal/interface/HGCalCluster_HW.h
	g++ $(LDFLAGS) $(CPPFLAGS) $(CPPFLAGSSTAGE2) EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* test-stage2_Nov24/compareCPPyCpp.cpp `root-config --libs --cflags` -L$(BOOST)/lib  -lboost_iostreams -lz -llzma -l yaml-cpp -o compareCPPyCpp.exe

testfp.exe:  test-stage2_Nov24/testfp.cpp
	g++ -IEMPTools/HLS_arbitrary_Precision_Types/include/ -I`root-config --incdir` test-stage2_Nov24/testfp.cpp `root-config --libs --cflags`  -o testfp.exe

compareStage2TWCP_FWvsEmul_allBxs.exe: test-stage2_Nov24/compareStage2TWCP_FWvsEmul_allBxs.cpp EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/
	g++  $(CPPFLAGSSTAGE2) EMPTools/CMSSWCode/L1Trigger/DemonstratorTools/src/* test-stage2_Nov24/compareStage2TWCP_FWvsEmul_allBxs.cpp `root-config --libs --cflags` -L$(BOOST)/lib  -lboost_iostreams -lz -llzma -o compareStage2TWCP_FWvsEmul_allBxs.exe

stage2SemiEmulationTCBits.exe: test-stage2_Nov24/stage2SemiEmulationTCBits.cpp inc/*.*  TPGStage2Emulation/*.hh
	g++ $(CPPFLAGSSTAGE2) test-stage2_Nov24/stage2SemiEmulationTCBits.cpp `root-config --libs --cflags` -o stage2SemiEmulationTCBits.exe -lEG $(LDFLAGS) $(CPPFLAGS) -l yaml-cpp

stage2SemiEmulationTCBitsAccumulatorFW.exe: test-stage2_Nov24/stage2SemiEmulationTCBitsAccumulatorFW.cpp inc/*.*  TPGStage2Emulation/*.hh
	g++ $(CPPFLAGSSTAGE2) test-stage2_Nov24/stage2SemiEmulationTCBitsAccumulatorFW.cpp `root-config --libs --cflags` -o stage2SemiEmulationTCBitsAccumulatorFW.exe -lEG $(LDFLAGS) $(CPPFLAGS) -l yaml-cpp

stage2SemiEmulator.exe: test-stage2_Nov24/stage2SemiEmulator.cpp inc/*.*  TPGStage2Emulation/*.hh
	g++ $(CPPFLAGSSTAGE2) test-stage2_Nov24/stage2SemiEmulator.cpp `root-config --libs --cflags` -o stage2SemiEmulator.exe -lEG $(LDFLAGS) $(CPPFLAGS) -l yaml-cpp

CalcResolution.exe: test-stage2_Nov24/CalcResolution.cpp inc/*.*  TPGStage2Emulation/*.hh
	g++ $(CPPFLAGSSTAGE2) test-stage2_Nov24/CalcResolution.cpp `root-config --libs --cflags` -o CalcResolution.exe -lEG $(LDFLAGS) $(CPPFLAGS) -l yaml-cpp

stage2BkgRate.exe: test-stage2_Nov24/stage2BkgRate.cpp inc/*.*  TPGStage2Emulation/*.hh
	g++ $(CPPFLAGSSTAGE2) test-stage2_Nov24/stage2BkgRate.cpp `root-config --libs --cflags` -o stage2BkgRate.exe -lEG $(LDFLAGS) $(CPPFLAGS) -l yaml-cpp

clean:
	rm *.exe
