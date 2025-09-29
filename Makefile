#all: emul_test-beam_Sep23 findEMax read_econt_Jul24 GenerateEmpRxFile
LDFLAGS=-L$(HOME)/Software/yaml-cpp/lib64
CPPFLAGS=-I$(HOME)/Software/yaml-cpp/include -I common/inc -I offline/inc -I inc -I TPGStage1Emulation/ -I TPGFEEmulation/ -I`root-config --incdir`

# LDFLAGS=-L/usr/lib64
# CPPFLAGS=-I/usr/include/yaml-cpp -I common/inc -I offline/inc -I inc -I TPGStage1Emulation/ -I TPGFEEmulation/ -I`root-config --incdir` 

all:  TpgTimingCheck.exe dump_event25.exe readTPGdata.exe #sipm-emul.exe findEMax.exe #findFixedpattern.exe findEMax.exe tpgdata_3T_tcproc.exe  tpgdata_3T_fe.exe  scanadc_Sep24.exe dump_event.exe tpgdata_2T_fe.exe emul_Sep24.exe emul_3T_Sep24.exe dqm-tpgdata.exe dqm-emulSep24.exe #loop_emul_Sep24.exe scanconfigval_Sep24.exe emul_Sep24.exe validateFixedADC.exe # findEMax.exe GenerateEmpRxFile.exe dump_event.exe emul_test-beam_Sep23.exe 

emul_test-beam_Sep23.exe:  test-beam_Sep23_macros/emul_test-beam_Sep23.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep23_macros/emul_test-beam_Sep23.cpp  -l yaml-cpp `root-config --libs --cflags` -o emul_test-beam_Sep23.exe

findEMax.exe: test-beam_Sep23_macros/findEMax.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++  $(LDFLAGS) $(CPPFLAGS)  test-beam_Sep23_macros/findEMax.cpp  -l yaml-cpp `root-config --libs --cflags` -o findEMax.exe

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

dqm-tpgdata.exe: test-beam_Sep24_macros/dqm-tpgdata.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep24_macros/dqm-tpgdata.cpp  -l yaml-cpp `root-config --libs --cflags` -o dqm-tpgdata.exe -lm

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

dqm-emulSep24.exe: test-beam_Sep24_macros/dqm-emulSep24.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep24_macros/dqm-emulSep24.cpp  -l yaml-cpp `root-config --libs --cflags` -o dqm-emulSep24.exe -lm

loop_emul_Sep24.exe: test-beam_Sep24_macros/loop_emul_Sep24.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep24_macros/loop_emul_Sep24.cpp  -l yaml-cpp `root-config --libs --cflags` -o loop_emul_Sep24.exe -lm

validateFixedADC.exe: test-beam_Sep24_macros/validateFixedADC.C inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep24_macros/validateFixedADC.C  -l yaml-cpp `root-config --libs --cflags` -o validateFixedADC.exe -lm

dump_event.exe: test-beam_Sep24_macros/dump_event.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep24_macros/dump_event.cpp  -l yaml-cpp `root-config --libs --cflags` -o dump_event.exe

GenerateEmpRxFile.exe: TPGStage1Emulation/GenerateEmpRxFile.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS)  TPGStage1Emulation/GenerateEmpRxFile.cpp  -l yaml-cpp `root-config --libs --cflags` -o GenerateEmpRxFile.exe

TestUnpackerTCProcInterface.exe: TPGStage1Emulation/TestUnpackerTCProcInterface.cpp TPGFEEmulation/*.hh *.h TPGStage1Emulation/*.hh common/inc/*.h
	g++ -I TPGStage1Emulation/ -I. TPGStage1Emulation/HGCalLayer1PhiOrderFwImpl.cc -I. -I common/inc -I inc -I TPGFEEmulation/ TPGStage1Emulation/TestUnpackerTCProcInterface.cpp -L /opt/local/lib  -l yaml-cpp `root-config --libs --cflags` -I`root-config --incdir` -o TestUnpackerTCProcInterface.exe

sipm-emul.exe: test-beam_Sep24_macros/sipm-emul.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h
	g++ $(LDFLAGS) $(CPPFLAGS) test-beam_Sep24_macros/sipm-emul.cpp  -l yaml-cpp `root-config --libs --cflags` -o sipm-emul.exe -lm

dump_event25.exe: test-beam_Sep25_macros/dump_event.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h offline/inc/*.h 
	g++ -std=c++11 $(LDFLAGS) $(CPPFLAGS) test-beam_Sep25_macros/dump_event.cpp  -l yaml-cpp `root-config --libs --cflags` -o dump_event25.exe

TpgTimingCheck.exe: test-beam_Sep25_macros/TpgTimingCheck.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h offline/inc/*.h 
	g++ -std=c++11 $(LDFLAGS) $(CPPFLAGS) test-beam_Sep25_macros/TpgTimingCheck.cpp  -l yaml-cpp `root-config --libs --cflags` -o TpgTimingCheck.exe

readTPGdata.exe: test-beam_Sep25_macros/readTPGdata.cpp inc/*.*  TPGFEEmulation/*.hh TPGStage1Emulation/*.hh common/inc/*.h offline/inc/*.h 
	g++ -std=c++11 $(LDFLAGS) $(CPPFLAGS) test-beam_Sep25_macros/readTPGdata.cpp  -l yaml-cpp `root-config --libs --cflags` -o readTPGdata.exe

clean:
	rm *.exe
