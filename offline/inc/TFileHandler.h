#ifndef TFileHandler_HH
#define TFileHandler_HH

#include <iostream>
#include <string>

#include "TFile.h"

class TFileHandler {
public:
  TFileHandler() : _rootFile(nullptr) {
  }
  
  TFileHandler(const std::string &fileName, bool w=true) {
    initialise(fileName,w);
  }
  /*  
  TFileHandler(const std::string &fileName, bool w=true) : _fileName(fileName), _doWrite(w) {
    if(_fileName.size()<5 ||
       _fileName.substr(fileName.size()-5,5)!=".root") _fileName+=".root";

    std::cout << "TFileHandler:: Opening ROOT file " << _fileName << std::endl;
    _rootFile=new TFile(_fileName.c_str(),"RECREATE");
  }
  */
  void initialise(const std::string &fileName, bool w=true) {
    _fileName=fileName;
    _doWrite=w;
    
    if(_fileName.size()<5 ||
       _fileName.substr(fileName.size()-5,5)!=".root") _fileName+=".root";

    std::cout << "TFileHandler:: Opening ROOT file " << _fileName << std::endl;
    _rootFile=new TFile(_fileName.c_str(),"RECREATE");
  }

  
  ~TFileHandler() {
    if(_doWrite) {
      std::cout << "TFileHandler:: Writing ROOT file " << _fileName << std::endl;
      _rootFile->Write();
    }

    std::cout << "TFileHandler:: Closing ROOT file " << _fileName << std::endl;
    _rootFile->Close();
  
    delete _rootFile;
  }

  void cd() {
    _rootFile->cd();
  }

private:
  std::string _fileName;
  bool _doWrite;
  TFile *_rootFile;
};

#endif
