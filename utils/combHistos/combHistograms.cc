#include <TROOT.h>
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>
#include <sys/stat.h>
#include <sys/unistd.h>

#include "root_utilities.h"
#include "string_utilities.h"

bool isThereOption(const char *);

bool isThereOption(const char *arg) 
{
  bool ans = false;
  
  if (std::string(arg) == std::string("-sb") ) {
    std::cout << "combHistograms> option -sb (signal+background) selected. "
	      << std::endl;
    ans = true;
    
  }
  else ans = false;
  
  return ans;
  
}

int main(int iargv, const char **argv) {

  std::string maindir = std::string("./output/");
  const char *option = argv[1];
  const char *targetdirName = argv[2];

  int minarg = 4;
  int minparam = 2;
  bool withSB = false;


  /////////////////////////
    // aoptions set the logarithmic scale if =1
  int aoption(0);
  aoption = 0;
    
  if(iargv == 1) {
    std::cout << "usage: " << "combHistograms [option -sb] target(dir) source1 source2 ..." 
	      << std::endl;
    std::cout << "at least two sources are needed (one if -sb option is given)." 
	      << std::endl;
    return 1;
  }
  
  /////////////////////////
  //look for the -sb option   
  if (isThereOption(option)) {withSB=true; ++minparam;}
  else targetdirName = argv[1];
  /////////////////////////////////////////////////////////////////////////////////
  
  if(iargv >= minarg) {
    
    std::cout << "combHistograms> resulting histograms will be stored in: "
	      << maindir 
	      << targetdirName 
	      << std::endl;
    std::cout << "combHistograms> andres@hep.man.ac.uk" << std::endl;
    
  }
  
  else {
    std::cout << "usage : " << "combHistograms [option -sb] target(dir) source1 source2 ..." 
	      << std::endl;
    std::cout << "at least two sources are needed if there is not an option." << std::endl;
    return 1;
  }
  
  int totalfiles = iargv - minparam;
  TString source[totalfiles];
  
  TList *FileList;
  FileList = new TList();
  
  ////////////////////////////////////////////////////////////////////////////////
  
  const char *cwd;
  unsigned int size = 256;
  char buffer[size];
  cwd = getcwd(buffer, size);
  
  TString dir = TString(cwd)  + TString("/") ;
  
  ////////////////////////////////////////////////////////////
  
  for (int i = 0; i < totalfiles; i++ ){
    source[i]= dir + TString(argv[i+minparam]);
    FileList->Add(TFile::Open(source[i].Data()));
  }
  
  ////////////////////////////////////////////////////////////
  
  std::string dirname = getBaseName(targetdirName);
  
  ////////////////////////////////////////////////////////////
  std::string outdir = maindir + dirname;
  mkdir(outdir.c_str(),S_IRWXU);
  chdir(outdir.c_str());
  
  ////////////////////////////////////////////////////////////
  TStyle *st1 = new TStyle("st1", "wwsAnalysis style");
  setStyleOptions(st1);
  st1->cd(); 
  
  TFile *basefile = (TFile*)FileList->First();
  
  ////////////////////////////////////////////////////////////////////

  
  std::string options = std::string(cwd) + std::string("/HistogramOptions.dat");
  
  std::ifstream *hOptionFile = new std::ifstream(options.c_str());
  if(!hOptionFile->good()) {
    std::cout << "combHistograms> could not open Option file> " << std::endl;
    return 1;
  }
  
  

  if (withSB) combineSignalBack(basefile, FileList, aoption, hOptionFile);
  
  else combineHistograms(basefile, FileList);
  
  ////////////////////////////////////////////////////////////////////
  
  std::cout << "combHistograms> done. " << std::endl;
  
  delete st1;
  delete FileList;
  delete hOptionFile;
  return 0;
  
}
