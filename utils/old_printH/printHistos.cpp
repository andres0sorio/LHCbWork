#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cmath>

#include <sys/stat.h>
#include <sys/unistd.h>

#include "utilities.hpp"
#include "plot_util.hpp"

int main(int iargv, const char **argv) {
  
  if(iargv == 1) {
    std::cout << "usage : " << "printHistos" << " (source file)" << std::endl;
    exit(1);
  }

  const char *sourceName = argv[1];
  
  std::string output = getBaseName(sourceName);
  
  std::string outputDir = std::string("./") + output;
    
  if(iargv == 2) {
    std::cout << "printHistos> EPS Histograms will be stored at: " 
	      << std::string(argv[1]) << "/" << std::endl;
    std::cout << "printHistos> andres@hep.man.ac.uk" << std::endl;
  }
  else if(iargv > 2) {
    std::cout << "usage : " << "printHistos" << " (source file)" << std::endl;
    exit(1);
  }
  
  //////////////////////////////////////
  const char *cwd;
  unsigned int size = 100;
  char buffer[size];
  cwd = getcwd(buffer, size);
  
  TString dir = TString(cwd)  + TString("/") ;
  
  ////////////////////////////////////////////////////////////
  TString source;
  
  source = dir + TString(argv[1]);
  
  TFile *inputFile = new TFile(source);
  
  mkdir(outputDir.c_str(),S_IRWXU);
  chdir(outputDir.c_str());
  
  //////////////////////////////////////
  
  std::cout << "printHistos> reading from file: " << std::endl;
  std::cout << source << std::endl;
  
  printHistograms(inputFile);
  
  std::cout << "printHistos> done. " << std::endl;
  
  chdir("../");
  
  inputFile->Close();
  
  return 0;
  
}
