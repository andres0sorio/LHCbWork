#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>

#include <sys/stat.h>
#include <sys/unistd.h>

#include "root_utilities.h"
#include "string_utilities.h"

int main(int iargv, const char **argv) {
  
  const char *targetfileName = argv[1];

  if(iargv > 3) {

    print_message("addHistograms> resulting histograms will be stored in: ");
    print_message(targetfileName);
    print_message("addHistograms> andres@hep.man.ac.uk");
  
  }
  
  else {
    print_message("usage : addHistograms target source1 source2 ...");
    print_message("at least two sources are needed.");
    return 1;
  }
  
  int totalfiles = iargv - 2;
  TString source[totalfiles];
  
  TList *FileList;
  TFile *Target;

  const char *cwd;
  unsigned int size = 200;
  char buffer[size];
  cwd = getcwd(buffer, size);
  
  TString outputfile = TString(cwd)
    + TString("/rootfiles/")
    + TString(targetfileName)
    + TString(".root");
  
  Target = new TFile(outputfile, "RECREATE");
  
  FileList = new TList();
  
  for (int i = 0; i < totalfiles; i++ ){
    source[i]= TString(argv[i+2]);
    FileList->Add(TFile::Open(source[i].Data()));
  }
  
  addRootHistos(Target, FileList);
  
  std::string maindir = std::string("./output/");
  std::string filename = getBaseName(targetfileName);
  
  std::string outdir = maindir + filename;
  mkdir(outdir.c_str(),S_IRWXU);
  chdir(outdir.c_str());
  
  //printHistograms(Target);
  
  print_message("addHistograms> done. ");
   
  return 0;
}
