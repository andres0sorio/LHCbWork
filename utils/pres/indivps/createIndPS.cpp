#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cmath>

#include <sys/stat.h>
#include <sys/unistd.h>

#include "ps_util2.hpp"

int main(int iargv, const char **argv) {
  
  const char *sourceName = argv[1];

  const char *histoName = argv[2];

  std::string output = std::string(histoName) + std::string(".ps");
  
  if(iargv > 1) {
    
    std::cout << "createIndPS> Selected histogram: " 
	      << histoName << std::endl
	      << "createIndPS> The PS file will be: " 
	      << output << std::endl;
    std::cout << "createIndPS> andres@hep.man.ac.uk" << std::endl;
    
  }
  else {
    std::cout << "usage : " << "createPS" << " source(dir) histogram" << std::endl;
    return 1;
  }
  
  std::string targetName = histoName + std::string(".tex");

  std::ofstream * outputFileName;
  
  outputFileName = new std::ofstream (targetName.c_str());
  
  ////////////////////////////////////////////
  // head of file
  std::ifstream is ("head_tex.dat");
  if(!is) {
    std::cout << "createIndPS> could not open input file> " << "head_tex.dat" << std::endl;
    exit(1);
  }
  
  char buffer[256];
  while (! is.eof())
    {
      is.getline (buffer,256);
      *outputFileName << buffer << std::endl; 
    }
  is.close();
  
  /////////////////////////////////////  
  chdir(sourceName);
  system("ls -p > list.txt");
  
  std::ifstream * cdir;
  cdir = new std::ifstream("list.txt", std::ifstream::in);
  
  findEPS ( cdir , outputFileName , output.c_str(), histoName );

  system("rm list.txt");

  chdir("../");

  //////////////////////////////////////
  // end of file
  
  *outputFileName << "\\end{document}" << std::endl;
  
  //////////////////////////////////////
  
  delete cdir;
  outputFileName->close();  
  delete outputFileName;
  
  std::string latexcmd = std::string("latex ") + targetName + std::string(" > createIndPS.log 2>&1");
  
  system(latexcmd.c_str());
  
  std::string dvipscmd = std::string("dvips -Ppdf  -t landscape ") + 
    histoName + std::string(".dvi >> createIndPS.log 2>&1");
  
  system(dvipscmd.c_str());
  
  system("rm *.aux *.dvi *.tex *.log");
  
  std::cout << "createIndPS> done. " << std::endl;
  
  return 0;
  
}
