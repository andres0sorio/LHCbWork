#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cmath>

#include <sys/stat.h>
#include <sys/unistd.h>

#include "ps_util.h"

int main(int iargv, const char **argv) {
  
  if (iargv == 1) {
    std::cout << "usage : " << "createPS" << " source(dir) " << std::endl;
    exit(1);
  }
  
  const char *sourceName = argv[1];
  std::string output = extractFileName(sourceName);

  if(iargv == 2) {
    std::cout << "createPS> The PS file will be: " 
	      << output << ".ps" <<std::endl;
    std::cout << "createPS> andres@hep.man.ac.uk" << std::endl;
  }
  else if(iargv > 2){
    std::cout << "usage : " << "createPS" << " source(dir) " << std::endl;
    exit(1);
  }
  
  std::string targetName = output + std::string(".tex");
  
  std::ofstream * outputFileName;
  
  outputFileName = new std::ofstream (targetName.c_str());


  ////////////////////////////////////////////
  // head of file
  std::ifstream is ("head_tex.dat");
  if(!is) {
    std::cout << "createPS> could not open input file> " << "head_tex.dat" << std::endl;
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
  
  int a = addEPS ( cdir , outputFileName , output.c_str() );
  
  updateTeXFile( outputFileName , a );
    
  system("rm list.txt");

  chdir("../");

  //////////////////////////////////////
  // end of file
  
  *outputFileName << "\\end{document}" << std::endl;
  
  //////////////////////////////////////
  
  delete cdir;
  outputFileName->close();  
  delete outputFileName;
  
  std::string latexcmd = std::string("latex ") + targetName + std::string(" > createPS.log 2>&1");
  
  std::cout << "createPS> making dvi ... " << std::endl;


#ifndef _DEBUG
  system(latexcmd.c_str());
  
  std::string dvipscmd = std::string("dvips -Ppdf -t a4 ") + 
    output + std::string(".dvi >> createPS.log 2>&1");
  
  std::cout << "createPS> making ps ... " << std::endl;
  
  system(dvipscmd.c_str());
  
  system("rm *.aux *.dvi *.tex *.log");
#endif

  std::cout << "createPS> done. " << std::endl;
  
  return 0;
  
}
