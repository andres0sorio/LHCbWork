#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cmath>

#include <sys/stat.h>
#include <sys/unistd.h>

#include "html_util.h"

int main(int iargv, const char **argv) {
  
  const char *sourceName = argv[1];

  std::string output = extractFileName(sourceName);
  
  if(iargv > 1) {
    
    std::cout << "createHTML> The HTML file will be: " 
	      << output << ".html" <<std::endl;
    std::cout << "createHTML> andres@hep.man.ac.uk" << std::endl;
    
  }
  else {
    std::cout << "usage : " << "createHTML" << " source(dir) " << std::endl;
    return 1;
  }
  
  std::string targetName = output + std::string(".html");
  
  std::ofstream * outputFileName;
  
  outputFileName = new std::ofstream (targetName.c_str());


  ////////////////////////////////////////////
  // head of file
  
  std::ifstream is ("head_html.dat");
  if(!is) {
    std::cout << "createHTML> could not open input file> " << "head_html.dat" << std::endl;
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
  // html body

  *outputFileName << "<b>" << output << "<br></b>" << std::endl; 

  chdir(sourceName);
  system("ls -p > list.txt");
  
  std::ifstream * cdir;
  cdir = new std::ifstream("list.txt", std::ifstream::in);

  addPictures ( cdir , outputFileName , output.c_str() );
  
  //////////////////////////////////////
  // end of file
  
  *outputFileName << "</body>"<< std::endl; 
  *outputFileName << "</html>" << std::endl; 
  *outputFileName << "<br>" << std::endl;
  
  //////////////////////////////////////
  
  std::cout << "createHTML> done. " << std::endl;
  
  delete cdir;
  
  system("rm list.txt");
  
  outputFileName->close();
  
  delete outputFileName;
  
  return 0;
  
}
