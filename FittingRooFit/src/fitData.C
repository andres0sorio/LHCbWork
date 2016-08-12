#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "OneAngleFit.h"

int main(int iargv, const char **argv) {
  
  if(iargv < 4 ) {
    std::cout << "usage : fitData [source file] [param file] [nevents]" 
	      << std::endl;
    return 0;
  }
  
  const char *infile   = argv[1];
  
  const char *param    = argv[2];
  
  int nevt = atoi(argv[3]);
  
  fitData(infile, param, nevt);  
  
  return 1;
  
}
