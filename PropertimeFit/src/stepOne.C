#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "UMLlh_Fit.h"

int main(int iargv, const char **argv) {
  
  if(iargv < 4 ) {
    std::cout << "usage : stepOne [source file] [nevents] [fit option]" << std::endl;
    std::cout << "stepOne: 1-no resolution \n"
              << "         2-with resolution \n"
              << "         3-external errors (no resolution)\n"
              << "         4-external errors (resolution)\n";
    return 0;
  }
  
  const char *in = argv[1];

  int nevt  = atoi(argv[2]);
  
  int fopt  = atoi(argv[3]);
  
  stepOne(in, nevt, fopt);  
  
  return 1;
  
}
