#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "UMLlh_Fit.h"

int main(int iargv, const char **argv) {
  
  if(iargv < 3 ) {
    std::cout << "usage : stepOne [source file] [nevents]" << std::endl;
    return 0;
  }
  
  const char *in = argv[1];

  int nevt = atoi(argv[2]);
  
  stepOne(in, nevt);  

  return 1;
  
}
