#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include <sys/stat.h>
#include <sys/unistd.h>

#include "putLabels.h"

int main(int iargv, const char **argv) {
  
  if(iargv < 2 ) {
    std::cout << "usage : putLabels [source file]" << std::endl;
    return 1;
  }
  
  const char *in = argv[1];

  putLabels *source = new putLabels("histo.opts", in );
  
  source->processFile();
  
  std::cout << "putLabels> done. " << std::endl;

  delete source;
  
  return 0;
  
}
