#include <string>

#include "produceData.h"

int main(int iargv, const char **argv) {
  
  if(iargv < 4 ) {
    std::cout << "usage : generate [input file] [opt] [nevents] [offset]" << std::endl;
    std::cout << "input file: set of input parameters" << std::endl;
    std::cout << "opt: an option" << std::endl;
    std::cout << "nevents: how many events to generate for each set" << std::endl;
    std::cout << "offset(optional): first job index" << std::endl;
    return 0;
  }
  
  const char *in  = argv[1];
  
  const char *opt = argv[2];
  
  int nevt        = atoi(argv[3]);
  
  int offset      = 0;
  
  if ( iargv == 5 ) {
    offset      = atoi(argv[4]);
  }
  else offset = 0;
  
  int maxdata = produceData(in,opt,nevt,offset);
  
  return maxdata;
  
}
