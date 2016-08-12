#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "SimdetLoader.h"

int main(int iargv, const char **argv) {

  int arguments(0);

  arguments = iargv;
  
  if (arguments < 3 || arguments > 3) {
    std::cout << "usage: doLeptonicAnalysis [file] [nevents]" << std::endl;
    exit(1);
  }
  
  const char *fileName = argv[1];
  const int nEvent = atoi(argv[2]);
  
  std::cout << "wwsAnalysis> andres@hep.man.ac.uk" << std::endl;
  std::cout << "wwsAnalysis> reading file: " << fileName << std::endl;
  std::cout << "wwsAnalysis> events: " << nEvent << std::endl;
  
  //////////////////////////////////////////////////
  //open input file
  SimdetLoader *sdl = new SimdetLoader(fileName);
  
  int i(0);
  
  while(i < nEvent) {
    
    sdl->next_event();
    
    ////
    int max = sdl->single_event->gpv.size();
    
    if( fmod(double(i+1),double(2500)) == 0.0 ) 
      std::cout << "max: " << max << " *" << std::endl;
    
    sdl->close_event();
    
    ++i;
    
  }
  
  delete sdl;
  
  if ( i == nEvent ) {
    std::cout << "simdetStruct> Total events: " << i << " achieved." << std::endl;
    std::cout << "simdetStruct> terminated. " << std::endl;
  }
  else {std::cout << "simdetStruct> terminated. " << std::endl;}
  
  return 0;
  
}
