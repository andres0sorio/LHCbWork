#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "Utilities.h"

int main(int iargv, const char **argv) {

  if(iargv < 3 ) {
    std::cout << "usage : genParFiles [prototype] [nfiles]" << std::endl;
    return 0;
  }
  
  
  const char *in = argv[1];
  int nf         = atoi(argv[2]);
  
  std::ofstream *os;
  std::vector<std::ofstream*> files;
  std::vector<std::ofstream*>::iterator itr;
  
  std::vector<double> v1;
  std::vector<double> v2;

  readData(in, v1, v2);
  
  int maxpars = v1.size();
  
  double val = 0.10;
  double dx  = 0.01;
  double init= val;
  
  for( int k = 0; k < nf; ++k) {
    
    val = init + dx * k;
    
    std::cout << "generating parameter file: " << k << std::endl;
    
    char filename[50];
    
    sprintf(filename,"initialSO_test_%d.dat",k);
    
    os = new std::ofstream(filename,ofstream::out);

    (*os) << "//initial values passed to Minuit" << std::endl;
    
    for ( int i = 0; i < maxpars; ++i) {
      
      if ( i != 2 ) (*os) << v1[i] << '\t' << (int)v2[i] << std::endl;
      else (*os) << val << '\t' << 1 << std::endl;
      
    }
    
    
    
    os->close();
  
  }
  
  
  for(itr = files.begin(); itr != files.end(); ++itr) {
    //(*itr)->close();
    delete (*itr);
  }
  
  return 1;
  
} 
