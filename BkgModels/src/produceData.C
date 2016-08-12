#include "produceData.h"

int produceData( const char *infile , const char *opt , int nevents, int offset){
  
  Int_t npar = 11;
  Double_t par[11];
  std::vector<double> values;
  
  readData(infile,values);
  
  int mxdata = values.size();
  
  if( (mxdata%npar) != 0 ) {
    std::cout << "produceData> Error: input data is not compatible" << std::endl;
    exit(1);
  }
  else std::cout << "produceData> number of data parameters is: " << npar << std::endl;
  
  int nsets = mxdata / npar;
  int k(0);
  
  char filename[50];
  
  for( int i = offset; i < (nsets+offset); ++i) {
    
    sprintf(filename,"phis_study_%s_%d.root", opt , i );
    
    for( int j = 0;  j < npar; ++j) {
      par[j] = values[k];
      ++k;
    }
    
    std::cout << "produceData> producing data ..." << std::endl;
    
    RandomDists *exp = new RandomDists(filename , nevents);
    exp->initialise(par,npar);
    exp->RunExperiment();
    
    delete exp;
    
  }
  
  std::cout << "produceData> Data was successfully produced." << std::endl;

  return nsets;
    
}


int produceDataWC( const char *infile , const char *opt , int nevents, int offset){
  
  Int_t npar = 18;
  Double_t par[18];
  std::vector<double> values;
  
  readData(infile,values);
  
  int mxdata = values.size();
  
  if( (mxdata%npar) != 0 ) {
    std::cout << "produceDataWC> Error: input data is not compatible" << std::endl;
    exit(1);
  }
  else std::cout << "produceDataWC> number of data parameters is: " << npar << std::endl;
  
  int nsets = mxdata / npar;
  int k(0);
  
  char filename[50];
  
  for( int i = offset; i < (nsets+offset); ++i) {
    
    sprintf(filename,"phis_study_%s_%d.root", opt , i );
    
    for( int j = 0;  j < npar; ++j) {
      par[j] = values[k];
      ++k;
    }
    
    std::cout << "produceData> producing data with resolution ..." << std::endl;
    
    ConvolGenerator *exp = new ConvolGenerator(filename , nevents);
    exp->initialise(par,npar);
    exp->RunExperiment();
    
    delete exp;
    
  }
  
  std::cout << "produceData> Data was successfully produced." << std::endl;

  return nsets;
  
}
