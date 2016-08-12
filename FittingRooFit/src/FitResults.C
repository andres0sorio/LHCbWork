// $Id: FitResults.C,v 1.5 2006/11/27 14:26:02 aosorio Exp $
// Include files 



// local
#include "FitResults.h"

//-----------------------------------------------------------------------------
// Implementation file for class : FitResults
//
// 2006-11-21 : Andres Felipe OSORIO OLIVEROS
//-----------------------------------------------------------------------------

std::istream& operator>>(std::istream &istr, Parameters &rhs) {
  
  double fitvalue(0.);
  double fiterror(0.);
  
  istr >> rhs.status >> rhs.maxparams;
  
  if ( rhs.status == 4 ) { 
    rhs.has_succeded = false;
    return istr;
  }
  
  rhs.has_succeded = true;
  
  for (int i = 0; i < rhs.maxparams; ++i) {
    istr >> fitvalue >> fiterror;
    rhs.param.push_back(fitvalue);
    rhs.errparam.push_back(fiterror);
  }
  
  std::string amin;
  
  istr >> amin;
  istr >> rhs.atmin;
  
  return istr;
}

std::ostream& operator<<(std::ostream &ostr, Parameters &rhs) {

  int status    = rhs.status;

  if ( rhs.status != 0 ) {
    std::cout << "MIGRAD returned: " << status << std::endl;
    return ostr;
  }

  int maxparams = rhs.maxparams;
  double fitvalue(0.);
  double fiterror(0.);
  
  for (int i = 0; i < maxparams; ++i) {
    fitvalue = rhs.param[i];
    fiterror = rhs.errparam[i];
    ostr << fitvalue << '\t' << fiterror << '\n';
  }
  
  ostr << "Atmin: " << rhs.atmin << std::endl;
  
  return ostr;
}

// Error matrix
std::istream& operator>>(std::istream &istr, ErrMatrix &rhs) {
  
  istr >> rhs.status >> rhs.maxparams;
  
  if ( rhs.status == 4 || rhs.maxparams == 1 )  {
    rhs.has_errmatrix = false;
    return istr;
  }

  rhs.has_errmatrix = true;
  
  for (int i = 0; i < rhs.maxparams; ++i) {
    for ( int j = 0; j < rhs.maxparams; ++j) {
      istr >> rhs.matrix[i][j];
    }
  }
  
  return istr;
}

std::ostream& operator<<(std::ostream &ostr, ErrMatrix &rhs) {
  
  std::cout << rhs.status << " " << rhs.maxparams << std::endl;
  
  if ( rhs.maxparams == 1 || rhs.status != 0 )  return ostr;
  
  for (int i = 0; i < rhs.maxparams; ++i) {
    for ( int j = 0; j < rhs.maxparams; ++j) {
      double val = rhs.matrix[i][j];
      ostr << val << '\t';
    }
    ostr << '\n';
  }
  
  return ostr;
}

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
FitResults::FitResults(  ) {

  n_succeded = 0;
  n_errmatrix = 0;
  
}

//=============================================================================
// Destructor
//=============================================================================
FitResults::~FitResults() {
  
  std::vector<Parameters*>::iterator itr1 = vparams.begin();
  while(itr1 != vparams.end()) {
    delete *itr1;
    itr1++;
  } 
  std::vector<ErrMatrix*>::iterator itr2 = verrmat.begin();
  while(itr2 != verrmat.end()) {
    delete *itr2;
    itr2++;
  }
}
//=============================================================================

void FitResults::AddResults( const char * infile ) 
{
  
  is = new std::ifstream(infile);
  if(!is) {
    std::cout << "addResults> could not open input file> " << std::string(infile) << std::endl;
    return ;
  }
  
  //std::cout << "addResults> input file open" << std::endl;
  
  Parameters *data = new Parameters();
  ErrMatrix  *errm = new ErrMatrix();
  
  char header[256];
  is->getline(header,256);
  
  (*is) >> *data; //get data
  (*is) >> *errm;
  
  vparams.push_back(data);
  verrmat.push_back(errm);
  
  if ( data->has_succeded  ) ++n_succeded;
  if ( errm->has_errmatrix ) ++n_errmatrix;
  
  //if ( ! data->has_succeded  ) {
  //  std::cout << *data;
  //  std::cout << *errm;
  // }
  
  is->close();
  delete is;
  
  return;

}
