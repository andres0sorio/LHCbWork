// $Id: FitResults.h,v 1.3 2006/11/27 14:26:02 aosorio Exp $
#ifndef FITRESULTS_H 
#define FITRESULTS_H 1

// Include files
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

/** @class FitResults FitResults.h
 *  
 *
 *  @author Andres Felipe OSORIO OLIVEROS
 *  @date   2006-11-21
 */

struct  Parameters {
  
  // io functions
  friend std::istream& operator>>(std::istream &istr, Parameters &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, Parameters &rhs);

  int status;
  int maxparams;
  std::vector<double> param;
  std::vector<double> errparam;
  double atmin;
  bool has_succeded;
  
};

struct  ErrMatrix {
  
  // io functions
  friend std::istream& operator>>(std::istream &istr, ErrMatrix &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, ErrMatrix &rhs);
  int status;
  int maxparams;
  double matrix[10][10];
  bool has_errmatrix;
  
};


class FitResults {
public: 
  /// Standard constructor
  FitResults( ); 
  
  FitResults(FitResults &e) {
    vparams= e.vparams;
    verrmat= e.verrmat;
  }
  
  virtual ~FitResults( ); ///< Destructor
  
  void AddResults( const char * );
  
  std::vector<Parameters*>::iterator FirstParameters() { return vparams.begin();}
  
  std::vector<Parameters*>::iterator LastParameters() { return vparams.end(); }
  
  std::vector<ErrMatrix*>::iterator FirstErrMatrix() { return verrmat.begin();}
  
  std::vector<ErrMatrix*>::iterator LastErrMatrix() { return verrmat.end(); }
  
  int Size() { return vparams.size();}

  int NSucceded() { return n_succeded; }
  
protected:
  
private:
  
  std::ifstream * is;
  
  std::vector<Parameters*> vparams;
  
  std::vector<ErrMatrix*> verrmat;

  int n_succeded;
  
  int n_errmatrix;
  
};
#endif // FITRESULTS_H
