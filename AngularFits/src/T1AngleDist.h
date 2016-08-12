// $Id: $
#ifndef T1ANGLEDIST_H 
#define T1ANGLEDIST_H 1

// Include files

#include "Parameters.h"

/** @class T1AngleDist T1AngleDist.h
 *  
 *
 *  @author Andres Osorio Oliveros
 *  @date   2006-10-21
 */

class T1AngleDist {

 public: 
  /// Standard constructor
  T1AngleDist( );
  
  T1AngleDist( const char * );  //import list of parameters from file
  
  virtual ~T1AngleDist( ); ///< Destructor
  
  void AddParameter( const char *, double );
  
  void SetParameter( const char *, double );

  void GetParameter( const char *, double &);

  double Evaluate(double, double);
  
  virtual double timePDF();
  
  virtual double thetaPDF();

  virtual double A0(double);
  
  virtual double Ap(double);
  
  virtual double At(double);

  
 protected:
  
  int maxpars;
  
  std::vector< Parameters* > parameters;

 private:

  
};
#endif // T1ANGLEDIST_H
