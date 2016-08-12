// $Id: AngularTerms.h,v 1.2 2007/04/09 11:49:08 aosorio Exp $
#ifndef ANGULARTERMS_H 
#define ANGULARTERMS_H 1

// Include files

/** @class AngularTerms AngularTerms.h
 *  
 *
 *  @author Andres Osorio Oliveros
 *  @date   2006-12-08
 */

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <cmath>
#include <vector>

/*
  
Definition in the transversity base

*/

inline  double f1( double x, double y, double z)
{
  return 2.0 * cos(y)*cos(y)*( 1.0 - sin(x)*sin(x)*cos(z)*cos(z) );
}

inline  double f2( double x, double y, double z)
{
  return sin(y)*sin(y)*( 1.0 - sin(x)*sin(x)*sin(z)*sin(z) );
}

inline  double f3( double x, double y, double z)
{
  return sin(y)*sin(y)*sin(x)*sin(x);
}

inline  double f4( double x, double y, double z)
{
  return sin(y)*sin(y)*sin(2.0*x)*sin(z);
}

inline  double f5( double x, double y, double z)
{
  double invsqrtwo = pow ( sqrt(2.0) , -1.0);
  return -1.0*invsqrtwo*sin(2.0*y)*sin(x)*sin(x)*sin(2.0*z);
}

inline  double f6( double x, double y, double z)
{
  double invsqrtwo = pow ( sqrt(2.0) , -1.0);
  return invsqrtwo*sin(2.0*y)*sin(2.0*x)*cos(z);
    
}

/*

Definition in the helicity basis

*/


inline  double h1( double x, double y, double z)
{
  return 4.0* sin(x)*sin(x)*cos(y)*cos(y);
}

inline  double h2( double x, double y, double z)
{
  return ((1.0+(cos(x)*cos(x)))*sin(y)*sin(y))-sin(x)*sin(x)*sin(y)*sin(y)*cos(2.0*z);
}

inline  double h3( double x, double y, double z)
{
  return ((1.0+(cos(x)*cos(x)))*sin(y)*sin(y))+sin(x)*sin(x)*sin(y)*sin(y)*cos(2.0*z);
}

inline  double h4( double x, double y, double z)
{
  return 2.0*sin(x)*sin(x)*sin(y)*sin(y)*sin(2.0*z);
}

inline  double h5( double x, double y, double z)
{
  double sqrtwo = sqrt(2.0);
  return -1.0*sqrtwo*sin(2.0*x)*sin(2.0*y)*cos(z);
}

inline  double h6( double x, double y, double z)
{
  double sqrtwo = sqrt(2.0);
  return sqrtwo*sin(2.0*x)*sin(2.0*y)*sin(z);
}



#endif // ANGULARTERMS_H
