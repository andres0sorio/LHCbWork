// $Id: AngularTerms.h,v 1.1 2006/12/09 18:19:18 aosorio Exp $
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
  return invsqrtwo*sin(2.0*y)*sin(x)*sin(x)*sin(2.0*z);
}

inline  double f6( double x, double y, double z)
{
  double invsqrtwo = pow ( sqrt(2.0) , -1.0);
  return invsqrtwo*sin(2.0*y)*sin(2.0*x)*cos(z);
    
}



#endif // ANGULARTERMS_H
