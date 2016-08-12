// $Id: ThreeAnglesCommon.h,v 1.2 2007/02/14 00:07:26 aosorio Exp $
#ifndef THREEANGLESCOMMON_H 
#define THREEANGLESCOMMON_H 1

/** ThreeAnglesCommon.h : common header files
 *  
 *
 *  @author Andres Osorio
 *  @date   2007-02-08
 */

// Include files

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <cmath>
#include <vector>

#include "Riostream.h"
#include <TROOT.h>
#include <TMath.h>


//..... inline functions definitions needed ....

inline double K0(double r1, double r2) {
  //return (1.0-r1)*(1.0-r2);
  //return 1.0 / (1.0 + r1 + r2);
  //return 1.0 - r1 - r2;
  return r2;
}

inline double Kt(double r1, double r2)
{
  //return r1;
  //return r1 / (1.0 + r1 + r2);
  return r1;
}

inline double Kp(double r1, double r2)
{
  //return (1.0-r1)*r2;
  //return r2 / (1.0 + r1 + r2);
  return 1.0 - r1 - r2;
}

inline double GHeavy( double gbar, double dgamma )
{
  return gbar - 0.5*dgamma; //default is -
}

inline double GLight( double gbar, double dgamma )
{
  return gbar + 0.5*dgamma; //default is +
}



#endif // THREEANGLESCOMMON_H
