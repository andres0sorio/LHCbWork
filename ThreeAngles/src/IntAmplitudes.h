// $Id: IntAmplitudes.h,v 1.2 2006/12/04 14:15:29 aosorio Exp $
#ifndef INTAMPLITUDES_H 
#define INTAMPLITUDES_H 1

// Include files

/** @class IntAmplitudes IntAmplitudes.h
 *  
 *
 *  @author Andres Felipe OSORIO OLIVEROS
 *  @date   2006-11-25
 */


#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <cmath>
#include <vector>


/* 
   
The following integrals are done in the range [tmin, tmax]
Mathematica assited with the integration
Andres Osorio Oliveros
Nov 25 2006 - CERN

*/

const double tmin = 0.0;
const double tmax = 20.0;


inline double A1def( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  
  double valB = 0.5*k0*((1.0+cos(phis))*(tauL)*exp(-(1.0/tauL)*tmax)
                        +(1.0-cos(phis))*(tauH)*exp(-(1.0/tauH)*tmax)
                        + 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)+Dms*cos(Dms*tmax))));
  
  double valA = 0.5*k0*((1.0+cos(phis))*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cos(phis))*(tauH)*exp(-(1.0/tauH)*tmin)
                        + 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)+Dms*cos(Dms*tmin))));
  
  return (valA-valB);
  
}


inline double A2def( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  
  double valB = 0.5*k0*((1.0+cos(phis))*(tauL)*exp(-(1.0/tauL)*tmax)
                        +(1.0-cos(phis))*(tauH)*exp(-(1.0/tauH)*tmax)
                        + 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)+Dms*cos(Dms*tmax))));
  
  double valA = 0.5*k0*((1.0+cos(phis))*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cos(phis))*(tauH)*exp(-(1.0/tauH)*tmin)
                        + 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)+Dms*cos(Dms*tmin))));
  
  return (valA-valB);
  
}

inline double A3def( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  
  double valB =  0.5*k0*((1.0-cos(phis))*(tauL)*exp(-(1.0/tauL)*tmax)
                         +(1.0+cos(phis))*(tauH)*exp(-(1.0/tauH)*tmax)
                         - 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)+Dms*cos(Dms*tmax))));
  
  double valA =  0.5*k0*((1.0-cos(phis))*(tauL)*exp(-(1.0/tauL)*tmin)
                         +(1.0+cos(phis))*(tauH)*exp(-(1.0/tauH)*tmin)
                         - 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)+Dms*cos(Dms*tmin))));
  
  return (valA-valB);
  
}


inline double B1def( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  
  double valB = 0.5*k0*((1.0+cos(phis))*(tauL)*exp(-(1.0/tauL)*tmax)
                        +(1.0-cos(phis))*(tauH)*exp(-(1.0/tauH)*tmax)
                        - 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)+Dms*cos(Dms*tmax))));
  
  double valA = 0.5*k0*((1.0+cos(phis))*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cos(phis))*(tauH)*exp(-(1.0/tauH)*tmin)
                        - 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)+Dms*cos(Dms*tmin))));
  
  return (valA-valB);
  
}


inline double B2def( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  
  double valB = 0.5*k0*((1.0+cos(phis))*(tauL)*exp(-(1.0/tauL)*tmax)
                        +(1.0-cos(phis))*(tauH)*exp(-(1.0/tauH)*tmax)
                        - 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)+Dms*cos(Dms*tmax))));
  
  double valA = 0.5*k0*((1.0+cos(phis))*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cos(phis))*(tauH)*exp(-(1.0/tauH)*tmin)
                        - 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)+Dms*cos(Dms*tmin))));
  
  return (valA-valB);
  
}


inline double B3def( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  
  double valB = 0.5*k0*((1.0-cos(phis))*(tauL)*exp(-(1.0/tauL)*tmax)
                        +(1.0+cos(phis))*(tauH)*exp(-(1.0/tauH)*tmax)
                        + 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)+Dms*cos(Dms*tmax))));
  
  double valA = 0.5*k0*((1.0-cos(phis))*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0+cos(phis))*(tauH)*exp(-(1.0/tauH)*tmin)
                        + 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)+Dms*cos(Dms*tmin))));
  
  return (valA-valB);
  
}


// + Interference terms

inline double A4def( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
  
  double valB = 0.5*k0*cos(tphase2-tphase1)* ( tauH*exp(-(1.0/tauH)*tmax)
                                               + tauL*exp(-(1.0/tauL)*tmax)
                                               - cos(phis)* ( tauH*exp(-(1.0/tauH)*tmax) - tauL*exp(-(1.0/tauL)*tmax) )
                                               + 2.0 * sin(phis) * exp (-(1.0/tauBar)*tmax) * gammaDms 
                                               * ( Dms* cos(Dms*tmax) + (1.0/tauBar)*sin(Dms*tmax)) );
  
  double valA = 0.5*k0*cos(tphase2-tphase1)* ( tauH*exp(-(1.0/tauH)*tmin)
                                               + tauL*exp(-(1.0/tauL)*tmin)
                                               - cos(phis)* ( tauH*exp(-(1.0/tauH)*tmin) - tauL*exp(-(1.0/tauL)*tmin) )
                                               + 2.0 * sin(phis) * exp (-(1.0/tauBar)*tmin) * gammaDms 
                                               * ( Dms* cos(Dms*tmin) + (1.0/tauBar)*sin(Dms*tmin)) );
  
  return (valA-valB);
  
}

inline double B4def( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
  
  double valB = 0.5*k0*cos(tphase2-tphase1)*( tauH*exp(-(1.0/tauH)*tmax)
                                    + tauL*exp(-(1.0/tauL)*tmax)
                                              - cos(phis)* ( tauH*exp(-(1.0/tauH)*tmax) - tauL*exp(-(1.0/tauL)*tmax) )
                                              - 2.0 * sin(phis) * exp (-(1.0/tauBar)*tmax) * gammaDms 
                                              * ( Dms* cos(Dms*tmax) + (1.0/tauBar)*sin(Dms*tmax)) );
  
  double valA = 0.5*k0*cos(tphase2-tphase1)*( tauH*exp(-(1.0/tauH)*tmin)
                                              + tauL*exp(-(1.0/tauL)*tmin)
                                              - cos(phis)* ( tauH*exp(-(1.0/tauH)*tmin) - tauL*exp(-(1.0/tauL)*tmin) )
                                              - 2.0 * sin(phis) * exp (-(1.0/tauBar)*tmin) * gammaDms 
                                              * ( Dms* cos(Dms*tmin) + (1.0/tauBar)*sin(Dms*tmin)) );
  
  return (valA-valB);
  
}


inline double A5def( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
  
  double valB = k0 * ( gammaDms*sin(tphase1)*exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*cos(Dms*tmax) - Dms*sin(Dms*tmax))
                       - gammaDms*cos(tphase1)*cos(phis)*exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)
                                                                                  + Dms*cos(Dms*tmax))
                       - 0.5*((tauH)*exp(-(1.0/tauH)*tmax)-(tauL)*exp(-(1.0/tauL)*tmax))*cos(tphase1)*sin(phis));
  
  double valA = k0 * ( gammaDms*sin(tphase1)*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*cos(Dms*tmin) - Dms*sin(Dms*tmin))
                       - gammaDms*cos(tphase1)*cos(phis)*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
                                                                                  + Dms*cos(Dms*tmin))
                       - 0.5*((tauH)*exp(-(1.0/tauH)*tmin)-(tauL)*exp(-(1.0/tauL)*tmin))*cos(tphase1)*sin(phis));
  
  return (valA-valB);
  
}

inline double B5def( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
  
  double valB = - k0*(gammaDms*sin(tphase1)*exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*cos(Dms*tmax)-Dms*sin(Dms*tmax))
                      - gammaDms*cos(tphase1)*cos(phis)*exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)
                                                                                 + Dms*cos(Dms*tmax))
                      + 0.5*((tauH)*exp(-(1.0/tauH)*tmax)-(tauL)*exp(-(1.0/tauL)*tmax))*cos(tphase1)*sin(phis));
  
  double valA = - k0*(gammaDms*sin(tphase1)*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*cos(Dms*tmin)-Dms*sin(Dms*tmin))
                      - gammaDms*cos(tphase1)*cos(phis)*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
                                                                                 + Dms*cos(Dms*tmin))
                      + 0.5*((tauH)*exp(-(1.0/tauH)*tmin)-(tauL)*exp(-(1.0/tauL)*tmin))*cos(tphase1)*sin(phis));
  
  return (valA-valB);
  
}

/* 
   
    The following integrals are done in the range [tmin, inf [
    Mathematica assited with the integration
    Andres Osorio
    Nov 25 2006 - CERN

*/


inline double A1( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  
  double valA = 0.5*k0*((1.0+cos(phis))*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cos(phis))*(tauH)*exp(-(1.0/tauH)*tmin)
                        + 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)+Dms*cos(Dms*tmin))));
  
  return valA;
  
}


inline double A2( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  
  double valA = 0.5*k0*((1.0+cos(phis))*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cos(phis))*(tauH)*exp(-(1.0/tauH)*tmin)
                        + 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)+Dms*cos(Dms*tmin))));
  
  return valA;
  
}

inline double A3( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  
  double valA =  0.5*k0*((1.0-cos(phis))*(tauL)*exp(-(1.0/tauL)*tmin)
                         +(1.0+cos(phis))*(tauH)*exp(-(1.0/tauH)*tmin)
                         - 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)+Dms*cos(Dms*tmin))));
  
  return valA;
  
}


inline double B1( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  
  double valA = 0.5*k0*((1.0+cos(phis))*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cos(phis))*(tauH)*exp(-(1.0/tauH)*tmin)
                        - 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)+Dms*cos(Dms*tmin))));
  
  return valA;
  
}


inline double B2( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  
  double valA = 0.5*k0*((1.0+cos(phis))*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cos(phis))*(tauH)*exp(-(1.0/tauH)*tmin)
                        - 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)+Dms*cos(Dms*tmin))));
  
  return valA;
  
}


inline double B3( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  
  double valA = 0.5*k0*((1.0-cos(phis))*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0+cos(phis))*(tauH)*exp(-(1.0/tauH)*tmin)
                        + 2.0*gammaDms*sin(phis)*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)+Dms*cos(Dms*tmin))));
  
  return valA;
  
}


// + Interference terms

inline double A4( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
  
  double valA = 0.5*k0*cos(tphase2-tphase1)* ( tauH*exp(-(1.0/tauH)*tmin)
                                     + tauL*exp(-(1.0/tauL)*tmin)
                                               - cos(phis)* ( tauH*exp(-(1.0/tauH)*tmin) - tauL*exp(-(1.0/tauL)*tmin) )
                                               + 2.0 * sin(phis) * exp (-(1.0/tauBar)*tmin) * gammaDms 
                                               * ( Dms* cos(Dms*tmin) + (1.0/tauBar)*sin(Dms*tmin)) );
  
  return valA;
  
}

inline double B4( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
  
  double valA = 0.5*k0*cos(tphase2-tphase1)* ( tauH*exp(-(1.0/tauH)*tmin)
                                               + tauL*exp(-(1.0/tauL)*tmin)
                                               - cos(phis)* ( tauH*exp(-(1.0/tauH)*tmin) - tauL*exp(-(1.0/tauL)*tmin) )
                                               - 2.0 * sin(phis) * exp (-(1.0/tauBar)*tmin) * gammaDms 
                                               * ( Dms* cos(Dms*tmin) + (1.0/tauBar)*sin(Dms*tmin)) );
  
  return valA;
  
}

inline double A5( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
  
  double valA = k0 * ( gammaDms*sin(tphase1)*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*cos(Dms*tmin) - Dms*sin(Dms*tmin))
                       - gammaDms*cos(tphase1)*cos(phis)*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
                                                                                  + Dms*cos(Dms*tmin))
                       - 0.5*((tauH)*exp(-(1.0/tauH)*tmin)-(tauL)*exp(-(1.0/tauL)*tmin))*cos(tphase1)*sin(phis));
  
  return valA;
  
}

inline double B5( double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
  
  double valA = - k0*(gammaDms*sin(tphase1)*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*cos(Dms*tmin)-Dms*sin(Dms*tmin))
                      - gammaDms*cos(tphase1)*cos(phis)*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
                                                                                 + Dms*cos(Dms*tmin))
                      + 0.5*((tauH)*exp(-(1.0/tauH)*tmin)-(tauL)*exp(-(1.0/tauL)*tmin))*cos(tphase1)*sin(phis));
  
  return valA;
  
}



#endif // INTAMPLITUDES_H
