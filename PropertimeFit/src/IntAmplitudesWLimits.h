// $Id: IntAmplitudesWLimits.h,v 1.1 2007/03/21 14:30:00 aosorio Exp $
#ifndef INTAMPLITUDESWLIMITS_H 
#define INTAMPLITUDESWLIMITS_H 1

// Include files

/** @class IntAmplitudesWLimits IntAmplitudesWLimits.h
 *  
 *
 *  @author Andres Felipe OSORIO OLIVEROS
 *  @date   2007-03-25
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

inline double A1def( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin, double tmax)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valB = 0.5*k0*((1.0+cosphis)*(tauL)*exp(-(1.0/tauL)*tmax)
                        +(1.0-cosphis)*(tauH)*exp(-(1.0/tauH)*tmax)
                        + 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)
                                                                         +Dms*cos(Dms*tmax))));
  
  double valA = 0.5*k0*((1.0+cosphis)*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cosphis)*(tauH)*exp(-(1.0/tauH)*tmin)
                        + 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
                                                                         +Dms*cos(Dms*tmin))));
  
  return (valA-valB);
  
}


inline double A2def( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin, double tmax)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valB = 0.5*k0*((1.0+cosphis)*(tauL)*exp(-(1.0/tauL)*tmax)
                        +(1.0-cosphis)*(tauH)*exp(-(1.0/tauH)*tmax)
                        + 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)
                                                                         +Dms*cos(Dms*tmax))));
  
  double valA = 0.5*k0*((1.0+cosphis)*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cosphis)*(tauH)*exp(-(1.0/tauH)*tmin)
                        + 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
                                                                         +Dms*cos(Dms*tmin))));
  
  return (valA-valB);
  
}

inline double A3def( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin, double tmax)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valB =  0.5*k0*((1.0-cosphis)*(tauL)*exp(-(1.0/tauL)*tmax)
                         +(1.0+cosphis)*(tauH)*exp(-(1.0/tauH)*tmax)
                         - 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)
                                                                          +Dms*cos(Dms*tmax))));
  
  double valA =  0.5*k0*((1.0-cosphis)*(tauL)*exp(-(1.0/tauL)*tmin)
                         +(1.0+cosphis)*(tauH)*exp(-(1.0/tauH)*tmin)
                         - 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
                                                                          +Dms*cos(Dms*tmin))));
  
  return (valA-valB);
  
}


inline double B1def( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin, double tmax)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valB = 0.5*k0*((1.0+cosphis)*(tauL)*exp(-(1.0/tauL)*tmax)
                        +(1.0-cosphis)*(tauH)*exp(-(1.0/tauH)*tmax)
                        - 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)
                                                                         +Dms*cos(Dms*tmax))));
  
  double valA = 0.5*k0*((1.0+cosphis)*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cosphis)*(tauH)*exp(-(1.0/tauH)*tmin)
                        - 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
                                                                         +Dms*cos(Dms*tmin))));
  
  return (valA-valB);
  
}


inline double B2def( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin, double tmax)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valB = 0.5*k0*((1.0+cosphis)*(tauL)*exp(-(1.0/tauL)*tmax)
                        +(1.0-cosphis)*(tauH)*exp(-(1.0/tauH)*tmax)
                        - 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)
                                                                         +Dms*cos(Dms*tmax))));
  
  double valA = 0.5*k0*((1.0+cosphis)*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cosphis)*(tauH)*exp(-(1.0/tauH)*tmin)
                        - 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
                                                                         +Dms*cos(Dms*tmin))));
  
  return (valA-valB);
  
}


inline double B3def( double k0, double tauL, double tauH, double tauBar,
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin, double tmax)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valB = 0.5*k0*((1.0-cosphis)*(tauL)*exp(-(1.0/tauL)*tmax)
                        +(1.0+cosphis)*(tauH)*exp(-(1.0/tauH)*tmax)
                        + 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*sin(Dms*tmax)
                                                                         +Dms*cos(Dms*tmax))));
  
  double valA = 0.5*k0*((1.0-cosphis)*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0+cosphis)*(tauH)*exp(-(1.0/tauH)*tmin)
                        + 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
                                                                         +Dms*cos(Dms*tmin))));
  
  return (valA-valB);
  
}


// + Interference terms

inline double A4def( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin, double tmax)
{
  
  double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valB = 0.5*k0*cos(tphase2-tphase1)* ( tauH*exp(-(1.0/tauH)*tmax)
                                               + tauL*exp(-(1.0/tauL)*tmax)
                                               - cosphis* ( tauH*exp(-(1.0/tauH)*tmax) 
                                                            - tauL*exp(-(1.0/tauL)*tmax) )
                                               + 2.0 * sinphis * exp (-(1.0/tauBar)*tmax) * gammaDms 
                                               * ( Dms* cos(Dms*tmax) + (1.0/tauBar)*sin(Dms*tmax)) );
  
  double valA = 0.5*k0*cos(tphase2-tphase1)* ( tauH*exp(-(1.0/tauH)*tmin)
                                               + tauL*exp(-(1.0/tauL)*tmin)
                                               - cosphis* ( tauH*exp(-(1.0/tauH)*tmin) 
                                                            - tauL*exp(-(1.0/tauL)*tmin) )
                                               + 2.0 * sinphis * exp (-(1.0/tauBar)*tmin) * gammaDms 
                                               * ( Dms* cos(Dms*tmin) + (1.0/tauBar)*sin(Dms*tmin)) );
  
  return (valA-valB);
  
}

inline double B4def( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin, double tmax)
{
  
  double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valB = 0.5*k0*cos(tphase2-tphase1)*( tauH*exp(-(1.0/tauH)*tmax)
                                              + tauL*exp(-(1.0/tauL)*tmax)
                                              - cosphis* ( tauH*exp(-(1.0/tauH)*tmax) 
                                                           - tauL*exp(-(1.0/tauL)*tmax) )
                                              - 2.0 * sinphis * exp (-(1.0/tauBar)*tmax) * gammaDms 
                                              * ( Dms* cos(Dms*tmax) + (1.0/tauBar)*sin(Dms*tmax)) );
  
  double valA = 0.5*k0*cos(tphase2-tphase1)*( tauH*exp(-(1.0/tauH)*tmin)
                                              + tauL*exp(-(1.0/tauL)*tmin)
                                              - cosphis* ( tauH*exp(-(1.0/tauH)*tmin) 
                                                           - tauL*exp(-(1.0/tauL)*tmin) )
                                              - 2.0 * sinphis * exp (-(1.0/tauBar)*tmin) * gammaDms 
                                              * ( Dms* cos(Dms*tmin) + (1.0/tauBar)*sin(Dms*tmin)) );
  
  return (valA-valB);
  
}


inline double A5def( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin, double tmax)
{
  
  double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valB = k0 * ( gammaDms*sin(tphase1)*exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*cos(Dms*tmax) 
                                                                      - Dms*sin(Dms*tmax))
                       - gammaDms*cos(tphase1)*cosphis*exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)
                                                                                *sin(Dms*tmax)
                                                                                + Dms*cos(Dms*tmax))
                       - 0.5*((tauH)*exp(-(1.0/tauH)*tmax)-(tauL)*exp(-(1.0/tauL)*tmax))
                       *cos(tphase1)*sinphis);
  
  double valA = k0 * ( gammaDms*sin(tphase1)*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*cos(Dms*tmin) 
                                                                      - Dms*sin(Dms*tmin))
                       - gammaDms*cos(tphase1)*cosphis*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)
                                                                                *sin(Dms*tmin)
                                                                                + Dms*cos(Dms*tmin))
                       - 0.5*((tauH)*exp(-(1.0/tauH)*tmin)-(tauL)*exp(-(1.0/tauL)*tmin))
                       *cos(tphase1)*sinphis);
  
  return (valA-valB);
  
}

inline double B5def( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin, double tmax)
{
  
  double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valB = - k0*(gammaDms*sin(tphase1)*exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)*cos(Dms*tmax)
                                                                     -Dms*sin(Dms*tmax))
                      - gammaDms*cos(tphase1)*cosphis*exp(-(1.0/tauBar)*tmax)*((1.0/tauBar)
                                                                               *sin(Dms*tmax)
                                                                               + Dms*cos(Dms*tmax))
                      + 0.5*((tauH)*exp(-(1.0/tauH)*tmax)-(tauL)*exp(-(1.0/tauL)*tmax))
                      *cos(tphase1)*sinphis);
  
  double valA = - k0*(gammaDms*sin(tphase1)*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*cos(Dms*tmin)
                                                                     -Dms*sin(Dms*tmin))
                      - gammaDms*cos(tphase1)*cosphis*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)
                                                                               *sin(Dms*tmin)
                                                                               + Dms*cos(Dms*tmin))
                      + 0.5*((tauH)*exp(-(1.0/tauH)*tmin)-(tauL)*exp(-(1.0/tauL)*tmin))
                      *cos(tphase1)*sinphis);
  
  return (valA-valB);
  
}

/* 
   
The following integrals are done in the range [tmin, inf [
Mathematica assited with the integration
Andres Osorio
Nov 25 2006 - CERN

*/


inline double A1lim( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valA = 0.5*k0*((1.0+cosphis)*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cosphis)*(tauH)*exp(-(1.0/tauH)*tmin)
                        + 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
                                                                         +Dms*cos(Dms*tmin))));
  
  return valA;
  
}


inline double A2lim( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valA = 0.5*k0*((1.0+cosphis)*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cosphis)*(tauH)*exp(-(1.0/tauH)*tmin)
                        + 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
                                                                         +Dms*cos(Dms*tmin))));
  
  return valA;
  
}

inline double A3lim( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valA =  0.5*k0*((1.0-cosphis)*(tauL)*exp(-(1.0/tauL)*tmin)
                         +(1.0+cosphis)*(tauH)*exp(-(1.0/tauH)*tmin)
                         - 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
                                                                          +Dms*cos(Dms*tmin))));
  
  return valA;
  
}


inline double B1lim( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valA = 0.5*k0*((1.0+cosphis)*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cosphis)*(tauH)*exp(-(1.0/tauH)*tmin)
                        - 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
                                                                         +Dms*cos(Dms*tmin))));
  
  return valA;
  
}


inline double B2lim( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valA = 0.5*k0*((1.0+cosphis)*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0-cosphis)*(tauH)*exp(-(1.0/tauH)*tmin)
                        - 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
                                                                         +Dms*cos(Dms*tmin))));
  
  return valA;
  
}


inline double B3lim( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin)
{
  
  double gammaDms = (1.0 / (((1.0/tauBar)*(1.0/tauBar)) + (Dms*Dms)) );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valA = 0.5*k0*((1.0-cosphis)*(tauL)*exp(-(1.0/tauL)*tmin)
                        +(1.0+cosphis)*(tauH)*exp(-(1.0/tauH)*tmin)
                        + 2.0*gammaDms*sinphis*(exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*sin(Dms*tmin)
                                                                         +Dms*cos(Dms*tmin))));
  
  return valA;
  
}


// + Interference terms

inline double A4lim( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin)
{
  
  double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valA = 0.5*k0*cos(tphase2-tphase1)* ( tauH*exp(-(1.0/tauH)*tmin)
                                               + tauL*exp(-(1.0/tauL)*tmin)
                                               - cosphis* ( tauH*exp(-(1.0/tauH)*tmin) 
                                                            - tauL*exp(-(1.0/tauL)*tmin) )
                                               + 2.0 * sinphis * exp (-(1.0/tauBar)*tmin) * gammaDms 
                                               * ( Dms* cos(Dms*tmin) + (1.0/tauBar)*sin(Dms*tmin)) );
  
  return valA;
  
}

inline double B4lim( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin)
{
  
  double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valA = 0.5*k0*cos(tphase2-tphase1)* ( tauH*exp(-(1.0/tauH)*tmin)
                                               + tauL*exp(-(1.0/tauL)*tmin)
                                               - cosphis* ( tauH*exp(-(1.0/tauH)*tmin) 
                                                            - tauL*exp(-(1.0/tauL)*tmin) )
                                               - 2.0 * sinphis * exp (-(1.0/tauBar)*tmin) * gammaDms 
                                               * ( Dms* cos(Dms*tmin) + (1.0/tauBar)*sin(Dms*tmin)) );
  
  return valA;
  
}

inline double A5lim( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin)
{
  
  double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valA = k0 * ( gammaDms*sin(tphase1)*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*cos(Dms*tmin) 
                                                                      - Dms*sin(Dms*tmin))
                       - gammaDms*cos(tphase1)*cosphis*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)
                                                                                *sin(Dms*tmin)
                                                                                + Dms*cos(Dms*tmin))
                       - 0.5*((tauH)*exp(-(1.0/tauH)*tmin)-(tauL)*exp(-(1.0/tauL)*tmin))
                       *cos(tphase1)*sinphis);
  
  return valA;
  
}

inline double B5lim( double k0, double tauL, double tauH, double tauBar, 
                     double Dms, double phis, double tphase1, double tphase2,
                     double tmin)
{
  
  double gammaDms = pow ( ((1.0/tauBar)*(1.0/tauBar))+(Dms*Dms), -1.0 );
  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double valA = - k0*(gammaDms*sin(tphase1)*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)*cos(Dms*tmin)
                                                                     -Dms*sin(Dms*tmin))
                      - gammaDms*cos(tphase1)*cosphis*exp(-(1.0/tauBar)*tmin)*((1.0/tauBar)
                                                                               *sin(Dms*tmin)
                                                                               + Dms*cos(Dms*tmin))
                      + 0.5*((tauH)*exp(-(1.0/tauH)*tmin)-(tauL)*exp(-(1.0/tauL)*tmin))
                      *cos(tphase1)*sinphis);
  
  return valA;
  
}


#endif // INTAMPLITUDESWLIMITS_H
