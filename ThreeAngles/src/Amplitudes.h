// $Id: Amplitudes.h,v 1.7 2006/11/27 14:26:02 aosorio Exp $
// Include files
#ifndef AMPLITUDES_H
#define AMPLITUDES_H 1


#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <cmath>
#include <vector>

const  double  tphase1 = 0.0;
const  double  tphase2= 3.1415926535897932385;


//Time dependent amplitudes

inline  double A1( double time, double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
   double val = 0.5*k0*((1.0+cos(phis))*exp(-(1.0/tauL)*time)
                       +(1.0-cos(phis))*exp(-(1.0/tauH)*time)
                       +2.0*exp(-(1.0/tauBar)*time)*sin(Dms*time)*sin(phis));
  return val;
  
}

inline  double A2( double time, double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{

   double val = 0.5*k0*((1.0+cos(phis))*exp(-(1.0/tauL)*time)
                       +(1.0-cos(phis))*exp(-(1.0/tauH)*time)
                       +2.0*exp(-(1.0/tauBar)*time)*sin(Dms*time)*sin(phis));
  
  return val;
  
}


inline  double A3( double time, double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
   double val = 0.5*k0*((1.0-cos(phis))*exp(-(1.0/tauL)*time)
                       +(1.0+cos(phis))*exp(-(1.0/tauH)*time)
                       -2.0*exp(-(1.0/tauBar)*time)*sin(Dms*time)*sin(phis)); 

  return val;
  
}

inline  double B1( double time, double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
   double val = 0.5*k0*((1.0+cos(phis))*exp(-(1.0/tauL)*time)
                       +(1.0-cos(phis))*exp(-(1.0/tauH)*time)
                       -2.0*exp(-(1.0/tauBar)*time)*sin(Dms*time)*sin(phis));
  
  return val;
  
}

inline  double B2( double time, double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{

   double val = 0.5*k0*((1.0+cos(phis))*exp(-(1.0/tauL)*time)
                       +(1.0-cos(phis))*exp(-(1.0/tauH)*time)
                       -2.0*exp(-(1.0/tauBar)*time)*sin(Dms*time)*sin(phis));
  
  return val;
  
}


inline double B3( double time, double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double val = 0.5*k0*((1.0-cos(phis))*exp(-(1.0/tauL)*time)
                       +(1.0+cos(phis))*exp(-(1.0/tauH)*time)
                       +2.0*exp(-(1.0/tauBar)*time)*sin(Dms*time)*sin(phis));
  
  return val;
  
}

// + Interference terms

inline double A4( double time, double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{

  double val =  0.5*k0*cos(tphase2-tphase1)*((1.0+cos(phis))*exp(-(1.0/tauL)*time)
                                             +(1.0-cos(phis))*exp(-(1.0/tauH)*time)
                                             +2.0*exp(-(1.0/tauBar)*time)*sin(Dms*time)*sin(phis));
  
  return val;
  
}

inline double A5( double time, double k0, double tauL, double tauH, double tauBar, double Dms, double phis)\
{
  
  double val = k0*(exp(-(1.0/tauBar)*time)*(sin(tphase1)*cos(Dms*time)
                                    -cos(tphase1)*sin(Dms*time)*cos(phis))
                   -0.5*(exp(-(1.0/tauH)*time)-exp(-(1.0/tauL)*time))*cos(tphase1)*sin(phis));
  
  return val;
  
}

inline double A6( double time, double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double val = k0*(exp(-(1.0/tauBar)*time)*(sin(tphase2)*cos(Dms*time)
                                    -cos(tphase2)*sin(Dms*time)*cos(phis))
                   -0.5*(exp(-(1.0/tauH)*time)-exp(-(1.0/tauL)*time))*cos(tphase2)*sin(phis));
  return val;
  
}

inline double B4( double time, double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double val = 0.5*k0*cos(tphase2-tphase1)*((1.0+cos(phis))*exp(-(1.0/tauL)*time)
                                            +(1.0-cos(phis))*exp(-(1.0/tauH)*time)
                                            -2.0*exp(-(1.0/tauBar)*time)*sin(Dms*time)*sin(phis));
  
  return val;
  
}

inline double B5( double time, double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double val = -k0*(exp(-(1.0/tauBar)*time)*(sin(tphase1)*cos(Dms*time)
                                     -cos(tphase1)*sin(Dms*time)*cos(phis))
                    +0.5*(exp(-(1.0/tauH)*time)-exp(-(1.0/tauL)*time))*cos(tphase1)*sin(phis));
  
  return val;
  
}


inline double B6( double time, double k0, double tauL, double tauH, double tauBar, double Dms, double phis)
{
  
  double val = -k0*(exp(-(1.0/tauBar)*time)*(sin(tphase2)*cos(Dms*time)
                                     -cos(tphase2)*sin(Dms*time)*cos(phis))
                    +0.5*(exp(-(1.0/tauH)*time)-exp(-(1.0/tauL)*time))*cos(tphase2)*sin(phis));
  
  return val;
  
}



#endif // AMPLITUDES_H
