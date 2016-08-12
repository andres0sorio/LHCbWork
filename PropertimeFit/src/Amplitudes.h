// $Id: Amplitudes.h,v 1.4 2007/02/23 10:12:07 aosorio Exp $
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

//const  double  tphase1 = 0.0;
//const  double  tphase2= 3.1415926535897932385;


//Time dependent amplitudes

inline  double A1( double time, double k0, double tauL, double tauH, double tauBar, 
                   double Dms, double phis, double tphase1, double tphase2)
{

  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double val = 0.5*k0*((1.0+cosphis)*exp(-(1.0/tauL)*time)
                       +(1.0-cosphis)*exp(-(1.0/tauH)*time)
                       +2.0*exp(-(1.0/tauBar)*time)*sin(Dms*time)*sinphis);
  return val;
  
}

inline  double A2( double time, double k0, double tauL, double tauH, double tauBar, 
                   double Dms, double phis, double tphase1, double tphase2)
{

  double cosphis = cos(phis);
  double sinphis = sin(phis);

  double val = 0.5*k0*((1.0+cosphis)*exp(-(1.0/tauL)*time)
                       +(1.0-cosphis)*exp(-(1.0/tauH)*time)
                       +2.0*exp(-(1.0/tauBar)*time)*sin(Dms*time)*sinphis);
  
  return val;
  
}


inline  double A3( double time, double k0, double tauL, double tauH, double tauBar, 
                   double Dms, double phis, double tphase1, double tphase2)
{

  double cosphis = cos(phis);
  double sinphis = sin(phis);

  double val = 0.5*k0*((1.0-cosphis)*exp(-(1.0/tauL)*time)
                        +(1.0+cosphis)*exp(-(1.0/tauH)*time)
                        -2.0*exp(-(1.0/tauBar)*time)*sin(Dms*time)*sinphis); 
   
   return val;
  
}

inline  double B1( double time, double k0, double tauL, double tauH, double tauBar, 
                   double Dms, double phis, double tphase1, double tphase2)
{
  
  double cosphis = cos(phis);
  double sinphis = sin(phis);

  double val = 0.5*k0*((1.0+cosphis)*exp(-(1.0/tauL)*time)
                       +(1.0-cosphis)*exp(-(1.0/tauH)*time)
                       -2.0*exp(-(1.0/tauBar)*time)*sin(Dms*time)*sinphis);
   
   return val;
  
}

inline  double B2( double time, double k0, double tauL, double tauH, double tauBar, 
                   double Dms, double phis, double tphase1, double tphase2)
{

  double cosphis = cos(phis);
  double sinphis = sin(phis);

  double val = 0.5*k0*((1.0+cosphis)*exp(-(1.0/tauL)*time)
                       +(1.0-cosphis)*exp(-(1.0/tauH)*time)
                       -2.0*exp(-(1.0/tauBar)*time)*sin(Dms*time)*sinphis);
  
  return val;
  
}


inline double B3( double time, double k0, double tauL, double tauH, double tauBar, 
                  double Dms, double phis, double tphase1, double tphase2)
{

  double cosphis = cos(phis);
  double sinphis = sin(phis);

  double val = 0.5*k0*((1.0-cosphis)*exp(-(1.0/tauL)*time)
                       +(1.0+cosphis)*exp(-(1.0/tauH)*time)
                       +2.0*exp(-(1.0/tauBar)*time)*sin(Dms*time)*sinphis);
  
  return val;
  
}

// + Interference terms

//A4 = Re{ A_{0} A_{||} }
inline double A4( double time, double k0, double tauL, double tauH, double tauBar, 
                  double Dms, double phis, double tphase1, double tphase2)
{

  double cosphis = cos(phis);
  double sinphis = sin(phis);

  double val =  0.5*k0*cos(tphase2-tphase1)*((1.0+cosphis)*exp(-(1.0/tauL)*time)
                                             +(1.0-cosphis)*exp(-(1.0/tauH)*time)
                                             +2.0*exp(-(1.0/tauBar)*time)*sin(Dms*time)*sinphis);
  
  return val;
  
}

//A5 = Im{ A_{||} A_{T} }
inline double A5( double time, double k0, double tauL, double tauH, double tauBar, 
                  double Dms, double phis, double tphase1, double tphase2)
{

  double cosphis = cos(phis);
  double sinphis = sin(phis);

  double val = k0*(exp(-(1.0/tauBar)*time)*(sin(tphase1)*cos(Dms*time)
                                            -cos(tphase1)*sin(Dms*time)*cosphis)
                   -0.5*(exp(-(1.0/tauH)*time)-exp(-(1.0/tauL)*time))*cos(tphase1)*sinphis);
  
  return val;
  
}

//A6 = Im{ A_{0} A_{T} }  
inline double A6( double time, double k0, double tauL, double tauH, double tauBar, 
                  double Dms, double phis, double tphase1, double tphase2)
{
  
  double cosphis = cos(phis);
  double sinphis = sin(phis);

  double val = k0*(exp(-(1.0/tauBar)*time)*(sin(tphase2)*cos(Dms*time)
                                            -cos(tphase2)*sin(Dms*time)*cosphis)
                   -0.5*(exp(-(1.0/tauH)*time)-exp(-(1.0/tauL)*time))*cos(tphase2)*sinphis);
  return val;
  
}

inline double B4( double time, double k0, double tauL, double tauH, double tauBar, 
                  double Dms, double phis, double tphase1, double tphase2)
{

  double cosphis = cos(phis);
  double sinphis = sin(phis);
  
  double val = 0.5*k0*cos(tphase2-tphase1)*((1.0+cosphis)*exp(-(1.0/tauL)*time)
                                            +(1.0-cosphis)*exp(-(1.0/tauH)*time)
                                            -2.0*exp(-(1.0/tauBar)*time)*sin(Dms*time)*sinphis);
  
  return val;
  
}

inline double B5( double time, double k0, double tauL, double tauH, double tauBar, 
                  double Dms, double phis, double tphase1, double tphase2)
{

  double cosphis = cos(phis);
  double sinphis = sin(phis);

  double val = -k0*(exp(-(1.0/tauBar)*time)*(sin(tphase1)*cos(Dms*time)
                                             -cos(tphase1)*sin(Dms*time)*cosphis)
                    +0.5*(exp(-(1.0/tauH)*time)-exp(-(1.0/tauL)*time))*cos(tphase1)*sinphis);
  
  return val;
  
}


inline double B6( double time, double k0, double tauL, double tauH, double tauBar, 
                  double Dms, double phis, double tphase1, double tphase2)
{

  double cosphis = cos(phis);
  double sinphis = sin(phis);

  double val = -k0*(exp(-(1.0/tauBar)*time)*(sin(tphase2)*cos(Dms*time)
                                             -cos(tphase2)*sin(Dms*time)*cosphis)
                    +0.5*(exp(-(1.0/tauH)*time)-exp(-(1.0/tauL)*time))*cos(tphase2)*sinphis);
  
  return val;
  
}

#endif // AMPLITUDES_H
