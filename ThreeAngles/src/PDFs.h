#ifndef PDFS_H
#define PDFS_H

//-----------------------------------------------------------------------------
// PDFs.C: collection of main pdfs for data generation
// 
// 2006-10-18 : Andres Osorio Oliveros
//-----------------------------------------------------------------------------

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
#include <TFile.h>
#include <TMath.h>

#include "ResModels.h"

#include "Amplitudes.h"
#include "IntAmplitudes.h"

#include "AngularTerms.h"
//#include "IntAngularTerms.h"

double normfactorWp( double *par ); 

double normfactorWm( double *par ); 

double nfactorjpsiphi( double *par );

double jpsiphipdf(double *x, double *par);

double jpsiphiWpPDF(double *x, double *par);

double jpsiphiWmPDF(double *x, double *par);

double angpdfwres(double *x, double *par);

double properTimeConv(double *x, double *par);



//Bs
double properTimeWpPDF(double *x, double *par);
double thetaWpPDF(double *x, double *par);
double psiWpPDF(double *x, double *par);
double phiWpPDF(double *x, double *par);

//Bs(bar)
double properTimeWmPDF(double *x, double *par);
double thetaWmPDF(double *x, double *par);
double psiWmPDF(double *x, double *par);
double phiWmPDF(double *x, double *par);


//Derivative of thetaPDF w.r.t. angle
double DthetaWpPDF( double x, void * params);
double DthetaWmPDF( double x, void * params);


//Derivative of psiPDF w.r.t. angle
double DpsiWpPDF  ( double , void *);
double DpsiWmPDF  ( double , void *);

//Derivative of phiPDF w.r.t. angle
double DphiWpPDF  ( double , void *);
double DphiWmPDF  ( double , void *);


//functions definitions
inline double K0(double r1, double r2) 
{
  return (1.0-r1)*(1.0-r2);
}

inline double Kt(double r1, double r2)
{
  return r1;
}

inline double Kp(double r1, double r2)
{
  return (1.0-r1)*r2;
}

inline double GHeavy( double gbar, double deltag )
{
  return gbar*(1.0 + 0.5*deltag);
}

inline double GLight( double gbar, double deltag )
{
  return gbar*(1.0 - 0.5*deltag);
}

///////////////////////////////////////////
struct pdf_params
{
  
  double gBar;
  double Dgamma;
  double ft;
  double fp;
  double dt1;
  double dt2;
  double phis;
  double Dms;
  double omega;
  double accko;
  double time;
  
  pdf_params() { };
  ~pdf_params() { };
  
  
};

#endif
