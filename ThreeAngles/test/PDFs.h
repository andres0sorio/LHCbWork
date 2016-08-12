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

double normfactorWp( double *par ); 
double normfactorWm( double *par ); 

double nfactorjpsiphi( double *par );

double jpsiphipdf(double *x, double *par);

double angpdfwres(double *x, double *par);

double properTimeConv(double *x, double *par);

/*
  integrated PDFs functions used to get the event variables
*/

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

///////////////////////////////////////////
struct pdf_params
{
  
  double gBar;
  double Dgamma;
  double Rt;
  double R0;
  double dt1;
  double dt2;
  double phis;
  double Dms;
  double omega;
  double accko;
  
  pdf_params() { };
  ~pdf_params() { };
  
  
};

//Derivative of psiPDF w.r.t. angle
//double thetaPDF( double , void *);

//Derivative of psiPDF w.r.t. angle
double DpsiWpPDF  ( double , void *);
double DpsiWmPDF  ( double , void *);

//Derivative of psiPDF w.r.t. angle
double DphiWpPDF  ( double , void *);
double DphiWmPDF  ( double , void *);


#endif
