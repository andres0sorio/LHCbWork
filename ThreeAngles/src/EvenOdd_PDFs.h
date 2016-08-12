#ifndef EVENODD_PDFS_H
#define EVENODD_PDFS_H

//-----------------------------------------------------------------------------
// EvenOdd_PDFs.C: collection of even and odd pdfs 
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

#include "PDFs.h"

double properTimePDFOdd(double *x, double *par);

double properTimePDFEven(double *x, double *par);

double thetaEven(double *x, double *par);

double thetaOdd(double *x, double *par);

double psiPDFEven(double *x, double *par);

double psiPDFOdd(double *x, double *par);

double phiPDFEven(double *x, double *par);

double phiPDFOdd(double *x, double *par);

double SinphiEven(double *x, double *par);

double SinpsiEven(double *x, double *par);

double SinphiOdd(double *x, double *par);

double SinpsiOdd(double *x, double *par);

double CosthetaTrEven(double *x, double *par);

double CosthetaTrOdd(double  *x, double *par);

double CosthetaTrTotal(double  *x, double *par);


#endif
