#ifndef PDFSTESTS_H
#define PDFSTEST_H

//-----------------------------------------------------------------------------
// PDFs.C: collection of pdfs
// 
// 2006-06-12 : Andres Osorio Oliveros
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

double jpsiphipdfA(double *x, double *par);

double jpsiphipdfB(double *x, double *par);

// decay rate as a function of phi_s
double decayratepdfw(double *x, double *par); 

#endif
