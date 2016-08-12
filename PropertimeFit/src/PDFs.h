#ifndef PDFS_H
#define PDFS_H

//-----------------------------------------------------------------------------
// PDFs.C: collection of main pdfs for data generation and fitting
//         without resolution
//
// 2006-10-18 : Andres Osorio Oliveros
//-----------------------------------------------------------------------------

// Include files
#include "ThreeAnglesCommon.h"
#include "Amplitudes.h"
#include "IntAmplitudes.h"
#include "IntAmplitudesWLimits.h"
#include "AngularTerms.h"
#include "PDFsBkg.h"

double totalpdf      (double *, double *);

double bsmassPDF     (double *, double *);

double sumofsignals  (double *, double *);

double normfactorWp  (double * ); 

double normfactorWm  (double * ); 

double nfactorjpsiphi(double * );

double jpsiphipdf    (double *, double *);

double jpsiphiWpPDF  (double *, double *);

double jpsiphiWmPDF  (double *, double *);

double properTimePDF( double *, double *);

double properTimeWpPDF(double *, double *);

double properTimeWmPDF(double *, double *);

//////////////////////
//
double properTimePDF2( double *, double *);
double nfactorjpsiphiWL(double * );

#endif
