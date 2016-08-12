// $Id: PDFsWRes.h,v 1.7 2007/03/22 13:05:33 aosorio Exp $
#ifndef PDFSWRES_H 
#define PDFSWRES_H 1

// Include files
#include "ResModels.h"
#include "PDFs.h"

/* @class PDFsWRes PDFsWRes.h
 *  
 *
 *  @author Andres Osorio
 *  @date   2007-02-08
 */

double convolveTotal ( double (*) (double (*) (double *, double *),
                                   double *, double *, double *),
                       double *, double *, double *);

double convolveSum   ( double (*) (double (*) (double *, double *),
                                   double *, double *, double *),
                       double *, double *, double *);

double convolveWp    ( double (*) (double (*) (double *, double *),
                                   double *, double *, double *),
                       double *, double *, double *);

double convolveWm    ( double (*) (double (*) (double *, double *),
                                   double *, double *, double *),
                       double *, double *, double *);

double totalpdfWRes     ( double *, double *);

double jpsiphipdfWRes   ( double *, double *);

double properTimeWRes   ( double *, double *);

double properTimeWRes2  ( double *, double *);


#endif // PDFSWRES_H
