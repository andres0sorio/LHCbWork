// $Id: $
#ifndef PDFSBKG_H 
#define PDFSBKG_H 1

// Include files

#include "ThreeAnglesCommon.h"
#include "ResModels.h"

/** @class PDFsBkg PDFsBkg.h
 * 
 *  Background pdfs implementation
 *
 *  @author Andres Osorio
 *  @date   2007-02-26
 */

double convolveBkg    ( double (*) (double (*) (double *, double *),
                                    double *, double *, double *),
                        double *, double *, double *);

double bkgPDF (double * , double * );

double expobkg( double * , double *);

double expobkgnorm( double * );


#endif // PDFSBKG_H
