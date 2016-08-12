// $Id: PDFsBkg.h,v 1.2 2007/04/03 00:31:08 aosorio Exp $
#ifndef PDFSBKG_H 
#define PDFSBKG_H 1

// Include files

#include "ThreeAnglesCommon.h"
#include "Amplitudes.h"
#include "IntAmplitudes.h"
#include "AngularTerms.h"

/** @class PDFsBkg PDFsBkg.h
 * 
 *  Background pdfs implementation
 *
 *  @author Andres Osorio
 *  @date   2007-02-26
 */

double bkgPDF (double * , double * );

double expobkg( double * , double *);

double expobkgnorm( double * );

double bkgmassPDF( double *, double *);


#endif // PDFSBKG_H
