// $Id: PDFsWacc.h,v 1.1 2006/11/19 17:21:59 aosorio Exp $
#ifndef PDFSWACC_H 
#define PDFSWACC_H 1

// Include files
#include "Distributions.h"
#include "PDFs.h"

/** @class PDFsWacc PDFsWacc.h
 *  
 *
 *  @author Andres Felipe OSORIO OLIVEROS
 *  @date   2006-11-18
 */

double TimeAcceptance  ( double   , double * );

double pTimeAcc        ( double * , double * );

double pTimeAccPDF     ( double * , double * );

double jpsiphiAccPDF   ( double * , double * );

//....Normalisation factors......


double accnormfactor ( double * pars );

double pTimeAcc      ( double , void * );


#endif // PDFSWACC_H
