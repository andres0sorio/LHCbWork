// $Id: ResModels.h,v 1.5 2007/02/11 23:50:38 aosorio Exp $
#ifndef RESMODELS_H 
#define RESMODELS_H 1

// Include files
#include "ThreeAnglesCommon.h"

#include <TRandom.h>
#include <TRandom3.h>

#include "GSLHelpers.h"

/** @ ResModels ResModels.h
 *  
 *
 *  @author Andres Osorio Oliveros
 *  @date   2006-08-15
 */

double convolve ( double (*) (double *, double *), 
                  double (*) (double (*) (double *, double *),
                              double *, double *, double *),
                  double *, double *, double *);

double withGaussian (double , void *);

double withGaussian (double (*) (double *, double *), 
                     double *, double *, double *);

double with2Gaussians ( double (*) (double *, double *), 
                        double *, double *, double *);

double smearValue( double (*) (double, double *),
                   double, double *);

double withGaussian( double, double * );

#endif // RESMODELS_H
