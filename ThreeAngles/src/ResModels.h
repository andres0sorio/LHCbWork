// $Id: ResModels.h,v 1.1.1.1 2006/08/22 14:07:20 aosorio Exp $
#ifndef RESMODELS_H 
#define RESMODELS_H 1

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
#include <TRandom.h>

#include "Integrate.h"


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

double withGaussian (double (*) (double *, double *), 
                     double *, double *, double *);

double with2Gaussians ( double (*) (double *, double *), 
			double *, double *, double *);

double smearValue( double (*) (double, double *),
		   double, double *);

double withGaussian( double, double * );

//double with2Gaussian2( double, double * );

#endif // RESMODELS_H
