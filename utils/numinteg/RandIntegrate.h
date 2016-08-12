// $Id: RandIntegrate.h,v 1.1 2006/08/05 14:39:44 aooliver Exp $
#ifndef NUMINTEG_RANDINTEGRATE_H 
#define NUMINTEG_RANDINTEGRATE_H 1

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

/** @ RandIntegrate RandIntegrate.h numinteg/RandIntegrate.h
 *  Monte Carlo Integration implementation
 *
 *  @author Andres Osorio Oliveros
 *  @date   2006-08-05
 */

typedef double (*Ptr2Function) (double *, double *);


double mcintegrate  ( Ptr2Function  , double * , double * , double * , double * , int );

//double mcintegrate2d( Ptr2Function  , double * , double * , double * , double * );

double mcintegrate3d( Ptr2Function  , double * , double * , double * , double * , int );



#endif // NUMINTEG_RANDINTEGRATE_H
