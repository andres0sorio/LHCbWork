// $Id: Distributions.h,v 1.1 2006/11/24 22:14:56 aosorio Exp $
#ifndef DISTRIBUTIONS_H 
#define DISTRIBUTIONS_H 1

// Include files
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

/** @class Distributions Distributions.h
 *  
 *
 *  @author Andres Felipe OSORIO OLIVEROS
 *  @date   2006-11-19
 */

//Local include

#include "Amplitudes.h"
#include "ResModels.h"

double jpsiphi(double *x, double *par);

double properTime(double *x, double *par);

//Not needed for the moment
//double theta(double *x, double *par);
//double psi(double *x, double *par);
//double phi(double *x, double *par);







#endif // DISTRIBUTIONS_H
