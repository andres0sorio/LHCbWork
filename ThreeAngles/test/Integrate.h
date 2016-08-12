#ifndef INTEGRATE_H
#define INTEGRATE_H

//-----------------------------------------------------------------------------
// Integrate.C: numerical integration
// Using Taka Yasuda simpson() implementation
// 2006-07-24 : Andres Osorio Oliveros
// 2006-08-02 : Added project1D and integrate3d: A. Osorio
//-----------------------------------------------------------------------------

// Include files
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <cmath>
#include <vector>

typedef double (*Ptr2Function) (double *, double *);

double integrate( Ptr2Function , double *, double *, double, double );

double simpson(std::vector<double>& , double, double );

double projectTo1D( Ptr2Function , double * , double * , double * , double * );

double integrate3d( Ptr2Function , double * , double * , double * , double * );

#endif
