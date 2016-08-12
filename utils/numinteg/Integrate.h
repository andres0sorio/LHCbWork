#ifndef INTEGRATE_H
#define INTEGRATE_H

//-----------------------------------------------------------------------------
// Integrate.C: numerical integration
// Using Taka Yasuda simpson()
// 2006-07-24 : Andres Osorio Oliveros
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

double integrate2d( Ptr2Function func , double * , double * , double * , double * );

double integrate3d( Ptr2Function func , double * , double * , double * , double * );

double simpson(std::vector<double>& , double, double );

double trapzd(std::vector<double>& , double , double );

#endif
