#ifndef TEST_H
#define TEST_H

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

#include "Riostream.h"
#include <TROOT.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom.h>
#include <TH1D.h>

#include "Integrate.h"
#include "RandIntegrate.h"

//////////////////////


double afunction( double *, double *);

double afunction2( double *, double *);

void test();

void testMCI();

#endif
