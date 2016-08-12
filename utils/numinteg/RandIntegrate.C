// $Id: RandIntegrate.C,v 1.1 2006/08/05 14:39:44 aooliver Exp $
// Include files 



// local
#include "RandIntegrate.h"

//-----------------------------------------------------------------------------
// Implementation file for: RandIntegrate
//
// 2006-08-05 : Andres Osorio Oliveros
//-----------------------------------------------------------------------------

//=============================================================================

double mcintegrate  ( Ptr2Function  func, double * x, double *par, double *xlo , double * xhi , int np)
{

  int seed = gRandom->GetSeed();
  gRandom->SetSeed(seed);

  double v = xhi[0]-xlo[0];
  double sum = 0.0;
  double xx[1];
  xx[0] = x[0];
  
  for(int i=0; i < np; ++i) {
    double r1 = gRandom->Uniform(1);
    xx[0] = xlo[0] + r1*(v);
    sum+= func(xx,par);
  }
  
  return v * (sum / (double) np);
  
}

double mcintegrate3d  ( Ptr2Function  func, double * x, double *par, double *xlo , double * xhi , int np)
{

  int seed = gRandom->GetSeed();
  gRandom->SetSeed(seed);

  double r1(0.);
  double r2(0.);
  double r3(0.);

  double v1 = xhi[0]-xlo[0];
  double v2 = xhi[1]-xlo[1];
  double v3 = xhi[2]-xlo[2];
  double sum = 0.0;
  double xx[3];
  xx[0] = x[0];
  xx[1] = x[1];
  xx[2] = x[2];
  
  for(int i=0; i < np; ++i) {
    r1 = gRandom->Uniform(1);
    r2 = gRandom->Uniform(1);
    r3 = gRandom->Uniform(1);
    xx[0] = xlo[0] + r1*v1;
    xx[1] = xlo[1] + r2*v2;
    xx[2] = xlo[2] + r3*v3;
    
    sum+= func(xx,par);
    
  }
  
  return v1 * v2 * v3 * (sum / (double) np);
  
}
