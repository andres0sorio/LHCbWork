// $Id: RandIntegrate.C,v 1.1 2006/08/31 15:31:36 aosorio Exp $
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

double projectTo1D( Ptr2Function  func, 
		    double * x, 
		    double *par, 
		    double *xlo , 
		    double *xhi ,
		    int np,
		    int nvar)
{
  
  /*
    projectTo1D: takes a function with nvar variables and projects to 1D (by integrating).
    Then remaining variable will always be considered x[0]
  */
  
  double r1(-1.0);
  double  v[5];
  double xx[5];
  
  for(int i=0; i< 5; ++i) {
    v[i] =-1.0;
    xx[i]=-1.0;
  }
   
  int seed = gRandom->GetSeed();
  gRandom->SetSeed(seed);
  
  for(int i=1; i < nvar; ++i){
    v[i]=xhi[i-1]-xlo[i-1];
  }
    
  double sum = 0.0;
  xx[0] = x[0]; // the one variable that doesnt get integrated
  
  for(int i=0; i < np; ++i) {
    
    for( int k=1; k < nvar; ++k) {
      r1 = gRandom->Uniform(1);
      xx[k] = xlo[k-1] + r1*v[k];
    }
    sum+= func(xx,par);
  }
  
  double vol = 1.0;
  for(int k=1; k < nvar; ++k) vol = v[k]*vol;
  
  return vol * (sum / (double) np);
  
}
