// $Id: GSLHelpers.h,v 1.5 2007/04/05 08:53:46 aosorio Exp $
#ifndef GSLHELPERS_H 
#define GSLHELPERS_H 1

// Include files
#include <Math/Integrator.h>
#include <Math/IntegrationTypes.h>
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/WrappedFunction.h"
#include "Math/RootFinder.h"
#include "Math/RootFinderAlgorithms.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>


///////////////////////////////////////////
struct pdf_params
{
  
  double gBar;
  double Dgamma;
  double Rt;
  double Rp;
  double dp1;
  double dp2;
  double phis;
  double Dms;
  double omega;
  double accko;
  double time;
  double qf;
  double fsig;
  double tbkg;
    
  pdf_params() { };
  ~pdf_params() { };
  
  
};


///////////////////////////////////////////
struct res_params
{
  
  double gBar;
  double Dgamma;
  double Rt;
  double Rp;
  double dp1;
  double dp2;
  double phis;
  double Dms;
  double omega;
  double accko;
  double sigma;
  double mu;
  double time;
  double theta;
  double psi;
  double phi;
  double mass;
  double qf;
  double norm;
  double fsig;
  double tbkg;
  double tmin;
  double tmax;
  
    
  double (*ptr2ff) (double *, double *);
  
  res_params() { };
  ~res_params() { };
  
  
};

#endif // GSLHELPERS_H
