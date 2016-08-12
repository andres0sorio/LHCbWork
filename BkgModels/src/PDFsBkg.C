// $Id: $
// Include files 



// local
#include "PDFsBkg.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PDFsBkg
//
// 2007-02-26 : Andres Osorio
//-----------------------------------------------------------------------------

double convolveBkg(  double (*ptr2model) (double (*ptr2ff) (double *, double *),
                                         double *, double *, double *),
                    double *fpar, double *x, double *cvpar ) 
{
  
  double f = convolve( expobkg, ptr2model , fpar, x, cvpar);
  
  //Convolution of the d(t) - prompt is simply the Gaussian evaluated at time t
  double fp = 0.633;
  double fg = 0.801;
  double s1 = 0.0584;
  double mu1= 0.0022;
  double s2 = 0.471;
  double mu2= -0.1114;
  
  //fp = 0.0;
  
  
  double g = fg*TMath::Gaus( x[0] , mu1 , s1, 1) 
    + (1.0-fg)*TMath::Gaus( x[0] , mu2 , s2, 1); 
  
  return ( fp*g + (1.0-fp)*f );
  
}

double bkgPDF( double *x, double *pars)
{
  
  //Simple one - no resolution
  
  double num = expobkg( x, pars );
  double den = expobkgnorm( pars );
  
  return (num / den);
  
}


double expobkg( double *x, double *pars)
{
  
  double time = x[0];
  
  if ( time < 0 ) return 0.0;
  
  double tau  = 0.555;
  double val  = exp( -(1.0/tau)*time );
  
  return val;
  
}


double expobkgnorm( double * pars )
{
  
  double factor =  1.0; 
  double tmin   =  0.0;
  double tmax   = 20.0;
  double tau    =  0.555;
  double val    = tau*exp(-(1.0/tau)*tmin)
    - tau*exp(-(1.0/tau)*tmax);
  
  return factor*val;
  
}


