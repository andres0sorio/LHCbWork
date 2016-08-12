// $Id: PDFsBkg.C,v 1.6 2007/04/03 14:10:08 aosorio Exp $
// Include files 



// local
#include "PDFsBkg.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PDFsBkg
//
// 2007-02-26 : Andres Osorio
//-----------------------------------------------------------------------------


double bkgPDF( double *x, double *pars)
{
  double factor =  8.0*3.14159265358979312; //integration over the three angles
  double num = expobkg( x, pars );
  double den = expobkgnorm( pars ) * factor;

  return (num / den);
  
}


double expobkg( double *x, double *pars)
{
  
  double time = x[0];
  
  if ( time < 0 ) return 0.0;
  
  double tau  = pars[17];
  double val  = exp( -(1.0/tau)*time );
  
  return val;
  
}


double expobkgnorm( double * pars )
{
 
  double tmin   =  0.0;
  double tau    =  pars[17];
  double val    =  tau*exp(-(1.0/tau)*tmin);
    
  return val;
  
}

double bkgmassPDF( double *x , double * pars)
{
  
  double kappa = 1.0300;
  double min   = 5.2;
  double max   = 5.5;
  double norm  = (kappa)*exp(-(1.0/kappa)*min) 
    - (kappa)*exp(-(1.0/kappa)*max);
  double val   = (1.0/norm) * exp(-(1.0/kappa)*x[4]);
  return val;
  
}


