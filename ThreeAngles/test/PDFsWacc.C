// $Id: PDFsWacc.C,v 1.1 2006/11/24 22:14:57 aosorio Exp $
// Include files 
#include <Math/Integrator.h>
#include <Math/IntegrationTypes.h>
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/WrappedFunction.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

// local
#include "PDFsWacc.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PDFsWacc
//
// 2006-11-18 : Andres Felipe OSORIO OLIVEROS
//-----------------------------------------------------------------------------





//=============================================================================


double TimeAcceptance  ( double time , double *par )
{

  double KK  = 1.1;
  double arg = KK*time;
  return (pow( arg,3.0))/(1.0 + pow( arg,3.0));
  
}

double pTimeAcc     ( double * xx , double *par )
{
  
  double tt = xx[0];

  double value = TimeAcceptance(tt,par) * properTime (xx ,par);
  
  return value;
  
}

double pTimeAccPDF (  double * xx , double *par )
{
  
  double value = pTimeAcc( xx , par ) ;
  
  //Normalisation factor:
  
  double nf = accnormfactor( par );
  
  return (1.0 / nf) * value;
  
}


double jpsiphiAccPDF   ( double *xx, double *par )
{
  
  double tt = xx[0];
  
  double value = TimeAcceptance(tt,par) * jpsiphi(xx ,par);
  
  //Normalisation factor:
  //
  //
  
  return value;
  
}


//--------------------------------------------------
// Normalisation - auxiliary functions definition




//...............  Numerical integration of pTimeAcc


double accnormfactor( double * pars )
{

  pdf_params * pms = new pdf_params();

  pms->gBar   = pars[0];
  pms->Dgamma = pars[1];
  pms->Rt     = pars[2];
  pms->R0     = pars[3];
  pms->dt1    = 0.0;
  pms->dt2    = 0.0;
  pms->phis   = pars[6];
  pms->Dms    = pars[7];
  pms->omega  = pars[8];

  ROOT::Math::Integrator::GSLFuncPointer ff = &pTimeAcc;
  
  ROOT::Math::Integrator * nminteg = new ROOT::Math::Integrator(ff,ROOT::Math::Integration::ADAPTIVE, 1.E-9, 1E-6, 1000 );
  
  //--method for integration
  nminteg->SetIntegrationRule(ROOT::Math::Integration::GAUSS21);

  //--range [a,inf]
  double eval =  nminteg->IntegralUp(ff,pms,0.0);

  delete nminteg;
  delete pms;
  
  
  return eval;
  
}


double pTimeAcc (double x, void * params) {
  
  struct pdf_params *p 
    = (struct pdf_params *) params;
  
  double par[20];
  par[0]  = p->gBar;
  par[1]  = p->Dgamma;
  par[2]  = p->Rt;
  par[3]  = p->R0;
  par[4]  = 0.0; 
  par[5]  = 0.0;
  par[6]  = p->phis;
  par[7]  = p->Dms;
  par[8]  = p->omega;
  par[9]  = 1.0;
  par[10] = 1.0;
  par[11] = 1.0;
  
  double v1 = pTimeAcc( &x , par);
  
  return v1;
  
}


