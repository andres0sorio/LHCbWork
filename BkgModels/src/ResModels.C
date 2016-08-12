// $Id: ResModels.C,v 1.8 2007/02/23 13:09:50 aosorio Exp $
// Include files 


// local
#include "ResModels.h"

//-----------------------------------------------------------------------------
// Implementation file for : ResModels
//
// 2006-08-15 : Andres Osorio Oliveros
//-----------------------------------------------------------------------------

double convolve ( double (*ptr2ff) (double *, double *), 
                  double (*ptr2model) (double (*ptr2ff) (double *, double *),
                                       double *, double *, double *),
                  double *fpar, double *x, double *par)
{
  
  /*
    convolve(...) 
    general function that interfaces with the different
    models for convolution
  */
  
  double result = (*ptr2model) ( ptr2ff , fpar , x , par );
  
  return result;
  
}

double withGaussian( double (*f) (double *, double *), 
                     double *fp, double *x, double *p )
{
  
  double sc       = 6.0;                      //convolution extends to +- gaussian sigma
  double result   = 0.0;
  double invsq2pi = 0.3989422804014;        //(2 pi)^(-1/2)
  double xmin     = 0.0;
  double xmax     = 0.0;
  double sigma    = p[0];
  
  res_params * pms = new res_params();
  
  pms->gBar    = fp[0];
  pms->Dgamma  = fp[1];
  pms->Rt      = fp[2];
  pms->Rp      = fp[3];
  pms->dp1     = fp[4];
  pms->dp2     = fp[5];
  pms->phis    = fp[6];
  pms->Dms     = fp[7];
  pms->omega   = fp[8];
  pms->sigma   = sigma;
  pms->mu      = 0.0; // <- has to be zero since by definition (t-ti-mu)^2
  pms->time    = x[0]; //time at wich the convolution is evaluated
  pms->theta   = x[1];
  pms->psi     = x[2];
  pms->phi     = x[3];
  pms->qf      = x[4];
  pms->accko   = 0.0;
  pms->ptr2ff  = f;   //pointer to the pdf
  
  // Range of the convolution integral
  xmin = x[0] - sc*sigma;
  xmax = x[0] + sc*sigma;
  
  ROOT::Math::Integrator::GSLFuncPointer ff = &withGaussian;
  
  ROOT::Math::Integrator * nminteg = new ROOT::Math::Integrator(ff,
                                                                ROOT::Math::Integration::ADAPTIVE, 
                                                                1.E-10, 1.E-7, 1000 );
  //--method for integration
  nminteg->SetIntegrationRule(ROOT::Math::Integration::GAUSS21);
  
  //--range [xmin,xmax]
  result = (invsq2pi / sigma) * nminteg->Integral(ff,pms,xmin,xmax);
  
  delete nminteg;
  delete pms;
  
  return result;
  
}

double with2Gaussians( double (*f) (double *, double *), 
                       double *fp, double *x, double *p )
{
  
  double sc        = 6.0;
  double invsq2pi  = 0.3989422804014;        //(2 pi)^(-1/2)
  double result1   = 0.0;
  double result2   = 0.0;
  double xmin      = 0.0;
  double xmax      = 0.0;
  
  double sigma     = p[0];
  double f1        = p[1];
  double mu1       = p[2];
  double s1        = p[3];
  double mu2       = p[4];
  double s2        = p[5];

  res_params * pms = new res_params();
  pms->gBar    = fp[0];
  pms->Dgamma  = fp[1];
  pms->Rt      = fp[2];
  pms->Rp      = fp[3];
  pms->dp1     = fp[4];
  pms->dp2     = fp[5];
  pms->phis    = fp[6];
  pms->Dms     = fp[7];
  pms->omega   = fp[8];
  pms->time    = x[0]; //time at wich the convolution is evaluated
  pms->theta   = x[1];
  pms->psi     = x[2];
  pms->phi     = x[3];
  pms->qf      = x[4]; //tag
  pms->accko   = 0.0;
  pms->sigma   = 1.0;
  pms->mu      = 0.0;
  
  pms->ptr2ff  = f;   //pointer to the pdf
  
  ROOT::Math::Integrator::GSLFuncPointer ff = &withGaussian;

  //simulate lifetime resolution (it is done on a event by event basis)
  //sigma = smearValue ( withGaussian , sigma, p );
  
  sigma = 1.0;
  
  //Range of the convolution integral
  xmin = (x[0] + mu1*sigma) - sc*s1*sigma;
  xmax = (x[0] + mu1*sigma) + sc*s1*sigma;

  pms->sigma   = s1*sigma;
  pms->mu      = mu1*sigma;
  
  ROOT::Math::Integrator * nminteg = new ROOT::Math::Integrator(ff,
                                                                ROOT::Math::Integration::ADAPTIVE, 
                                                                1.E-10, 1.E-7, 1000 );
  //--method for integration
  nminteg->SetIntegrationRule(ROOT::Math::Integration::GAUSS21);
  
  //--range [xmin,xmax]
  result1 = (invsq2pi/s1) * nminteg->Integral(ff,pms,xmin,xmax); 
  
  delete nminteg;
  
  xmin = (x[0] + mu2*sigma) - sc*s2*sigma;
  xmax = (x[0] + mu2*sigma) + sc*s2*sigma;
  
  pms->sigma   = s2*sigma;
  pms->mu      = mu2*sigma;
  
  nminteg = new ROOT::Math::Integrator(ff,
                                       ROOT::Math::Integration::ADAPTIVE, 
                                       1.E-9, 1E-6, 1000 );
  
  result2 = (invsq2pi/s2) * nminteg->Integral(ff,pms,xmin,xmax); 
  
  delete nminteg;
  delete pms;
  
  return (f1* result1 + (1-f1)*result2);
  
}

double withGaussian (double xx, void *pars)
{
  
  //recover first the parameters
  struct res_params *p 
    = (struct res_params *) pars;
  
  double par[20];
  par[0]  = p->gBar;
  par[1]  = p->Dgamma;
  par[2]  = p->Rt;
  par[3]  = p->Rp;
  par[4]  = p->dp1;
  par[5]  = p->dp2;
  par[6]  = p->phis;
  par[7]  = p->Dms;
  par[8]  = p->omega;
  
  double sigma = p->sigma;
  double mu    = p->mu;
  double ti    = p->time;
  double tt[5];
  
  tt[0] = xx;
  tt[1] = p->theta;
  tt[2] = p->psi;
  tt[3] = p->phi;
  tt[4] = p->qf;
  
  //convolution:
  double val = (*p->ptr2ff)(tt,par) * TMath::Gaus(xx - ti, mu , sigma);
  
  return val;
  
}

double smearValue( double (*ptr2meth) (double, double *),
		   double x , double *par)
{
  
  double result(0.);
  result = (*ptr2meth)( x , par );
  return result;
  
}

double withGaussian( double x, double *p )
{

  //smear value - fast implementation to account detector effects
  double mean = x;
  double sigma = p[0];
  TRandom3 * rnd = new TRandom3(0);
  double r1 = rnd->Gaus(mean,sigma);
  delete rnd;
  return r1;
  
}
