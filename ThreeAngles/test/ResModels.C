// $Id: ResModels.C,v 1.1 2006/10/18 09:45:19 aosorio Exp $
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
    convolve(...) - general function(method) that interfaces with the different
    models for convolution
  */
  
  double result = (*ptr2model) ( ptr2ff , fpar , x , par );
  
  return result;
  
}

double withGaussian( double (*f) (double *, double *), 
                     double *fp, double *x, double *p )
{
  
  double ti       = x[0];
  double np       = 100.0;                  //number of convolutions steps
  double sc       = 5;                      //convolution extends to +- gauassian sigma
  double invsq2pi = 0.3989422804014;        //(2 pi)^(-1/2)
  double xmin     = 0.0;
  double xmax     = 0.0;
  double step     = 0.0; 
  double xx       = 0.0;
  double tt[1];
  
  double sigma    = p[0];
  double result   = 0.0;
  
  // Range of the convolution integral
  xmin = ti - sc*sigma;
  xmax = ti + sc*sigma;
  step = (xmax - xmin)/ np;
  
  std::vector<double> fvalues;
  
  for( int i=0; i < np; ++i) {
    xx = xmin + i * step;
    tt[0] = xx;
    double val = (*f)(tt,fp) * TMath::Gaus(xx,ti,sigma);
    fvalues.push_back(val);
  }
  
  result = simpson( fvalues, xmin, xmax );
  result = result * (invsq2pi / sigma);

  return result;
  
}

double with2Gaussians( double (*f) (double *, double *), 
		       double *fp, double *x, double *p )
{
  
  double ti       = x[0];
  double np       = 100.0;                  //number of convolutions steps
  double sc       = 5;                      //convolution extends to +- gauassian sigma
  double invsq2pi = 0.3989422804014;        //(2 pi)^(-1/2)
  double xmin     = 0.0;
  double xmax     = 0.0;
  double step     = 0.0; 
  double xx       = 0.0;
  double tt[1];
  
  double sigma    = p[0];
  double f1       = p[1];
  double mu1      = p[2];
  double s1       = p[3];
  double mu2      = p[4];
  double s2       = p[5];
  double norm     = 1.0 / (f1*sigma*(s1-s2) + s2*sigma);
  double result1   = 0.0;
  double result2   = 0.0;
  
  // Range of the convolution integral
  xmin = (ti+mu1*sigma) - sc*s1*sigma;
  xmax = (ti+mu1*sigma) + sc*s1*sigma;
  step = (xmax - xmin)/ np;
  
  std::vector<double> fvalues;
  
  for( int i=0; i < np; ++i) {
    xx = xmin + i * step;
    tt[0] = xx;
    double val = (*f)(tt,fp) * TMath::Gaus(xx-ti,mu1*sigma, s1*sigma);
    fvalues.push_back(val);
  }
  
  result1 = simpson( fvalues, xmin, xmax );
  
  fvalues.clear();
  xmin = (ti+mu2*sigma) - sc*s2*sigma;
  xmax = (ti+mu2*sigma) + sc*s2*sigma;
  step = (xmax - xmin)/ np;
  
  for( int i=0; i < np; ++i) {
    xx = xmin + i * step;
    tt[0] = xx;
    double val = (*f)(tt,fp) *  TMath::Gaus(xx-ti,mu2*sigma, s2*sigma);
    fvalues.push_back(val);
  }
  
  result2 = simpson( fvalues, xmin, xmax );
  
  return (f1* result1 + (1-f1)*result2)*( invsq2pi * norm );
  
}



double smearValue( double (*ptr2meth) (double, double *),
		   double x , double *par)
{
  
  double result(0.);
  result = (*ptr2meth)( x , par);
  return result;
  
}

double withGaussian( double x, double *p )
{
  
  //smear value - fast implementation to account detector effects
  double mean = x;
  double sigma = p[0];
  double r1(0.);
  int seed = gRandom->GetSeed();
  gRandom->SetSeed(seed);
  r1 = gRandom->Gaus(mean,sigma);
  return r1;
  
}
