// $Id: PDFsWRes.C,v 1.10 2007/03/22 13:05:33 aosorio Exp $
// Include files 


// local
#include "PDFsWRes.h"

/* -----------------------------------------------------------------------------
   Implementation file for class : PDFsWRes
   Convolution of PDFs
   2007-02-08 : Andres Osorio
   ----------------------------------------------------------------------------- */

// Bs->Jpsi+phi _pdf (x) + Background resolution model
double convolveTotal(  double (*ptr2model) (double (*ptr2ff) (double *, double *),
                                            double *, double *, double *),
                       double *fpar, double *x, double *cvpar ) 
{
  
  double f = convolve( totalpdf, ptr2model , fpar, x, cvpar);
  
  return f;
  
}

// Bs->Jpsi+phi _pdf (x) resolution model
double convolveSum(  double (*ptr2model) (double (*ptr2ff) (double *, double *),
                                          double *, double *, double *),
                     double *fpar, double *x, double *cvpar ) 
{
  
  double f = convolve( sumofsignals, ptr2model , fpar, x, cvpar);
  
  return f;
  
}

double convolveWp(  double (*ptr2model) (double (*ptr2ff) (double *, double *),
                                         double *, double *, double *),
                    double *fpar, double *x, double *cvpar ) 
{
  
  double f = convolve( jpsiphiWpPDF, ptr2model , fpar, x, cvpar);

  return f;
  
}

// _B_s->Jpsi+phi _pdf (x) resolution model
double convolveWm(  double (*ptr2model) (double (*ptr2ff) (double *, double *),
                                         double *, double *, double *),
                    double *fpar, double *x, double *cvpar ) 
{
  
  double f = convolve( jpsiphiWmPDF, ptr2model, fpar, x, cvpar);
  return f;
  
}

//......................................................

double totalpdfWRes(double *x, double *par)
{
  
  double convpars[6];
  convpars[0] = par[9];   //time resolution
  convpars[1] = par[10];  //f1
  convpars[2] = par[11];  //mu1
  convpars[3] = par[12];  //s1
  convpars[4] = par[13];  //mu2
  convpars[5] = par[14];  //s2
  
  double f = convolve( totalpdf , withGaussian, par, x, convpars);
  
  return f;
  
}

double jpsiphipdfWRes(double *x, double *par)
{
  
  double convpars[6];
  convpars[0] = par[9];   //time resolution
  convpars[1] = par[10];  //f1
  convpars[2] = par[11];  //mu1
  convpars[3] = par[12];  //s1
  convpars[4] = par[13];  //mu2
  convpars[5] = par[14];  //s2
  
  double f = convolve( sumofsignals , withGaussian, par, x, convpars);
  
  return f;
  
}

double properTimeWRes(double *x, double *par)
{
  
  double convpars[6];
  convpars[0] = par[9];   //time resolution
  convpars[1] = par[10];  //f1
  convpars[2] = par[11];  //mu1
  convpars[3] = par[12];  //s1
  convpars[4] = par[13];  //mu2
  convpars[5] = par[14];  //s2

  double f = convolve( properTimePDF , withGaussian, par, x , convpars);
  
  return f;
  
}

double properTimeWRes2(double *x, double *par)
{
  
  double convpars[6];
  convpars[0] = par[9];   //time resolution
  convpars[1] = par[10];  //f1
  convpars[2] = par[11];  //mu1
  convpars[3] = par[12];  //s1
  convpars[4] = par[13];  //mu2
  convpars[5] = par[14];  //s2

  double xx[8];
  xx[0]=x[0];
  xx[1]=0.0;
  xx[2]=0.0;
  xx[3]=0.0;
  xx[4]=0.0;
  xx[5]=0.0;
  xx[6]=0.0;
  xx[7]=convpars[0];

  double f = convolve( properTimePDF2 , withGaussian, par, xx , convpars);
  
  return f;
  
}
