#include "Test.h"

double afunction( double *x, double *coef) {
  
  //1 var function
  
  double f = TMath::Exp(-coef[0]*x[0])+coef[1];
  return f;
    
}

double afunction2( double *x, double *coef) {
  
  //1 var function
  
  double f = coef[0] / (1.0 + x[0]*x[0]);
  return f;
  
}

double afunction3( double *x, double *coef) {
  
  //1 var function
  
  double f = coef[0] / (1.0 + x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  return f;
  
}


double a2dfunction( double *x, double *coef) {

  //2 var function
  
  double f = coef[0]*TMath::Cos(x[0])*TMath::Sin(x[1])*TMath::Sin(x[1]);
  return f;

}

double a3dfunction( double *x, double *coef) {

  //2 var function
  
  double f = coef[0]*TMath::Cos(x[0])*TMath::Sin(x[1])*TMath::Sin(x[1])
    + TMath::Cos(x[2])*TMath::Sin(x[2]);

  return f;

}

void test() {

  
  double xx[3];
  xx[0] = 0.0;
  xx[1] = 0.0;
  xx[2] = 0.0;

  double pp[2];
  pp[0] = 0.8;
  pp[1] = 1.2;

  double area = integrate(afunction, xx, pp, 0.0, 10.0 );
  
  std::cout << area << std::endl;

  double inf[3];
  double sup[3];
  inf[0] = 0.0;
  inf[1] = 0.0;
  sup[0] = 10.0;
  sup[1] = 10.0;
  inf[2] =  0.0;
  sup[2] = 10.0;

  double result1 = integrate2d (a2dfunction, xx, pp, inf, sup );
  double result2 = integrate3d (a3dfunction, xx, pp, inf, sup );
  
  std::cout << result1 << std::endl;
  std::cout << result2 << std::endl;
  
}

void testMCI ()
{
  
  double inf[1];
  double sup[1];
  inf[0] = 0.0;
  sup[0] = 1.0;

  double x[1];
  double p[1];
  x[0]=0.0;
  p[0]=1.0;
  

  TH1D *h1 = new TH1D("h1","integration results", 50, 0.7, 0.8);

  double result1=0.0;
  
  for (int i = 0; i < 1000; ++i) {
    result1 = mcintegrate(afunction2, x , p , inf, sup, 100000);
    h1->Fill(result1);
  }
  
  h1->Draw();
  
  std::cout << result1 << std::endl;

}

void testMCI3d ()
{
  
  double inf[3];
  double sup[3];
  inf[0] = 0.0;
  sup[0] = 1.0;
  inf[1] = 0.0;
  sup[1] = 1.0;
  inf[2] = 0.0;
  sup[2] = 1.0;

  double x[3];
  double p[1];
  x[0]=0.0;
  x[1]=0.0;
  x[2]=0.0;

  p[0]=1.0;
  
  TH1D *h1 = new TH1D("h1","integration results", 50, 0.4, 0.6);
  
  double result1=0.0;
  
  for (int i = 0; i < 1000; ++i) {
    result1 = mcintegrate3d(afunction3, x , p , inf, sup, 100000);
    h1->Fill(result1);
  }
  
  h1->Draw();
  
  std::cout << result1 << std::endl;

  double result2 = integrate3d(afunction3, x , p , inf, sup);

  std::cout << result2 << std::endl;

}
