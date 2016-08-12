#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>

#include "PDFs.h"
#include "PDFsTest.h"
#include "Amplitudes.h"
#include "PDFsWacc.h"

int main(int iargv, const char **argv) {
  
  if(iargv < 2 ) {
    std::cout << "usage : pdfValues [input file]" << std::endl;
    return 0;
  }
  
  
  const char *infile  = argv[1];
  
  double pars[15];

  pars[0]= 0.700;
  pars[1]= 0.200;
  pars[2]= 0.200;
  pars[3]= 0.600;
  pars[4]= 0.000;
  pars[5]= 3.14159265358979312;
  pars[6]= -0.04;
  pars[7]= 20;
  pars[8]= 0.30;
  pars[9]= 0.00;
  pars[10]= 1.0;
  pars[11]= 0.0;
  pars[12]= 0.0;
  pars[13]= 0.0;
  pars[14]= 0.0;
  

  double Kt    = pars[2];
  double K0    = pars[3]; 
  double Kp    = 1.0 - Kt - K0;
  
  std::cout << setprecision(10) << "Kt:  " << Kt
            << " K0: " << K0
            << " Kp: " << Kp << std::endl;
  

  double K0Kp=sqrt(K0*Kp);
  double K0Kt=sqrt(K0*Kt);
  double KpKt=sqrt(Kp*Kt);

  double Gamma = pars[0];
  double DGrate=pars[1];
  double G_H=pars[0]*(1.0-DGrate/2.);
  double G_L=pars[0]*(1.0+DGrate/2.);
  
  double tauL = (1.0/G_L);
  double tauH = (1.0/G_H);
  double tauBar = (1.0/Gamma);
  
  double a1(0.), a2(0.), a3(0.);
  double a4(0.), a5(0.), a6(0.);

  double b1(0.), b2(0.), b3(0.);
  double b4(0.), b5(0.), b6(0.);
  
  double xx[4];

  ifstream in;
  in.open(infile);
  
  if(!in.is_open()) {
    std::cout << "could not open input file." << std::endl;
    return 0;
  }
  
  std::cout << "opening files:" << std::endl;
  
  char comment[256];
  in.getline (comment,256);
  std::cout << std::string(comment) << std::endl;
  
  double value(0.0);
  double norm(0.0);
  double v1(0.0);
  double v2(0.0);

  while (1) {
    in >> xx[0] >> xx[1] >> xx[2] >> xx[3];
    if (!in.good()) break;
    value = jpsiphipdf( xx, pars);
    norm  = nfactorjpsiphi( pars);
    v1    = jpsiphiPlus( xx, pars);
    v2    = jpsiphiMinus( xx, pars);

    a1    = A1(xx[0],K0,tauL,tauH,tauBar,pars[7],pars[6]);
    a2    = A2(xx[0],Kp,tauL,tauH,tauBar,pars[7],pars[6]);
    a3    = A3(xx[0],Kt,tauL,tauH,tauBar,pars[7],pars[6]);
    a4    = A4(xx[0],K0Kp,tauL,tauH,tauBar,pars[7],pars[6]);
    a5    = A5(xx[0],KpKt,tauL,tauH,tauBar,pars[7],pars[6]);  
    a6    = A6(xx[0],K0Kt,tauL,tauH,tauBar,pars[7],pars[6]);

    b1    = B1(xx[0],K0,tauL,tauH,tauBar,pars[7],pars[6]);
    b2    = B2(xx[0],Kp,tauL,tauH,tauBar,pars[7],pars[6]);
    b3    = B3(xx[0],Kt,tauL,tauH,tauBar,pars[7],pars[6]);
    b4    = B4(xx[0],K0Kp,tauL,tauH,tauBar,pars[7],pars[6]);
    b5    = B5(xx[0],KpKt,tauL,tauH,tauBar,pars[7],pars[6]);  
    b6    = B6(xx[0],K0Kt,tauL,tauH,tauBar,pars[7],pars[6]);
    
    std::cout << "normalisation: " << norm  << '\t';
    std::cout << "function value: " << value << '\t';
    std::cout << "W+: " << v1 << '\t';
    std::cout << "W-: " << v2 << std::endl;
  
    std::cout << "a1 " << a1 << std::endl;
    std::cout << "a2 " << a2 << std::endl;
    std::cout << "a3 " << a3 << std::endl;
    std::cout << "a4 " << a4 << std::endl;
    std::cout << "a5 " << a5 << std::endl;
    std::cout << "a6 " << a6 << std::endl;
    std::cout << "b1 " << b1 << std::endl;
    std::cout << "b2 " << b2 << std::endl;
    std::cout << "b3 " << b3 << std::endl;
    std::cout << "b4 " << b4 << std::endl;
    std::cout << "b5 " << b5 << std::endl;
    std::cout << "b6 " << b6 << std::endl;
    
    a1    = A1def(K0,tauL,tauH,tauBar,pars[7],pars[6]);
    a2    = A2def(Kp,tauL,tauH,tauBar,pars[7],pars[6]);
    a3    = A3def(Kt,tauL,tauH,tauBar,pars[7],pars[6]);
    a4    = A4def(K0Kp,tauL,tauH,tauBar,pars[7],pars[6]);
    a5    = A5def(KpKt,tauL,tauH,tauBar,pars[7],pars[6]);
    
    b1    = B1def(K0,tauL,tauH,tauBar,pars[7],pars[6]);
    b2    = B2def(Kp,tauL,tauH,tauBar,pars[7],pars[6]);
    b3    = B3def(Kt,tauL,tauH,tauBar,pars[7],pars[6]);
    b4    = B4def(K0Kp,tauL,tauH,tauBar,pars[7],pars[6]);  
    b5    = B5def(KpKt,tauL,tauH,tauBar,pars[7],pars[6]);  

    std::cout << "Now their integrated values" << std::endl;
        
    std::cout << "a1 " << a1 << std::endl;
    std::cout << "a2 " << a2 << std::endl;
    std::cout << "a3 " << a3 << std::endl;
    std::cout << "a4 " << a4 << std::endl;
    std::cout << "a5 " << a5 << std::endl;
    
    std::cout << "b1 " << b1 << std::endl;
    std::cout << "b2 " << b2 << std::endl;
    std::cout << "b3 " << b3 << std::endl;
    std::cout << "b4 " << b4 << std::endl;
    std::cout << "b5 " << b5 << std::endl;
    
    std::cout << "The distributions: PT, theta, psi, phi" << std::endl;

    std::cout << "Bs:" << std::endl;
    
    std::cout << properTimeWpPDF(&xx[0],pars) << '\t'
              << thetaWpPDF(&xx[1],pars) << '\t'
              << psiWpPDF(&xx[2],pars) << '\t'
              << phiWpPDF(&xx[3],pars) << std::endl;

    std::cout << "Bs_bar:" << std::endl;
    
    std::cout << properTimeWmPDF(&xx[0],pars) << '\t'
              << thetaWmPDF(&xx[1],pars) << '\t'
              << psiWmPDF(&xx[2],pars) << '\t'
              << phiWmPDF(&xx[3],pars) << std::endl;
    
    



    //Testing numerical integration:
    
    std::cout << "pTimeAcc: " << pTimeAcc( xx , pars ) << std::endl;
    std::cout << "normalisation: " << accnormfactor(pars) << std::endl;
    std::cout << "pTimePDF: " << pTimeAccPDF(xx, pars) << std::endl;
    
    
    
    






    
  }    
  
  in.close();
  
  return 0;
  
}
