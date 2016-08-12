// Include files
// local

#include "PDFs.h"

//-----------------------------------------------------------------------------
// 2006-11-28 : Andres Osorio Oliveros
//-----------------------------------------------------------------------------


double TotalPDF(double *x, double *par)
{
  //pm[0]=w;
  //pm[1]=q;
  //pm[2]=Amix;
  //pm[3]=DGamma;
  //pm[4]=GammaBAR;
  //pm[5]=deltaM;
  //pm[6]=mean;
  //pm[7]=width;
  //pm[8]=sigmaPT;
  //pm[9]=Rt;
  //pm[10]=f1;
  //pm[11]=s1;
  //pm[12]=s2;
  //pm[13]=mu1;
  //pm[14]=mu2;
 
  double time     = x[0];
  double tr       = x[1];
  
  double w        = par[0];
  double q        = par[1];
  double Amix     = par[2];
  double DGamma   = par[3];
  double GammaBar = par[4];
  double deltam   = par[5];
  double CPratio  = par[9];
   
  double GammaL   = GammaBar +  0.5*DGamma;
  double GammaH   = GammaBar -  0.5*DGamma;
  
  double signalPtODD  = 0.0;
  double signalPtEVEN = 0.0;
  double normterm     = 0.0;
  double part1        = 0.0;
  double part2        = 0.0;
  
  double AH = TMath::Exp(-1.0 * GammaH * time);
  double AL = TMath::Exp(-1.0 * GammaL * time);

  double A2 = (1.0 - 2.0 * w) * q * Amix;
  double A3 = TMath::Exp(-1.0 * GammaBar*time)*TMath::Sin(deltam*time);
    
  signalPtODD  = AH + (A2*A3);

  signalPtEVEN = AL - (A2*A3);
  
  //Normalisation factor only for pdf
  
  double tmin    = 0.20;
  double tmax    = 20.0;
  
  double gammaDm = pow ( (GammaBar*GammaBar + deltam*deltam) , -1.0 );
  
  double intOddA = (1/GammaH)*exp(-GammaH*tmin)
    + (1.0-2.0*w)*q*Amix*gammaDm*exp(-GammaBar*tmin)*(GammaBar*sin(deltam*tmin)
                                                      + deltam*cos(deltam*tmin));
  
  double intEvenA = (1/GammaL)*exp(-GammaL*tmin)
    - (1.0-2.0*w)*q*Amix*gammaDm*exp(-GammaBar*tmin)*(GammaBar*sin(deltam*tmin)
                                                      + deltam*cos(deltam*tmin));
  
  double intOddB = (1/GammaH)*exp(-GammaH*tmax)
    + (1.0-2.0*w)*q*Amix*gammaDm*exp(-GammaBar*tmax)*(GammaBar*sin(deltam*tmax)
                                                      + deltam*cos(deltam*tmax));
  
  double intEvenB = (1/GammaL)*exp(-GammaL*tmax)
    - (1.0-2.0*w)*q*Amix*gammaDm*exp(-GammaBar*tmax)*(GammaBar*sin(deltam*tmax)
                                                      + deltam*cos(deltam*tmax));
  
  normterm = (CPratio)*(1.33333)*(intOddA-intOddB) 
    + (1.0-CPratio)*(1.333333)*(intEvenA-intEvenB);
  
  part1 = (CPratio)*(1.0-cos(tr)*cos(tr))*signalPtODD;
  
  part2 = 0.5*(1.0 - CPratio)*(1.0+cos(tr)*cos(tr))*signalPtEVEN;
  
  return (1.0/normterm) * (part1 + part2);
  
}
