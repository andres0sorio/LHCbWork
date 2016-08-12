// $Id: PDFsTest.C,v 1.6 2006/12/04 14:16:16 aosorio Exp $
// Include files
// local
#include "PDFsTest.h"
#include "Amplitudes.h"


//-----------------------------------------------------------------------------
// Implementation file for class : PDFs
//
// 2006-06-12 : Andres Osorio Oliveros
//-----------------------------------------------------------------------------
  
double jpsiphiPlus(double *x, double *par)
{
  
  //Diff. decay rate - definition as a function of t and three angles* + parameters
  //from CERN-TH 2000-101.... pg.42.
  //* the definition is in terms of dcos(theta') dcos(theta")
  
  double sqrtwo = TMath::Sqrt(2.0); // Sqrt[2]
  
  double time      = x[0];
  double theta     = x[1];   // l+ polar angle in Jpsi rest frame
  double psi       = x[2];   // K+ polar angle in Phi rest frame
  double phi       = x[3];   // sum of azimuthal angles
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  double wphase   = par[6];
  double xs       = par[7];
  
  double dms      = xs;
  Double_t G_H=Gamma*(1-DGrate/2.);
  Double_t G_L=Gamma*(1+DGrate/2.);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double K0    = R0; 
  double Kt    = Rt;
  double Kp    = 1.0 - Rt - R0;
  
  double K0Kp=sqrt(K0)*sqrt(Kp);
  double K0Kt=sqrt(K0)*sqrt(Kt);
  double KpKt=sqrt(Kp)*sqrt(Kt);
  
  //Bs
  double A0A0   = A1(time, K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2(time, Kp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3(time, Kt,   tauL, tauH, tauBar, dms, wphase);
  double ReA0Ap = A4(time, K0Kp, tauL, tauH, tauBar, dms, wphase);
  double ImApAt = A5(time, KpKt, tauL, tauH, tauBar, dms, wphase);
  double ImA0At = A6(time, K0Kt, tauL, tauH, tauBar, dms, wphase);
  
  std::cout << "A1: " << A0A0 << '\t' << ApAp << '\t' << AtAt << '\n';
  


  double sqcostheta = cos(theta)*cos(theta);
  double sqsintheta = sin(theta)*sin(theta);
  double sqcospsi   = cos(psi)*cos(psi);
  double sqsinpsi   = sin(psi)*sin(psi);
  double cos2phi    = cos(2.0*phi);
  double sin2phi    = sin(2.0*phi);
  double sin2theta  = sin(2.0*theta);
  double sin2psi    = sin(2.0*psi);
  
  /////////////////////////////////////////
  //W+
  Double_t v1 =  (4*A0A0*sqsintheta*sqcospsi 
                  + ApAp* ( (1+sqcostheta)*sqsinpsi - sqsintheta*sqsinpsi*cos2phi)
                  + AtAt* ( (1+sqcostheta)*sqsinpsi + sqsintheta*sqsinpsi*cos2phi)
                  + 2.0*ReA0Ap * sqsintheta*sqsinpsi*sin2phi
                  - sqrtwo*ImApAt * sin2theta * sin2psi * cos(phi)
                  + sqrtwo*ImA0At * sin2theta * sin2psi * sin(phi))*sin(theta)*sin(psi);
  
  return v1 ;
  
}

double jpsiphiMinus(double *x, double *par)
{
  
  //Diff. decay rate - definition as a function of t and three angles* + parameters
  //from CERN-TH 2000-101.... pg.42.
  //* the definition is in terms of dcos(theta') dcos(theta")
  
  //double factor = 0.04476232774; // 9 / (64 pi)

  double sqrtwo = TMath::Sqrt(2.0); // Sqrt[2]
  
  double time      = x[0];
  double theta     = x[1];   // l+ polar angle in Jpsi rest frame
  double psi       = x[2];   // K+ polar angle in Phi rest frame
  double phi       = x[3];   // sum of azimuthal angles
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  double wphase   = par[6];
  double xs       = par[7];
  
  double dms      = xs;
  Double_t G_H=Gamma*(1-DGrate/2.);
  Double_t G_L=Gamma*(1+DGrate/2.);
  
  double K0    = R0; 
  double Kt    = Rt;
  double Kp    = 1.0 - R0 - Rt;
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double K0Kp=sqrt(K0*Kp);
  double K0Kt=sqrt(K0*Kt);
  double KpKt=sqrt(Kp*Kt);

  //Bar (Bs)
  double B0B0   =  B1(time, K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   =  B2(time, Kp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   =  B3(time, Kt,   tauL, tauH, tauBar, dms, wphase);
  double ReB0Bp =  B4(time, K0Kp, tauL, tauH, tauBar, dms, wphase);
  double ImBpBt =  B5(time, KpKt, tauL, tauH, tauBar, dms, wphase);
  double ImB0Bt =  B6(time, K0Kt, tauL, tauH, tauBar, dms, wphase);

  double sqcostheta = cos(theta)*cos(theta);
  double sqsintheta = sin(theta)*sin(theta);
  double sqcospsi   = cos(psi)*cos(psi);
  double sqsinpsi   = sin(psi)*sin(psi);
  double cos2phi    = cos(2.0*phi);
  double sin2phi    = sin(2.0*phi);
  double sin2theta  = sin(2.0*theta);
  double sin2psi    = sin(2.0*psi);

  //W-
  Double_t v2 =  (4*B0B0*sqsintheta*sqcospsi 
                  + BpBp* ( (1+sqcostheta)*sqsinpsi - sqsintheta*sqsinpsi*cos2phi)
                  + BtBt* ( (1+sqcostheta)*sqsinpsi + sqsintheta*sqsinpsi*cos2phi)
                  + 2.0*ReB0Bp * sqsintheta*sqsinpsi*sin2phi
                  - sqrtwo*ImBpBt * sin2theta * sin2psi * cos(phi)
                  + sqrtwo*ImB0Bt * sin2theta * sin2psi * sin(phi))*sin(theta)*sin(psi);
  
  return v2;
  
}

