// $Id: PDFs.C,v 1.8 2006/12/04 16:03:41 aosorio Exp $
// Include files
// local
#include "PDFs.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PDFs
//
// 2006-06-12 : Andres Osorio Oliveros
//-----------------------------------------------------------------------------


double angpdfwres(double *x, double *par) {
  
  double convpars[6];
  convpars[0] = par[9];
  convpars[1] = par[10];   //f1
  convpars[2] = par[11];   //mu1
  convpars[3] = par[12];   //s1
  convpars[4] = par[13];   //mu2
  convpars[5] = par[14];   //s2
  
  double f = convolve( jpsiphipdf, with2Gaussians, par, x, convpars);
  
  return f;
  
} 

double properTimeConv(double *x, double *par) {
  
  double convpars[6];
  convpars[0] = par[11];
  convpars[1] = par[12];   //f1
  convpars[2] = par[13];   //mu1
  convpars[3] = par[14];   //s1
  convpars[4] = par[15];   //mu2
  convpars[5] = par[16];   //s2
  
  double f = convolve( properTimeWpPDF, with2Gaussians, par, x, convpars);
   
  return f;
  
}

double jpsiphipdf(double *x, double *par)
{
  
  //Diff. decay rate - definition as a function of t and three angles* + parameters
  //from CERN-TH 2000-101.... pg.42.
  //* the definition is in terms of dcos(theta') dcos(theta") phi
  
  double sqrtwo = TMath::Sqrt(2.0); // Sqrt[2]
  
  double time      = x[0];
  double theta     = x[1];   // l+ polar angle in Jpsi rest frame
  double psi       = x[2];   // K+ polar angle in Phi rest frame
  double phi       = x[3];   // sum of azimuthal angles
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2]; // transverse component (CP Even + CP Odd)
  double R0       = par[3]; // total longitudinal component (CP even)
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  
  double dms      = xs;
  
  double G_H=Gamma*(1.0-0.5*DGrate);
  double G_L=Gamma*(1.0+0.5*DGrate);

  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  

  double Kt    = Rt;
  double K0    = R0;
  double Kp    = 1.0 - Kt - K0;

  double K0Kp=sqrt(K0)*sqrt(Kp);
  double K0Kt=sqrt(K0)*sqrt(Kt);
  double KpKt=sqrt(Kp)*sqrt(Kt);
  
  double sqcostheta = cos(theta)*cos(theta);
  double sqsintheta = sin(theta)*sin(theta);
  double sqcospsi   = cos(psi)*cos(psi);
  double sqsinpsi   = sin(psi)*sin(psi);
  double cos2phi    = cos(2.0*phi);
  double sin2phi    = sin(2.0*phi);
  
  //Bs
  double A0A0   = A1(time, K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2(time, Kp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3(time, Kt,   tauL, tauH, tauBar, dms, wphase);
  double ReA0Ap = A4(time, K0Kp, tauL, tauH, tauBar, dms, wphase);
  double ImApAt = A5(time, KpKt, tauL, tauH, tauBar, dms, wphase);
  double ImA0At = A6(time, K0Kt, tauL, tauH, tauBar, dms, wphase);
  
  //Bar (Bs)
  double B0B0   =  B1(time, K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   =  B2(time, Kp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   =  B3(time, Kt,   tauL, tauH, tauBar, dms, wphase);
  double ReB0Bp =  B4(time, K0Kp, tauL, tauH, tauBar, dms, wphase);
  double ImBpBt =  B5(time, KpKt, tauL, tauH, tauBar, dms, wphase);
  double ImB0Bt =  B6(time, K0Kt, tauL, tauH, tauBar, dms, wphase);

  /////////////////////////////////////////
  //W+
  double v1 = ( 4.0*A0A0*sqsintheta*sqcospsi
                + ApAp* ( (1+sqcostheta)*sqsinpsi - sqsintheta*sqsinpsi*cos2phi)
                + AtAt* ( (1+sqcostheta)*sqsinpsi + sqsintheta*sqsinpsi*cos2phi)
                + 2.0 * ReA0Ap  * sqsintheta*sqsinpsi*sin2phi
                - sqrtwo*ImApAt * sin(2.0*theta) * sin(2.0*psi) * cos(phi) 
                + sqrtwo*ImA0At * sin(2.0*theta) * sin(2.0*psi) * sin(phi) ) * sin(theta) * sin(psi);

  
  //W-
  double v2 = ( 4.0*B0B0*sqsintheta*sqcospsi
                + BpBp* ( (1+sqcostheta)*sqsinpsi - sqsintheta*sqsinpsi*cos2phi)
                + BtBt* ( (1+sqcostheta)*sqsinpsi + sqsintheta*sqsinpsi*cos2phi)
                + 2.0*  ReB0Bp  * sqsintheta*sqsinpsi*sin2phi 
                - sqrtwo*ImBpBt * sin(2.0*theta) * sin(2.0*psi) * cos(phi)
                + sqrtwo*ImB0Bt * sin(2.0*theta) * sin(2.0*psi) * sin(phi) ) * sin(theta) * sin(phi);
  
  double norm = nfactorjpsiphi( par );
  
  return (1.0 / norm ) * ( (1.0-omega)*v1 + (omega*v2) );
  
  
}

double nfactorjpsiphi(double *par) 
{
  
  double factor =  (64.0/9.0)*TMath::Pi();
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double dms      = xs;
  
  double G_H=Gamma*(1.0-0.5*DGrate);
  double G_L=Gamma*(1.0+0.5*DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double Kt    = Rt;
  double K0    = R0;
  double Kp    = 1.- Kt -K0;
  
  double A0A0   = A1def(K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2def(Kp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3def(Kt,   tauL, tauH, tauBar, dms, wphase);
  //
  double B0B0   = B1def(K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   = B2def(Kp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   = B3def(Kt,   tauL, tauH, tauBar, dms, wphase);
  
  double v1 = A0A0 + ApAp + AtAt;
  double v2 = B0B0 + BpBp + BtBt;
  

  return factor * ( (1.0-omega)*v1 + (omega*v2) );
  
}

double normfactorWp(double *par) 
{
  
  double factor =  (64.0/9.0)*TMath::Pi();
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double dms      = xs;
  
  double G_H=Gamma*(1.0-0.5*DGrate);
  double G_L=Gamma*(1.0+0.5*DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
  double A0A0   = A1def(K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2def(Kp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3def(Kt,   tauL, tauH, tauBar, dms, wphase);
  
  return factor * (A0A0 + ApAp + AtAt);
  
}

double normfactorWm(double *par) 
{
  
  double factor =  (64.0/9.0)*TMath::Pi();
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double dms      = xs;
  
  double G_H=Gamma*(1.0-0.5*DGrate);
  double G_L=Gamma*(1.0+0.5*DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
  double B0B0   = B1def(K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   = B2def(Kp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   = B3def(Kt,   tauL, tauH, tauBar, dms, wphase);
  
  return factor * (B0B0 + BpBp + BtBt);
  
}

/*
  
PDF describing the single observables
Mathematica assisted the analytical integration
Andres Osorio - 14 Aug 2006 / 19 Nov 2006

*/

double properTimeWpPDF(double *x, double *par)
{
  double time      = x[0];
  
  if( time < 0 ) return 0.0;
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double dms      = xs;
  
  double G_H=Gamma*(1.0-0.5*DGrate);
  double G_L=Gamma*(1.0+0.5*DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);

  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
  //Bs
  double A0A0   = A1(time, K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2(time, Kp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3(time, Kt,   tauL, tauH, tauBar, dms, wphase);
  
  double factor =  22.3402144255274173; //Integral of the pdf w.r.t. three angles
  
  /////////////////////////////////////////
  //W+
  double v1 = A0A0 + ApAp + AtAt;
  double norm = normfactorWp( par );
  
  return (1.0/norm) * factor * v1;
  
}

double properTimeWmPDF(double *x, double *par)
{
  double time      = x[0];
  
  if( time < 0 ) return 0.0;
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double dms      = xs;
  
  double G_H=Gamma*(1.0-0.5*DGrate);
  double G_L=Gamma*(1.0+0.5*DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
  //Bar (Bs)
  double B0B0   =  B1(time, K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   =  B2(time, Kp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   =  B3(time, Kt,   tauL, tauH, tauBar, dms, wphase);
  
  double factor =  22.3402144255274173; //Integral of the pdf w.r.t. three angles
  
  /////////////////////////////////////////
  //W-
  double v2 = B0B0 + BpBp + BtBt;
  double norm = normfactorWm( par );
  
  return (1.0/norm) * factor * v2;
  
}

double thetaWpPDF(double *x, double *par)
{
  
  double theta    = x[0];
  double fourothreepi = (4.0/3.0)*TMath::Pi();
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double dms      = xs;
  
  double G_H=Gamma*(1.0-0.5*DGrate);
  double G_L=Gamma*(1.0+0.5*DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
  //Bs
  double A0A0   = A1def(K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2def(Kp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3def(Kt,   tauL, tauH, tauBar, dms, wphase);
  
  /////////////////////////////////////////
  //W+
  double v1 = fourothreepi * sin(theta) * ( (ApAp + AtAt)*(3.0+cos(2.0*theta))
                                            + 4.0*A0A0*sin(theta)*sin(theta) );
  
  double norm = normfactorWp( par );
  
  return (1.0/norm) * v1;
  
}

double thetaWmPDF(double *x, double *par)
{
  
  double theta    = x[0];
  double fourothreepi = (4.0/3.0)*TMath::Pi();
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double dms      = xs;
  
  double G_H=Gamma*(1.0-0.5*DGrate);
  double G_L=Gamma*(1.0+0.5*DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);

  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
  //Bs
  double B0B0   = B1def(K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   = B2def(Kp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   = B3def(Kt,   tauL, tauH, tauBar, dms, wphase);
  /////////////////////////////////////////
  //W-
  double v2 = fourothreepi * sin(theta) * ( (BpBp + BtBt)*(3.0+cos(2.0*theta))
                                            + 4.0 * B0B0*sin(theta)*sin(theta) );
  
  double norm = normfactorWm( par );
  
  return (1.0/norm) * v2;
  
}


double psiWpPDF(double *x, double *par)
{
  
  double psi      = x[0];
  double sixteenothreepi = (16.0/3.0)*TMath::Pi(); // (16/3)*pi
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double dms      = xs;
  
  double G_H=Gamma*(1.0-0.5*DGrate);
  double G_L=Gamma*(1.0+0.5*DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);

  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
  //Time depency integrated over
  //Bs
  double A0A0   = A1def(K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2def(Kp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3def(Kt,   tauL, tauH, tauBar, dms, wphase);
  //W+
  double v1 = sixteenothreepi * sin(psi)* ( 2.0*A0A0*cos(psi)*cos(psi)
                                            + (ApAp + AtAt)*sin(psi)*sin(psi));

  double norm = normfactorWp( par );

  return (1.0/norm) * v1;
  
}

double psiWmPDF(double *x, double *par)
{
  
  double psi      = x[0];
  double sixteenothreepi = (16.0/3.0)*TMath::Pi(); // (16/3)*pi
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double dms      = xs;
  
  double G_H=Gamma*(1.0-0.5*DGrate);
  double G_L=Gamma*(1.0+0.5*DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);

  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
  //Bar (Bs)
  double B0B0   = B1def(K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   = B2def(Kp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   = B3def(Kt,   tauL, tauH, tauBar, dms, wphase);
  
  //W-
  double v2 = sixteenothreepi * sin(psi)* ( 2.0*B0B0*cos(psi)*cos(psi)
                                            + (BpBp + BtBt)*sin(psi)*sin(psi));
  double norm = normfactorWm( par );
  
  return (1.0/norm) * v2;
  
}

double phiWpPDF(double *x, double *par)
{
  
  double phi      = x[0];
  double sixteenonine = (16.0/9.0);
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double dms      = xs;
  
  double G_H=Gamma*(1.0-0.5*DGrate);
  double G_L=Gamma*(1.0+0.5*DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
   
  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
  double K0Kp=sqrt(K0)*sqrt(Kp);
  
  //Time depency integrated over
  //Bs
  double A0A0   = A1def(K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2def(Kp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3def(Kt,   tauL, tauH, tauBar, dms, wphase);
  double ReA0Ap = A4def(K0Kp, tauL, tauH, tauBar, dms, wphase);
  
  //W+
  double twophi = 2.0*phi;
  double sin2phi = TMath::Sin(twophi);
  double cos2phi = TMath::Cos(twophi);
  
  double v1 = sixteenonine * ( ( AtAt - ApAp ) * cos2phi 
                               + 2.0*( A0A0 + ApAp + AtAt + (ReA0Ap * sin2phi) ));
  
  double norm = normfactorWp( par );
  
  return (1.0/norm) * v1;
  
}

double phiWmPDF(double *x, double *par)
{
  
  double phi      = x[0];
  double sixteenonine = (16.0/9.0);
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double dms      = xs;
  
  double G_H=Gamma*(1.0-0.5*DGrate);
  double G_L=Gamma*(1.0+0.5*DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
  double K0Kp=sqrt(K0)*sqrt(Kp);
  
  //Bar (Bs)
  double B0B0   = B1def(K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   = B2def(Kp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   = B3def(Kt,   tauL, tauH, tauBar, dms, wphase);
  double ReB0Bp = B4def(K0Kp, tauL, tauH, tauBar, dms, wphase);
  
  double twophi = 2.0*phi;
  double sin2phi = TMath::Sin(twophi);
  double cos2phi = TMath::Cos(twophi);
  
  //W-
  double v2 =  sixteenonine * ( ( BtBt - BpBp ) * cos2phi 
                                + 2.0*( B0B0 + BpBp + BtBt + (ReB0Bp * sin2phi) ));
  
  double norm = normfactorWm( par );
  
  return (1.0/norm) * v2;
  
}

///////////////////////////////////////////
//Differential function definitions

double DpsiWpPDF( double x, void * params)
{
  
  struct pdf_params *p 
    = (struct pdf_params *) params;

  const double pi = TMath::Pi();
  double eightothree = (8.0/3.0);
  
  double R0     = p->R0;
  double Rt     = p->Rt;
  double wphase = p->phis;
  double Gamma  = p->gBar;
  double DGrate = p->Dgamma;
  double dms    = p->Dms;
  
  double G_H=Gamma*(1.0-0.5*DGrate);
  double G_L=Gamma*(1.0+0.5*DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;

  double A0A0   = A1def(K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2def(Kp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3def(Kt,   tauL, tauH, tauBar, dms, wphase);

  double dv = eightothree * pi *cos(x)*(-2.0*A0A0
                                        + 3.0*(AtAt+ApAp)
                                        +(6.0*A0A0 - 3.0*(AtAt+ApAp))*cos(2.0*x));
  
  return dv;
  
}

double DpsiWmPDF( double x, void * params)
{
  
  struct pdf_params *p 
    = (struct pdf_params *) params;
  
  const double pi = TMath::Pi();
  double eightothree = (16.0/3.0);
  
  double R0     = p->R0;
  double Rt     = p->Rt;
  double wphase = p->phis;
  double Gamma  = p->gBar;
  double DGrate = p->Dgamma;
  double dms    = p->Dms;
  
  double G_H=Gamma*(1.0-0.5*DGrate);
  double G_L=Gamma*(1.0+0.5*DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
  double B0B0   = B1def(K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   = B2def(Kp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   = B3def(Kt,   tauL, tauH, tauBar, dms, wphase);
  
  double dv = eightothree * pi *cos(x)*(-2.0*B0B0
                                            + 3.0*(BtBt+BpBp)
                                            +(6.0*B0B0 - 3.0*(BtBt+BpBp))*cos(2.0*x));
  
  return dv;
  
}


double DphiWpPDF( double x, void * params)
{
  
  struct pdf_params *p 
    = (struct pdf_params *) params;

  double sixteenonine = (16.0/9.0);
  
  double R0     = p->R0;
  double Rt     = p->Rt;
  double wphase = p->phis;
  double Gamma  = p->gBar;
  double DGrate = p->Dgamma;
  double dms    = p->Dms;
  
  double G_H=Gamma*(1.0-0.5*DGrate);
  double G_L=Gamma*(1.0+0.5*DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double Kt    = Rt;
  double K0    = R0;
  double Kp    = 1.- K0 -Kt;
  double K0Kp  = sqrt(K0)*sqrt(Kp);
  
  double ApAp   = A2def(Kp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3def(Kt,   tauL, tauH, tauBar, dms, wphase);
  double ReA0Ap = A4def(K0Kp, tauL, tauH, tauBar, dms, wphase);
  
  double dv = sixteenonine * (4.0*ReA0Ap*cos(2.0*x)
                              -2.0*(AtAt-ApAp)*sin(2.0*x));
  return dv;
  
}

double DphiWmPDF( double x, void * params)
{
  
  struct pdf_params *p 
    = (struct pdf_params *) params;
  
  double sixteenonine = (16.0/9.0);

  double R0     = p->R0;
  double Rt     = p->Rt;
  double wphase = p->phis;
  double Gamma  = p->gBar;
  double DGrate = p->Dgamma;
  double dms    = p->Dms;
  
  double G_H=Gamma*(1.0-0.5*DGrate);
  double G_L=Gamma*(1.0+0.5*DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double Kt    = Rt;
  double K0    = R0;
  double Kp    = 1.- R0 -Rt;
  double K0Kp  = sqrt(K0)*sqrt(Kp);
  
  double BpBp   = B2def(Kp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   = B3def(Kt,   tauL, tauH, tauBar, dms, wphase);
  double ReB0Bp = B4def(K0Kp, tauL, tauH, tauBar, dms, wphase);
  
  double dv = sixteenonine * (4.0*ReB0Bp*cos(2.0*x)
                              -2.0*(BtBt-BpBp)*sin(2.0*x));
  
  return dv;
  
}

