// $Id: PDFs.C,v 1.16 2006/12/09 18:17:55 aosorio Exp $
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
  double ft       = par[2];
  double fp       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  
  double dms      = xs;
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double KtKt  = Kt(ft,fp);
  double KpKp  = Kp(ft,fp);
  double K0K0  = K0(ft,fp);
  
  double K0Kp  = sqrt(K0K0)*sqrt(KpKp);
  double K0Kt  = sqrt(K0K0)*sqrt(KtKt);
  double KpKt  = sqrt(KpKp)*sqrt(KtKt);
  
  double sqcostheta = cos(theta)*cos(theta);
  double sqsintheta = sin(theta)*sin(theta);
  double sqcospsi   = cos(psi)*cos(psi);
  double sqsinpsi   = sin(psi)*sin(psi);
  double cos2phi    = cos(2.0*phi);
  double sin2phi    = sin(2.0*phi);
  
  //Bs
  double A0A0   = A1(time, K0K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2(time, KpKp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3(time, KtKt,   tauL, tauH, tauBar, dms, wphase);
  double ReA0Ap = A4(time, K0Kp,   tauL, tauH, tauBar, dms, wphase);
  double ImApAt = A5(time, KpKt,   tauL, tauH, tauBar, dms, wphase);
  double ImA0At = A6(time, K0Kt,   tauL, tauH, tauBar, dms, wphase);
  
  //Bar (Bs)
  double B0B0   =  B1(time, K0K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   =  B2(time, KpKp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   =  B3(time, KtKt,   tauL, tauH, tauBar, dms, wphase);
  double ReB0Bp =  B4(time, K0Kp, tauL, tauH, tauBar, dms, wphase);
  double ImBpBt =  B5(time, KpKt, tauL, tauH, tauBar, dms, wphase);
  double ImB0Bt =  B6(time, K0Kt, tauL, tauH, tauBar, dms, wphase);

  /////////////////////////////////////////
  //W+
  double v1 = ( 4.0*A0A0*sqsintheta*sqcospsi
                + ApAp* ( (1+sqcostheta)*sqsinpsi - sqsintheta*sqsinpsi*cos2phi)
                + AtAt* ( (1+sqcostheta)*sqsinpsi + sqsintheta*sqsinpsi*cos2phi)
                + 2.0*ImApAt * sqsintheta*sqsinpsi*sin2phi
                - sqrtwo*ReA0Ap * sin(2.0*theta) * sin(2.0*psi) * cos(phi)
                + sqrtwo*ImA0At * sin(2.0*theta) * sin(2.0*psi) * sin(phi) ) ;
  
  //W-
  double v2 = ( 4.0*B0B0*sqsintheta*sqcospsi
                + BpBp* ( (1+sqcostheta)*sqsinpsi - sqsintheta*sqsinpsi*cos2phi)
                + BtBt* ( (1+sqcostheta)*sqsinpsi + sqsintheta*sqsinpsi*cos2phi)
                + 2.0*ImBpBt * sqsintheta*sqsinpsi*sin2phi
                - sqrtwo*ReB0Bp * sin(2.0*theta) * sin(2.0*psi) * cos(phi)
                + sqrtwo*ImB0Bt * sin(2.0*theta) * sin(2.0*psi) * sin(phi) );
  
  double norm = nfactorjpsiphi( par );
  
  return (1.0 / norm ) * ( (1.0-omega)*v1 + (omega*v2) );
  
  
}

double nfactorjpsiphi(double *par) 
{
  
  double factor =  (64.0/9.0)*TMath::Pi();
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double ft       = par[2];
  double fp       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double dms      = xs;
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double KtKt  = Kt(ft,fp);
  double KpKp  = Kp(ft,fp);
  double K0K0  = K0(ft,fp);

  double A0A0   = A1def(K0K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2def(KpKp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3def(KtKt,   tauL, tauH, tauBar, dms, wphase);
  //
  double B0B0   = B1def(K0K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   = B2def(KpKp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   = B3def(KtKt,   tauL, tauH, tauBar, dms, wphase);
  
  double v1 = A0A0 + ApAp + AtAt;
  double v2 = B0B0 + BpBp + BtBt;
  
  return factor * ( (1.0-omega)*v1 + (omega*v2) );
  
}

double jpsiphiWpPDF(double *x, double *par)
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
  double ft       = par[2];
  double fp       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double dms      = xs;
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double KtKt  = Kt(ft,fp);
  double KpKp  = Kp(ft,fp);
  double K0K0  = K0(ft,fp);

  double K0Kp  = sqrt(K0K0)*sqrt(KpKp);
  double K0Kt  = sqrt(K0K0)*sqrt(KtKt);
  double KpKt  = sqrt(KpKp)*sqrt(KtKt);
  
  double sqcostheta = cos(theta)*cos(theta);
  double sqsintheta = sin(theta)*sin(theta);
  double sqcospsi   = cos(psi)*cos(psi);
  double sqsinpsi   = sin(psi)*sin(psi);
  double cos2phi    = cos(2.0*phi);
  double sin2phi    = sin(2.0*phi);
  
  //Bs
  double A0A0   = A1(time, K0K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2(time, KpKp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3(time, KtKt,   tauL, tauH, tauBar, dms, wphase);
  double ReA0Ap = A4(time, K0Kp, tauL, tauH, tauBar, dms, wphase);
  double ImApAt = A5(time, KpKt, tauL, tauH, tauBar, dms, wphase);
  double ImA0At = A6(time, K0Kt, tauL, tauH, tauBar, dms, wphase);
  
  /////////////////////////////////////////
  //W+
  double v1 = ( 4.0*A0A0*sqsintheta*sqcospsi
                + ApAp* ( (1+sqcostheta)*sqsinpsi - sqsintheta*sqsinpsi*cos2phi)
                + AtAt* ( (1+sqcostheta)*sqsinpsi + sqsintheta*sqsinpsi*cos2phi)
                + 2.0*ImApAt * sqsintheta*sqsinpsi*sin2phi
                - sqrtwo*ReA0Ap * sin(2.0*theta) * sin(2.0*psi) * cos(phi)
                + sqrtwo*ImA0At * sin(2.0*theta) * sin(2.0*psi) * sin(phi)) *
    sin(theta)*sin(psi);
  
  double norm = normfactorWp( par );
  
  return (1.0 / norm ) * v1 ;
  
  
}


double jpsiphiWmPDF(double *x, double *par)
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
  double ft       = par[2];
  double fp       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double dms      = xs;
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double KtKt  = Kt(ft,fp);
  double KpKp  = Kp(ft,fp);
  double K0K0  = K0(ft,fp);
  
  double K0Kp  = sqrt(K0K0)*sqrt(KpKp);
  double K0Kt  = sqrt(K0K0)*sqrt(KtKt);
  double KpKt  = sqrt(KpKp)*sqrt(KtKt);
  
  double sqcostheta = cos(theta)*cos(theta);
  double sqsintheta = sin(theta)*sin(theta);
  double sqcospsi   = cos(psi)*cos(psi);
  double sqsinpsi   = sin(psi)*sin(psi);
  double cos2phi    = cos(2.0*phi);
  double sin2phi    = sin(2.0*phi);
  
  //Bar (Bs)
  double B0B0   =  B1(time, K0K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   =  B2(time, KpKp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   =  B3(time, KtKt,   tauL, tauH, tauBar, dms, wphase);
  double ReB0Bp =  B4(time, K0Kp, tauL, tauH, tauBar, dms, wphase);
  double ImBpBt =  B5(time, KpKt, tauL, tauH, tauBar, dms, wphase);
  double ImB0Bt =  B6(time, K0Kt, tauL, tauH, tauBar, dms, wphase);

  /////////////////////////////////////////
  //W-
  double v2 = ( 4.0*B0B0*sqsintheta*sqcospsi
                + BpBp* ( (1+sqcostheta)*sqsinpsi - sqsintheta*sqsinpsi*cos2phi)
                + BtBt* ( (1+sqcostheta)*sqsinpsi + sqsintheta*sqsinpsi*cos2phi)
                + 2.0*ImBpBt * sqsintheta*sqsinpsi*sin2phi
                - sqrtwo*ReB0Bp * sin(2.0*theta) * sin(2.0*psi) * cos(phi)
                + sqrtwo*ImB0Bt * sin(2.0*theta) * sin(2.0*psi) * sin(phi) )*
    sin(theta)*sin(psi);
  
  double norm = normfactorWm( par );
  
  return (1.0 / norm ) * v2;

}


double normfactorWp(double *par) 
{
  
  double factor =  (64.0/9.0)*TMath::Pi();
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double ft       = par[2];
  double fp       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double dms      = xs;
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double KtKt  = Kt(ft,fp);
  double KpKp  = Kp(ft,fp);
  double K0K0  = K0(ft,fp);
  
  double A0A0   = A1def(K0K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2def(KpKp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3def(KtKt,   tauL, tauH, tauBar, dms, wphase);
  
  return factor * (A0A0 + ApAp + AtAt);
  
}

double normfactorWm(double *par) 
{
  
  double factor =  (64.0/9.0)*TMath::Pi();
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double ft       = par[2];
  double fp       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double dms      = xs;
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);

  double KtKt  = Kt(ft,fp);
  double KpKp  = Kp(ft,fp);
  double K0K0  = K0(ft,fp);

  double B0B0   = B1def(K0K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   = B2def(KpKp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   = B3def(KtKt,   tauL, tauH, tauBar, dms, wphase);
  
  return factor * (B0B0 + BpBp + BtBt);
  
}
