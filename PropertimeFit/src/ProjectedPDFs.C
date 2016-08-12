// $Id: ProjectedPDFs.C,v 1.6 2007/04/09 11:49:08 aosorio Exp $
// Include files 

// local

#include "ProjectedPDFs.h"

//-----------------------------------------------------------------------------
// Implementation file for class : ProjectedPDFs
//
// 2006-12-09 : Andres Osorio
//-----------------------------------------------------------------------------

/*
  
PDF describing the single observables
Mathematica assisted the analytical integration
Andres Osorio - 14 Aug 2006 / 19 Nov 2006

*/
double properTimeLog( double *x, double *par)
{
  
  double val = properTimePDF( x, par);
  
  return par[18]*val;
    
}


double properTimePDF( double *x, double *par) 
{

  int    q        = (int)par[15];
  double epsilon[3];
  double omega    = par[8];
  
  epsilon[0] = omega;
  epsilon[1] = 0.5;
  epsilon[2] = (1.0 - omega);
  double w1  = epsilon[q+1];
  double w2  = 1.0-w1;

  double norm = w1*normfactorWpWL( par ) + w2*normfactorWmWL( par );
  
  return ( 1.0/norm )*( w1*properTimeWpPDF(x,par) 
                        + w2*properTimeWmPDF(x,par) );
  
  
}

double thetaPDF( double *x, double *par) 
{

  int    q        = (int)par[15];
  double epsilon[3];
  double omega    = par[8];
  
  epsilon[0] = omega;
  epsilon[1] = 0.5;
  epsilon[2] = (1.0 - omega);
  double w1  = epsilon[q+1];
  double w2  = 1.0-w1;
  
  double norm = w1*normfactorWpWL( par ) + w2*normfactorWmWL( par );
  
  return ( 1.0/norm )*( w1*thetaWpPDF(x,par) 
                        + w2*thetaWmPDF(x,par) );
  
}

double psiPDF( double *x, double *par) 
{

  int    q        = (int)par[15];
  double epsilon[3];
  double omega    = par[8];
  
  epsilon[0] = omega;
  epsilon[1] = 0.5;
  epsilon[2] = (1.0 - omega);
  double w1  = epsilon[q+1];
  double w2  = 1.0-w1;
  
  double norm = w1*normfactorWpWL( par ) + w2*normfactorWmWL( par );
  
  return ( 1.0/norm )*( w1*psiWpPDF(x,par) 
                        + w2*psiWmPDF(x,par) );
  
}

double phiPDF( double *x, double *par) 
{

  int    q        = (int)par[15];
  double epsilon[3];
  double omega    = par[8];
  
  epsilon[0] = omega;
  epsilon[1] = 0.5;
  epsilon[2] = (1.0 - omega);
  double w1  = epsilon[q+1];
  double w2  = 1.0-w1;
  
  double norm = w1*normfactorWpWL( par ) + w2*normfactorWmWL( par );
  
  return ( 1.0/norm )*( w1*phiWpPDF(x,par) 
                        + w2*phiWmPDF(x,par) );
  
}

double CosthetaTrPDF( double *x, double *par) 
{
  
  int    q        = (int)par[15];
  double epsilon[3];
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double ft       = par[2];
  double fp       = par[3];
  double dp1      = par[4];
  double dp2      = par[5];
  double wphase   = par[6];
  double dms      = par[7];
  double omega    = par[8];
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double KtKt  = Kt(ft,fp);
  double KpKp  = Kp(ft,fp);
  double K0K0  = K0(ft,fp);
  
  double tlow     = par[16];
  double thigh    = par[17];
 
  epsilon[0] = omega;
  epsilon[1] = 0.5;
  epsilon[2] = (1.0 - omega);
  double w1  = epsilon[q+1];
  double w2  = 1.0-w1;
  
  double costr  = x[0];
  
  double A0A0   = A1def(K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  double ApAp   = A2def(KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  double AtAt   = A3def(KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  
  double B0B0   = B1def(K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  double BpBp   = B2def(KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  double BtBt   = B3def(KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  
  double norm = 0.5*( w1*normfactorWpWL( par ) + w2*normfactorWmWL( par ) );
  
  double Wp = (2.0/3.0)*TMath::Pi()*( 3.0 *( A0A0 + ApAp ) 
                                      + (2.0*(costr*costr)-1.0)*( A0A0 + ApAp ) 
                                      + 2.0 * AtAt 
                                      - 2.0 * (2.0*(costr*costr)-1.0)* AtAt );
  
  double Wm = (2.0/3.0)*TMath::Pi()*( 3.0 *( B0B0 + BpBp ) 
                                      + (2.0 *(costr*costr)-1.0)*( B0B0 + BpBp ) 
                                      + 2.0 * BtBt 
                                      - 2.0 *(2.0*(costr*costr)-1.0)* BtBt );
  
  return ( 1.0/norm )*( w1*Wp + w2*Wm );
  
}

double CospsiTrPDF( double *x, double *par) 
{
  
  int    q        = (int)par[15];
  double epsilon[3];
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double ft       = par[2];
  double fp       = par[3];
  double dp1      = par[4];
  double dp2      = par[5];
  double wphase   = par[6];
  double dms      = par[7];
  double omega    = par[8];
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double KtKt  = Kt(ft,fp);
  double KpKp  = Kp(ft,fp);
  double K0K0  = K0(ft,fp);
  
  double tlow     = par[16];
  double thigh    = par[17];
 
  epsilon[0] = omega;
  epsilon[1] = 0.5;
  epsilon[2] = (1.0 - omega);
  double w1  = epsilon[q+1];
  double w2  = 1.0-w1;
  
  double costr  = x[0];
  
  double A0A0   = A1def(K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  double ApAp   = A2def(KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  double AtAt   = A3def(KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  
  double B0B0   = B1def(K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  double BpBp   = B2def(KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  double BtBt   = B3def(KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  
  double norm = 0.5*( w1*normfactorWpWL( par ) + w2*normfactorWmWL( par ) );
  
  double Wp = (4.0/3.0)*TMath::Pi()*( 2.0*A0A0 + ApAp 
                                      + 2.0*A0A0*(2.0*(costr*costr)-1.0) - ApAp*(2.0*(costr*costr)-1.0)
                                      + AtAt
                                      - AtAt*(2.0*(costr*costr)-1.0)
                                      );
  
  
  double Wm = (4.0/3.0)*TMath::Pi()*( 2.0*B0B0 + BpBp 
                                      + 2.0*B0B0*(2.0*(costr*costr)-1.0) - BpBp*(2.0*(costr*costr)-1.0)
                                      + BtBt
                                      - BtBt*(2.0*(costr*costr)-1.0) 
                                      );
  
  return ( 1.0/norm )*( w1*Wp + w2*Wm );
  
}

double properTimeWpPDF(double *x, double *par)
{
  
  double time      = x[0];
  
  if( time < 0 ) return 0.0;
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double ft       = par[2];
  double fp       = par[3];
  double dp1      = par[4];
  double dp2      = par[5];
  double wphase   = par[6];
  double dms      = par[7];
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);

  double KtKt  = Kt(ft,fp);
  double KpKp  = Kp(ft,fp);
  double K0K0  = K0(ft,fp);

  //Bs
  double A0A0   = A1(time, K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ApAp   = A2(time, KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double AtAt   = A3(time, KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  
  double factor =  22.3402144255274173; //Integral of the pdf w.r.t. three angles
  
  /////////////////////////////////////////
  //W+
  double v1 = A0A0 + ApAp + AtAt;
  
  return factor * v1;
  
}

double properTimeWmPDF(double *x, double *par)
{
  double time      = x[0];
  
  if( time < 0 ) return 0.0;
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double ft       = par[2];
  double fp       = par[3];
  double dp1      = par[4];
  double dp2      = par[5];
  double wphase   = par[6];
  double dms      = par[7];
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double KtKt  = Kt(ft,fp);
  double KpKp  = Kp(ft,fp);
  double K0K0  = K0(ft,fp);

  //Bar (Bs)
  double B0B0   =  B1(time, K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BpBp   =  B2(time, KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BtBt   =  B3(time, KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  
  double factor =  22.3402144255274173; //Integral of the pdf w.r.t. three angles
  
  /////////////////////////////////////////
  //W-
  double v2 = B0B0 + BpBp + BtBt;
  
  return factor * v2;
  
}


double thetaWpPDF(double *x, double *par)
{
  
  double theta    = x[0];
  double fourothreepi = (4.0/3.0)*TMath::Pi();

  double Gamma    = par[0];
  double DGrate   = par[1];
  double ft       = par[2];
  double fp       = par[3];
  double dp1      = par[4];
  double dp2      = par[5];
  double wphase   = par[6];
  double dms      = par[7];
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double KtKt  = Kt(ft,fp);
  double KpKp  = Kp(ft,fp);
  double K0K0  = K0(ft,fp);

  //Bs
  double A0A0   = A1def(K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ApAp   = A2def(KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double AtAt   = A3def(KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  
  /////////////////////////////////////////
  //W+
  double v1 = fourothreepi * sin(theta) * ( (ApAp + AtAt)*(3.0+cos(2.0*theta))
                                            + 4.0*A0A0*sin(theta)*sin(theta) );
  
  return v1;
  
}

double thetaWmPDF(double *x, double *par)
{
  
  double theta    = x[0];
  double fourothreepi = (4.0/3.0)*TMath::Pi();
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double ft       = par[2];
  double fp       = par[3];
  double dp1      = par[4];
  double dp2      = par[5];
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

  //Bs
  double B0B0   = B1def(K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BpBp   = B2def(KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BtBt   = B3def(KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  /////////////////////////////////////////
  //W-
  double v2 = fourothreepi * sin(theta) * ( (BpBp + BtBt)*(3.0+cos(2.0*theta))
                                            + 4.0 * B0B0*sin(theta)*sin(theta) );
  
  return v2;
  
}


double psiWpPDF(double *x, double *par)
{
  
  double psi      = x[0];
  double sixteenothreepi = (16.0/3.0)*TMath::Pi(); // (16/3)*pi
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double ft       = par[2];
  double fp       = par[3];
  double dp1      = par[4];
  double dp2      = par[5];
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

  //Time depency integrated over
  //Bs
  double A0A0   = A1def(K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ApAp   = A2def(KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double AtAt   = A3def(KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  //W+
  double v1 = sixteenothreepi * sin(psi)* ( 2.0*A0A0*cos(psi)*cos(psi)
                                            + (ApAp + AtAt)*sin(psi)*sin(psi));

  return v1;
  
}

double psiWmPDF(double *x, double *par)
{
  
  double psi      = x[0];
  double sixteenothreepi = (16.0/3.0)*TMath::Pi(); // (16/3)*pi
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double ft       = par[2];
  double fp       = par[3];
  double dp1      = par[4];
  double dp2      = par[5];
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

  //Bar (Bs)
  double B0B0   = B1def(K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BpBp   = B2def(KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BtBt   = B3def(KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  
  //W-
  double v2 = sixteenothreepi * sin(psi)* ( 2.0*B0B0*cos(psi)*cos(psi)
                                            + (BpBp + BtBt)*sin(psi)*sin(psi));
  
  return v2;
  
}

double phiWpPDF(double *x, double *par)
{
  
  double phi      = x[0];
  double factor2       = 0.02083333333; // (1/48)
  double factor3       = 170.666666666; // 512/3.0;
  double factor4       =  85.333333333; // 256/3.0;
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double ft       = par[2];
  double fp       = par[3];
  double dp1      = par[4];
  double dp2      = par[5];
  double wphase   = par[6];
  double dms      = par[7];
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double KtKt  = Kt(ft,fp);
  double KpKp  = Kp(ft,fp);
  double K0K0  = K0(ft,fp);
  double KpKt  = sqrt(KpKp)*sqrt(KtKt);
    
  //Time depency integrated over
  //Bs
  double A0A0   = A1def(K0K0, tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ApAp   = A2def(KpKp, tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double AtAt   = A3def(KtKt, tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ImApAt = A5def(KpKt, tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  
  //W+
  double twophi  = 2.0*phi;
  double sin2phi = TMath::Sin(twophi);
  double cos2phi = TMath::Cos(twophi);
  
  double v1 = factor2 * ( factor3 * ( A0A0 + ApAp + AtAt )
                          - factor4 * cos2phi * ( ApAp - AtAt )
                          + factor3 * sin2phi * ImApAt );
  
  return v1;
  
}

double phiWmPDF(double *x, double *par)
{
  
  double phi      = x[0];
  double factor2       = 0.02083333333; // (1/48)
  double factor3       = 170.666666666; // 512/3.0;
  double factor4       =  85.333333333; // 256/3.0;
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double ft       = par[2];
  double fp       = par[3];
  double dp1      = par[4];
  double dp2      = par[5];
  double wphase   = par[6];
  double dms      = par[7];
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double KtKt  = Kt(ft,fp);
  double KpKp  = Kp(ft,fp);
  double K0K0  = K0(ft,fp);
  double KpKt  = sqrt(KpKp)*sqrt(KtKt);
  
  //Bar (Bs)
  double B0B0   = B1def(K0K0, tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BpBp   = B2def(KpKp, tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BtBt   = B3def(KtKt, tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ImBpBt = B5def(KpKt, tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  
  double twophi  = 2.0*phi;
  double sin2phi = TMath::Sin(twophi);
  double cos2phi = TMath::Cos(twophi);
  
  //W-
  double v2 =  factor2 * ( factor3 * ( B0B0 + BpBp + BtBt )
                           - factor4 * cos2phi * ( BpBp - BtBt )
                           + factor3 * sin2phi * ImBpBt );
  return v2;
  
}


double normfactorWpWL(double *par) 
{
  
  double factor   = (64.0/9.0)*TMath::Pi();
  double Gamma    = par[0];
  double DGrate   = par[1];
  double ft       = par[2];
  double fp       = par[3];
  double dp1      = par[4];
  double dp2      = par[5];
  double wphase   = par[6];
  double dms      = par[7];

  double tlow     = par[16];
  double thigh    = par[17];
    
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double KtKt   = Kt(ft,fp);
  double KpKp   = Kp(ft,fp);
  double K0K0   = K0(ft,fp);
  
  double A0A0   = A1def(K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  double ApAp   = A2def(KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  double AtAt   = A3def(KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  
  return factor * (A0A0 + ApAp + AtAt);
  
}

double normfactorWmWL(double *par) 
{
  
  double factor   = (64.0/9.0)*TMath::Pi();
  double Gamma    = par[0];
  double DGrate   = par[1];
  double ft       = par[2];
  double fp       = par[3];
  double dp1      = par[4];
  double dp2      = par[5];
  double wphase   = par[6];
  double dms      = par[7];
    
  double tlow     = par[16];
  double thigh    = par[17];
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);

  double KtKt   = Kt(ft,fp);
  double KpKp   = Kp(ft,fp);
  double K0K0   = K0(ft,fp);

  double B0B0   = B1def(K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  double BpBp   = B2def(KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  double BtBt   = B3def(KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, tlow, thigh);
  
  return factor * (B0B0 + BpBp + BtBt);
  
}

