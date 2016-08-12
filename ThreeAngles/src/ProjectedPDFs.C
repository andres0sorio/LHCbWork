// $Id: ProjectedPDFs.C,v 1.1 2006/12/09 18:19:18 aosorio Exp $
// Include files 

// local
#include "PDFs.h"

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

double properTimeWpPDF(double *x, double *par)
{
  double time      = x[0];
  
  if( time < 0 ) return 0.0;
  
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

  //Bs
  double A0A0   = A1(time, K0K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2(time, KpKp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3(time, KtKt,   tauL, tauH, tauBar, dms, wphase);
  
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

  //Bar (Bs)
  double B0B0   =  B1(time, K0K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   =  B2(time, KpKp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   =  B3(time, KtKt,   tauL, tauH, tauBar, dms, wphase);
  
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

  //Bs
  double A0A0   = A1def(K0K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2def(KpKp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3def(KtKt,   tauL, tauH, tauBar, dms, wphase);
  
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

  //Bs
  double B0B0   = B1def(K0K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   = B2def(KpKp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   = B3def(KtKt,   tauL, tauH, tauBar, dms, wphase);
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

  //Time depency integrated over
  //Bs
  double A0A0   = A1def(K0K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2def(KpKp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3def(KtKt,   tauL, tauH, tauBar, dms, wphase);
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

  //Bar (Bs)
  double B0B0   = B1def(K0K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   = B2def(KpKp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   = B3def(KtKt,   tauL, tauH, tauBar, dms, wphase);
  
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
  double KpKt  = sqrt(KpKp)*sqrt(KtKt);
    
  //Time depency integrated over
  //Bs
  double A0A0   = A1def(K0K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2def(KpKp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3def(KtKt,   tauL, tauH, tauBar, dms, wphase);
  double ImApAt = A5def(KpKt, tauL, tauH, tauBar, dms, wphase);
  
  //W+
  double twophi = 2.0*phi;
  double sin2phi = TMath::Sin(twophi);
  double cos2phi = TMath::Cos(twophi);
  
  double v1 = sixteenonine * (( AtAt - ApAp ) * cos2phi 
                              + 2.0*(A0A0 + ApAp + AtAt 
                                     + ( ImApAt * sin2phi )));
  double norm = normfactorWp( par );
  
  return (1.0/norm) * v1;
  
}

double phiWmPDF(double *x, double *par)
{
  
  double phi      = x[0];
  double sixteenonine = (16.0/9.0);
  
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
  
  double KpKt  = sqrt(KpKp)*sqrt(KtKt);
  
  //Bar (Bs)
  double B0B0   = B1def(K0K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   = B2def(KpKp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   = B3def(KtKt,   tauL, tauH, tauBar, dms, wphase);
  double ImBpBt = B5def(KpKt, tauL, tauH, tauBar, dms, wphase);
  
  double twophi = 2.0*phi;
  double sin2phi = TMath::Sin(twophi);
  double cos2phi = TMath::Cos(twophi);
  
  //W-
  double v2 =  sixteenonine * (( BtBt - BpBp ) * cos2phi
                               + 2.0*(B0B0 + BpBp + BtBt 
                                      + ( ImBpBt * sin2phi )));
  
  double norm = normfactorWm( par );
  
  return (1.0/norm) * v2;
  
}
