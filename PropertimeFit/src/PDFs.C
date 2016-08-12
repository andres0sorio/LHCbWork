// $Id: PDFs.C,v 1.17 2007/04/05 08:53:46 aosorio Exp $
// Include files
// local
#include "PDFs.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PDFs
//
// 2006-06-12 : Andres Osorio Oliveros
//-----------------------------------------------------------------------------

double totalpdf ( double *x, double *par )
{
  //total pdf != to the three-angle version
  
  double fsig = par[16];
  double arg1 = properTimePDF( x , par );
  double arg2 = bkgPDF( x , par );

  return fsig*arg1 + (1.0-fsig)*arg2;
  
}

double bsmassPDF( double *x, double *par )
{
  
  double mass  = x[4];
  double bmass = 5.369;
  double mres  = 0.010;
  double arg   = TMath::Gaus(mass,bmass,mres,1);
  
  return arg;
  
}

double sumofsignals( double *x, double *par )
{
  
  double norm  = x[6];

  //In transversity basis
  norm = 0.5 * norm;
  
  return (1.0/norm) * jpsiphipdf (x , par);
  
}

double jpsiphipdf(double *x, double *par)
{
  
  //Diff. decay rate - definition as a function of t and three angles* + parameters
  //from CERN-TH 2000-101.... pg.42.
  //* the definition is in terms of dcos(theta') dcos(theta") phi
 
  double time      = x[0];
  double theta     = x[1];   // l+ polar angle in Jpsi rest frame
  double psi       = x[2];   // K+ polar angle in Phi rest frame
  double phi       = x[3];   // sum of azimuthal angles
  int    q         = (int)x[5];

  if ( time < 0 ) return 0.0;
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double ft       = par[2];
  double fp       = par[3];
  double dp1      = par[4];
  double dp2      = par[5];
  double wphase   = par[6];
  double dms      = par[7];
  double omega    = par[8];

  double epsilon[3];
  epsilon[0] = omega;
  epsilon[1] = 0.5;
  epsilon[2] = (1.0 - omega);
  double w1  = epsilon[q+1];
  double w2  = 1.0-w1;
      
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
  
  //Bs
  double A0A0   = A1(time, K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ApAp   = A2(time, KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double AtAt   = A3(time, KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ReA0Ap = A4(time, K0Kp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ImApAt = A5(time, KpKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ImA0At = A6(time, K0Kt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  
  //Bar (Bs)
  double B0B0   =  B1(time, K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BpBp   =  B2(time, KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BtBt   =  B3(time, KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ReB0Bp =  B4(time, K0Kp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ImBpBt =  B5(time, KpKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ImB0Bt =  B6(time, K0Kt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  
  /////////////////////////////////////////
  //W+
  double v1 = ( A0A0*f1(theta,psi,phi) 
                + ApAp*f2(theta,psi,phi) 
                + AtAt*f3(theta,psi,phi)
                + ImApAt*f4(theta,psi,phi) 
                + ReA0Ap*f5(theta,psi,phi)
                + ImA0At*f6(theta,psi,phi) );
  
  //W-
  double v2 = ( B0B0*f1(theta,psi,phi)
                + BpBp*f2(theta,psi,phi)
                + BtBt*f3(theta,psi,phi)
                + ImBpBt*f4(theta,psi,phi)
                + ReB0Bp*f5(theta,psi,phi)
                + ImB0Bt*f6(theta,psi,phi) );
  
  return ( w1*v1 + w2*v2 );
  
}

double nfactorjpsiphi(double *par) 
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
  double omega    = par[8];
  int    q        = (int)par[15];

  double epsilon[3];
  epsilon[0] = omega;
  epsilon[1] = 0.5;
  epsilon[2] = (1.0 - omega);
  double w1  = epsilon[q+1];
  double w2  = 1.0-w1;
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double KtKt   = Kt(ft,fp);
  double KpKp   = Kp(ft,fp);
  double K0K0   = K0(ft,fp);
  
  double A0A0   = A1def(K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ApAp   = A2def(KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double AtAt   = A3def(KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  //
  double B0B0   = B1def(K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BpBp   = B2def(KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BtBt   = B3def(KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  
  double v1 = A0A0 + ApAp + AtAt;
  double v2 = B0B0 + BpBp + BtBt;

  return factor * ( w1*v1 + w2*v2 );
  
}

double jpsiphiWpPDF(double *x, double *par)
{
  
  //Diff. decay rate - definition as a function of t and three angles* + parameters
  //from CERN-TH 2000-101.... pg.42.
  //* the definition is in terms of dcos(theta') dcos(theta") phi
  
  double time      = x[0];
  double theta     = x[1];   // l+ polar angle in Jpsi rest frame
  double psi       = x[2];   // K+ polar angle in Phi rest frame
  double phi       = x[3];   // sum of azimuthal angles
  
  if ( time < 0 ) return 0.0;
    
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
  
  double KtKt   = Kt(ft,fp);
  double KpKp   = Kp(ft,fp);
  double K0K0   = K0(ft,fp);
  
  double K0Kp   = sqrt(K0K0)*sqrt(KpKp);
  double K0Kt   = sqrt(K0K0)*sqrt(KtKt);
  double KpKt   = sqrt(KpKp)*sqrt(KtKt);
  
  //Bs
  double A0A0   = A1(time, K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ApAp   = A2(time, KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double AtAt   = A3(time, KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ReA0Ap = A4(time, K0Kp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ImApAt = A5(time, KpKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ImA0At = A6(time, K0Kt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  
  /////////////////////////////////////////
  //W+
  double v1 = ( A0A0*h1(theta,psi,phi) 
                + ApAp*h2(theta,psi,phi) 
                + AtAt*h3(theta,psi,phi)
                + ImApAt*h4(theta,psi,phi) 
                + ReA0Ap*h5(theta,psi,phi)
                + ImA0At*h6(theta,psi,phi) );
  
  return v1;
  
  //double norm = normfactorWp( par );
  //return (1.0 / norm ) * v1;

}


double jpsiphiWmPDF(double *x, double *par)
{
  
  //Diff. decay rate - definition as a function of t and three angles* + parameters
  //from CERN-TH 2000-101.... pg.42.
  //* the definition is in terms of dcos(theta') dcos(theta") phi
  
  double time      = x[0];
  double theta     = x[1];   // l+ polar angle in Jpsi rest frame
  double psi       = x[2];   // K+ polar angle in Phi rest frame
  double phi       = x[3];   // sum of azimuthal angles
  
  if ( time < 0 ) return 0.0;
  
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
  
  double KtKt   = Kt(ft,fp);
  double KpKp   = Kp(ft,fp);
  double K0K0   = K0(ft,fp);
  
  double K0Kp   = sqrt(K0K0)*sqrt(KpKp);
  double K0Kt   = sqrt(K0K0)*sqrt(KtKt);
  double KpKt   = sqrt(KpKp)*sqrt(KtKt);
  
  //Bar (Bs)
  double B0B0   =  B1(time, K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BpBp   =  B2(time, KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BtBt   =  B3(time, KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ReB0Bp =  B4(time, K0Kp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ImBpBt =  B5(time, KpKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ImB0Bt =  B6(time, K0Kt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  
  /////////////////////////////////////////
  //W-
  double v2 = ( B0B0*h1(theta,psi,phi)
                + BpBp*h2(theta,psi,phi)
                + BtBt*h3(theta,psi,phi)
                + ImBpBt*h4(theta,psi,phi)
                + ReB0Bp*h5(theta,psi,phi)
                + ImB0Bt*h6(theta,psi,phi) );

  return v2;

  //double norm = normfactorWm( par );
  //return (1.0 / norm ) * v2;
  
}


double normfactorWp(double *par) 
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
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double KtKt   = Kt(ft,fp);
  double KpKp   = Kp(ft,fp);
  double K0K0   = K0(ft,fp);
  
  double A0A0   = A1def(K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ApAp   = A2def(KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double AtAt   = A3def(KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  
  return factor * (A0A0 + ApAp + AtAt);
  
}

double normfactorWm(double *par) 
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
    
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);

  double KtKt   = Kt(ft,fp);
  double KpKp   = Kp(ft,fp);
  double K0K0   = K0(ft,fp);

  double B0B0   = B1def(K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BpBp   = B2def(KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BtBt   = B3def(KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  
  return factor * (B0B0 + BpBp + BtBt);
  
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

  double norm  = x[6];
  //In transversity basis
  norm = 0.5*norm;

  double term1 = properTimeWpPDF(x,par);
  double term2 = properTimeWmPDF(x,par);
  
  double val = ( 1.0/norm )*( w1*term1 + w2*term2 );
  
  //std::cout << val << " w1 " << w1 << " w2 " << w2 <<  " n: " <<  norm << std::endl;
  //std::cout << term1 << " " << term2 << std::endl;
  //std::cout << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << std::endl;
  
  return val;
  
  
  
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
  double A0A0   = A1(time, 1.0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double ApAp   = A2(time, 1.0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double AtAt   = A3(time, 1.0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  
  double factor =  22.3402144255274173; //Integral of the pdf w.r.t. three angles
  //In transversity basis
  factor = 0.5*factor;
  
  /////////////////////////////////////////
  //W+
  double v1 = (A0A0 + ApAp)*(1.0-KtKt) + AtAt*(KtKt);
  
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
  double B0B0   =  B1(time, 1.0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BpBp   =  B2(time, 1.0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  double BtBt   =  B3(time, 1.0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2);
  
  double factor =  22.3402144255274173; //Integral of the pdf w.r.t. three angles
  //In transversity basis
  factor = 0.5*factor;
    
  /////////////////////////////////////////
  //W-
  double v2 = (B0B0 + BpBp)*(1.0-KtKt) + BtBt*KtKt;
  
  return factor * v2;
  
}

double properTimePDF2( double *x, double *par) 
{
  
  int    q        = (int)par[15];
  double epsilon[3];
  double omega    = par[8];
  
  epsilon[0] = omega;
  epsilon[1] = 0.5;
  epsilon[2] = (1.0 - omega);
  double w1  = epsilon[q+1];
  double w2  = 1.0-w1;

  double fsig = par[16];

  double norm  = x[6];
  //In transversity basis
  norm = 0.5 * nfactorjpsiphiWL( par );
  
  double term1 = properTimeWpPDF(x,par);
  double term2 = properTimeWmPDF(x,par);
  
  double signal = ( 1.0/norm )*( w1*term1 + w2*term2 );
  double bkg    = bkgPDF( x , par );
  double val    = fsig * signal + (1.0-fsig)* bkg;
    
  //std::cout << val << " w1 " << w1 << " w2 " << w2 <<  " n: " <<  norm << std::endl;
  //std::cout << term1 << " " << term2 << std::endl;
  //std::cout << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << std::endl;
  
  return val;
  
  
  
}


double nfactorjpsiphiWL(double *par) 
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
  double omega    = par[8];
  int    q        = (int)par[15];

  double epsilon[3];
  epsilon[0] = omega;
  epsilon[1] = 0.5;
  epsilon[2] = (1.0 - omega);
  double w1  = epsilon[q+1];
  double w2  = 1.0-w1;
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double KtKt   = Kt(ft,fp);
  double KpKp   = Kp(ft,fp);
  double K0K0   = K0(ft,fp);
  
  double A0A0   = A1def(K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, 0.0, par[19]);
  double ApAp   = A2def(KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, 0.0, par[19]);
  double AtAt   = A3def(KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, 0.0, par[19]);
  //
  double B0B0   = B1def(K0K0,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, 0.0, par[19]);
  double BpBp   = B2def(KpKp,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, 0.0, par[19]);
  double BtBt   = B3def(KtKt,   tauL, tauH, tauBar, dms, wphase, dp1, dp2, 0.0, par[19]);
  
  double v1 = A0A0 + ApAp + AtAt;
  double v2 = B0B0 + BpBp + BtBt;

  return factor * ( w1*v1 + w2*v2 );
  
}
