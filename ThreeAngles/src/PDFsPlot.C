// $Id: PDFsPlot.C,v 1.2 2006/11/27 14:26:03 aosorio Exp $
// Include files 



// local

#include "PDFsPlot.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PDFsPlot
//
// 2006-11-15 : Andres Felipe OSORIO OLIVEROS
//-----------------------------------------------------------------------------

double jpsiphipdfVar(double *x, double *par)
{

  //Diff. decay rate - definition as a function of t and three angles* + parameters
  //from CERN-TH .... pg.42.
  //* the definition of the diff. rate is in terms of dcos(theta') dcos(theta")
  //need to take this into account when doing the INTEGRATION over
  //theta' and theta"
  
  double time   = x[0];
  double theta  = x[1];   // cosine of l+ polar angle in Jpsi rest frame [-1;1]
  double psi    = x[2];   // cosine of K+ polar angle in Phi rest frame [-1;1]
  double phi    = x[3];   // sum of azimuthal angles
  
  int plotopt = (int)par[9];
  
  if( plotopt == 2) {
    time   = x[1];
    theta  = x[0];   // cosine of l+ polar angle in Jpsi rest frame [-1;1]
    psi    = x[2];   // cosine of K+ polar angle in Phi rest frame [-1;1]
    phi    = x[3];   // sum of azimuthal angles
  } else if ( plotopt == 3) {
    time      = x[1];
    theta     = x[2];   // cosine of l+ polar angle in Jpsi rest frame [-1;1]
    psi       = x[0];   // cosine of K+ polar angle in Phi rest frame [-1;1]
    phi       = x[3];    
  } else if ( plotopt == 4 ) {
    time      = x[1];
    theta     = x[2];   // cosine of l+ polar angle in Jpsi rest frame [-1;1]
    psi       = x[3];   // cosine of K+ polar angle in Phi rest frame [-1;1]
    phi       = x[0];
  } else {}

  double sqrtwo = TMath::Sqrt(2.0); // Sqrt[2]
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2]; // transverse component (CP Even + CP Odd)
  double Rp       = par[3]; // total longitudinal component (CP even)
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  
  double dms      = xs;
  
  double G_H=Gamma*(1.0-DGrate/2.);
  double G_L=Gamma*(1.0+DGrate/2.);

  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  

  double Kt    = Rt;
  double Kp    = (1.-Rt)*Rp;
  double K0    = (1.0-Rt)*(1.0-Rp);
  
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
                - 2.0*ImApAt * sqsintheta*sqsinpsi*sin2phi
                + sqrtwo*ReA0Ap * sin(2.0*theta) * sin(2.0*psi) * cos(phi)
                + sqrtwo*ImA0At * sin(2.0*theta) * sin(2.0*psi) * sin(phi) ) 
    * sin(theta) * sin(psi);
  
  //W-
  double v2 = ( 4.0*B0B0*sqsintheta*sqcospsi
                + BpBp* ( (1+sqcostheta)*sqsinpsi - sqsintheta*sqsinpsi*cos2phi)
                + BtBt* ( (1+sqcostheta)*sqsinpsi + sqsintheta*sqsinpsi*cos2phi)
                - 2.0*ImBpBt * sqsintheta*sqsinpsi*sin2phi
                + sqrtwo*ReB0Bp * sin(2.0*theta) * sin(2.0*psi) * cos(phi)
                + sqrtwo*ImB0Bt * sin(2.0*theta) * sin(2.0*psi) * sin(phi) )
    * sin(theta) * sin(psi);
  
  
  double norm = nfactor( par );
  
  
  return (1.0/norm) * ( ((1.0 - omega)*v1) + (omega)*v2);
  
}

double nfactor(double *par) 
{
  
  double factor =  (64.0/9.0)*TMath::Pi();
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double dms      = xs;
  
  double G_H=Gamma*(1.0-DGrate/2.);
  double G_L=Gamma*(1.0+DGrate/2.);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);

  double Kt    = Rt;
  double Kp    = (1.-Rt)*Rp;
  double K0    = (1.0-Rt)*(1.0-Rp);
  
  double A0A0   = A1(K0,   tauL, tauH, tauBar, dms, wphase);
  double ApAp   = A2(Kp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3(Kt,   tauL, tauH, tauBar, dms, wphase);
  //
  double B0B0   = B1(K0,   tauL, tauH, tauBar, dms, wphase);
  double BpBp   = B2(Kp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   = B3(Kt,   tauL, tauH, tauBar, dms, wphase);
  
  double v1 = A0A0 + ApAp + AtAt;
  double v2 = B0B0 + BpBp + BtBt;
  
  // return factor * ( (omega*v1) + (1.0-omega)*v2 );
  return factor * ( (1.0-omega)*v1 + (omega*v2) );
  
}
