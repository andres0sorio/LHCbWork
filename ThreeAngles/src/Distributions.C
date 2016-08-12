// $Id: Distributions.C,v 1.2 2006/11/20 17:27:43 aosorio Exp $
// Include files 



// local
#include "Distributions.h"

//-----------------------------------------------------------------------------
// Implementation file for class : Distributions
//
// 2006-11-19 : Andres Felipe OSORIO OLIVEROS
//-----------------------------------------------------------------------------

double jpsiphi(double *x, double *par)
{

  //Diff. decay rate - definition as a function of t and three angles* + parameters
  //from CERN-TH 2000-101.... pg.42.
  //* the definition is in terms of dcos(theta') dcos(theta") phi
  
  double sqrtwo = 1.41421356237; // Sqrt[2]
  
  double time      = x[0];
  double theta     = x[1];   // l+ polar angle in Jpsi rest frame
  double psi       = x[2];   // K+ polar angle in Phi rest frame
  double phi       = x[3];   // sum of azimuthal angles
  
  double Gamma    = par[0];
  double DGamma   = par[1];
  double Rt       = par[2]; // transverse component (CP Even + CP Odd)
  double R0       = par[3]; // total longitudinal component (CP even)
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  
  double dms      = xs;
  
  double DGrate=DGamma/Gamma;
  double G_H=Gamma*(1.0-DGrate/2.);
  double G_L=Gamma*(1.0+DGrate/2.);
  
  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
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
  double A0A0   = A1(time, K0,   G_L, G_H, Gamma, dms, wphase);
  double ApAp   = A2(time, Kp,   G_L, G_H, Gamma, dms, wphase);
  double AtAt   = A3(time, Kt,   G_L, G_H, Gamma, dms, wphase);
  double ReA0Ap = A4(time, K0Kp, G_L, G_H, Gamma, dms, wphase);
  double ImApAt = A5(time, KpKt, G_L, G_H, Gamma, dms, wphase);
  double ImA0At = A6(time, K0Kt, G_L, G_H, Gamma, dms, wphase);
  
  //Bar (Bs)
  double B0B0   =  B1(time, K0,   G_L, G_H, Gamma, dms, wphase);
  double BpBp   =  B2(time, Kp,   G_L, G_H, Gamma, dms, wphase);
  double BtBt   =  B3(time, Kt,   G_L, G_H, Gamma, dms, wphase);
  double ReB0Bp =  B4(time, K0Kp, G_L, G_H, Gamma, dms, wphase);
  double ImBpBt =  B5(time, KpKt, G_L, G_H, Gamma, dms, wphase);
  double ImB0Bt =  B6(time, K0Kt, G_L, G_H, Gamma, dms, wphase);
  
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

  
  return ( (1.0-omega)*v1 + (omega*v2) );
  
}


double properTime(double *x, double *par)
{

  //jpsiphi pdf three angles: this is the propertime component
  //the three angles are integrated
  //The MC integral matches this function
  
  double time      = x[0];
  
  if( time < 0 ) return 0.0;
  
  double Gamma    = par[0];
  double DGamma   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double scale    = par[10];
  
  double dms      = xs;
  
  double DGrate   = DGamma/Gamma;
  double G_H=Gamma*(1.0-DGrate/2.);
  double G_L=Gamma*(1.0+DGrate/2.);
  
  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;

  scale   = 1.0;
  
  //Bs
  double A0A0   = A1(time, K0,   G_L, G_H, Gamma, dms, wphase);
  double ApAp   = A2(time, Kp,   G_L, G_H, Gamma, dms, wphase);
  double AtAt   = A3(time, Kt,   G_L, G_H, Gamma, dms, wphase);
  
  //Bar (Bs)
  double B0B0   =  B1(time, K0,   G_L, G_H, Gamma, dms, wphase);
  double BpBp   =  B2(time, Kp,   G_L, G_H, Gamma, dms, wphase);
  double BtBt   =  B3(time, Kt,   G_L, G_H, Gamma, dms, wphase);
  
  double factor =  22.3402144255274173; //Integral of the pdf w.r.t. three angles

  /////////////////////////////////////////
  //W+
  double v1 = A0A0 + ApAp + AtAt;
  //W-
  double v2 = B0B0 + BpBp + BtBt;
  
  return scale * factor * ( (1.0-omega)*v1 + (omega*v2) );
  
}

//=============================================================================
