// $Id: PDFs_Helicity.C,v 1.3 2006/11/24 22:12:55 aosorio Exp $
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
  
  double f = convolve( properTimePDF, with2Gaussians, par, x, convpars);
   
  return f;
  
}

double jpsiphipdf(double *x, double *par)
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
  double DGrate   = par[1];
  double Rt       = par[2]; // transverse component (CP Even + CP Odd)
  double R0       = par[3]; // total longitudinal component (CP even)
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  
  double dms      = xs;
  
  double G_H=Gamma*(1.0-DGrate/2.);
  double G_L=Gamma*(1.0+DGrate/2.);

  // double tauH   = (1.0 / G_H);
  // double tauL   = (1.0 / G_L);
  // double tauBar = (1.0 / Gamma);
  
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
  
  double norm = normfactor( par );
  
  //return (1.0 / norm ) * ( (omega*v1) + (1.0-omega)*v2);
  return (1.0 / norm ) * ( (1.0-omega)*v1 + (omega*v2) );
  
  
}

double normfactor(double *par) 
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
  
  double G_H=Gamma*(1.0-DGrate/2.);
  double G_L=Gamma*(1.0+DGrate/2.);
  
  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
  double A0A0   = A1(K0,   G_L, G_H, Gamma, dms, wphase);
  double ApAp   = A2(Kp,   G_L, G_H, Gamma, dms, wphase);
  double AtAt   = A3(Kt,   G_L, G_H, Gamma, dms, wphase);
  //
  double B0B0   = B1(K0,   G_L, G_H, Gamma, dms, wphase);
  double BpBp   = B2(Kp,   G_L, G_H, Gamma, dms, wphase);
  double BtBt   = B3(Kt,   G_L, G_H, Gamma, dms, wphase);
  
  double v1 = A0A0 + ApAp + AtAt;
  double v2 = B0B0 + BpBp + BtBt;
  
  // return factor * ( (omega*v1) + (1.0-omega)*v2 );
  return factor * ( (1.0-omega)*v1 + (omega*v2) );
  
}

/*
  
PDF describing the single observables
Mathematica assisted the analytical integration
Andres Osorio - 14 Aug 2006 / 19 Nov 2006

*/

double properTimePDF(double *x, double *par)
{
  
  //jpsiphi pdf three angles: this is the propertime component
  //the three angles are integrated
  //The MC integral matches this function
  
  double time      = x[0];
  
  if( time < 0 ) return 0.0;
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double scale    = par[10];
  double dms      = xs;
  
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
  
  double norm = normfactor( par );
  
  //return (1.0/norm) * factor * scale * (omega * v1 + (1.0-omega)*v2);
  return scale * (1.0/norm) * factor * ( (1.0-omega)*v1 + (omega*v2) );
  
}


double thetaPDF(double *x, double *par)
{
  
  //jpsiphi pdf three angles: this is the theta component
  //time + two other angles are integrated
  
  double fourothreepi = (4.0/3.0)*TMath::Pi();
  
  double theta    = x[0];
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double scale    = par[10];
  double dms      = xs;
  
  double G_H=Gamma*(1.0-DGrate/2.);
  double G_L=Gamma*(1.0+DGrate/2.);
  
  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;

  scale   = 1.0;
    
  //Bs
  double A0A0   = A1(K0,   G_L, G_H, Gamma, dms, wphase);
  double ApAp   = A2(Kp,   G_L, G_H, Gamma, dms, wphase);
  double AtAt   = A3(Kt,   G_L, G_H, Gamma, dms, wphase);
  
  //Bar (Bs)
  double B0B0   =  B1(K0,   G_L, G_H, Gamma, dms, wphase);
  double BpBp   =  B2(Kp,   G_L, G_H, Gamma, dms, wphase);
  double BtBt   =  B3(Kt,   G_L, G_H, Gamma, dms, wphase);
  
  /////////////////////////////////////////
  //W+
  double v1 = fourothreepi * sin(theta) * ( (ApAp + AtAt)*(3.0+cos(2.0*theta))
                                            + 4.0*A0A0*sin(theta)*sin(theta) );
    
  //W-
  double v2 = fourothreepi * sin(theta) * ( (BpBp + BtBt)*(3.0+cos(2.0*theta))
                                            + 4.0 * B0B0*sin(theta)*sin(theta) );
  
  
  double norm = normfactor( par );
  
  //return -1.0  * (1.0/norm) * scale * ((omega*v1)+(1.0-omega)*v2);
  return scale * (1.0/norm) * (((1.0-omega)*v1)+(omega*v2));
  
}

double psiPDF(double *x, double *par)
{
  
  //jpsiphi pdf three angles: this is the psi component
  //time + two other angles are integrated
  
  double sixteenothreepi = (16.0/3.0)*TMath::Pi(); // (16/3)*pi
  
  double psi      = x[0];
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double scale    = par[10];
  double dms      = xs;
  
  double G_H=Gamma*(1.0-DGrate/2.);
  double G_L=Gamma*(1.0+DGrate/2.);
  
  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
  scale   = 1.0;
  
  //Time depency integrated over
  //Bs
  double A0A0   = A1(K0,   G_L, G_H, Gamma, dms, wphase);
  double ApAp   = A2(Kp,   G_L, G_H, Gamma, dms, wphase);
  double AtAt   = A3(Kt,   G_L, G_H, Gamma, dms, wphase);
  
  //Bar (Bs)
  double B0B0   = B1(K0,   G_L, G_H, Gamma, dms, wphase);
  double BpBp   = B2(Kp,   G_L, G_H, Gamma, dms, wphase);
  double BtBt   = B3(Kt,   G_L, G_H, Gamma, dms, wphase);
  
  //W+
  double v1 = sixteenothreepi * sin(psi)* ( 2.0*A0A0*cos(psi)*cos(psi)
                                            + (ApAp + AtAt)*sin(psi)*sin(psi));
  //W-
  double v2 = sixteenothreepi * sin(psi)* ( 2.0*B0B0*cos(psi)*cos(psi)
                                            + (BpBp + BtBt)*sin(psi)*sin(psi));
  double norm = normfactor( par );

  //return -1.0 * (1.0/norm) * scale * ((omega*v1)+(1.0-omega)*v2);
  return (1.0/norm) * scale * (((1.0-omega)*v1)+(omega*v2));
  
}

double phiPDF(double *x, double *par)
{
  

  //jpsiphi pdf three angles: this is the phi component
  //time + two other angles are integrated
  double sixteenonine = (16.0/9.0);
  
  double phi      = x[0];
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double R0       = par[3];
  
  //par[4] and par[5] free - they where tphase1 tphase2
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double scale    = par[10];
  double dms      = xs;
  
  double G_H=Gamma*(1.0-DGrate/2.);
  double G_L=Gamma*(1.0+DGrate/2.);
  
  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
  double KpKt=sqrt(Kp)*sqrt(Kt);
  
  scale   = 1.0;
  
  //Time depency integrated over
  //Bs
  double A0A0   = A1(K0,   G_L, G_H, Gamma, dms, wphase);
  double ApAp   = A2(Kp,   G_L, G_H, Gamma, dms, wphase);
  double AtAt   = A3(Kt,   G_L, G_H, Gamma, dms, wphase);
  double ImApAt = A5(KpKt, G_L, G_H, Gamma, dms, wphase);
  
  //Bar (Bs)
  double B0B0   = B1(K0,   G_L, G_H, Gamma, dms, wphase);
  double BpBp   = B2(Kp,   G_L, G_H, Gamma, dms, wphase);
  double BtBt   = B3(Kt,   G_L, G_H, Gamma, dms, wphase);
  double ImBpBt = B5(KpKt, G_L, G_H, Gamma, dms, wphase);
  
  //W+
  double twophi = 2.0*phi;
  
  double sin2phi = TMath::Sin(twophi);
  double cos2phi = TMath::Cos(twophi);
  
  double v1 = sixteenonine * (( AtAt - ApAp ) * cos2phi 
                              + 2.0*(A0A0 + ApAp + AtAt 
                                     + ( ImApAt * sin2phi )));
  //W-
  double v2 =  sixteenonine * (( BtBt - BpBp ) * cos2phi
                               + 2.0*(B0B0 + BpBp + BtBt 
                                      + ( ImBpBt * sin2phi )));
  
  
  double norm = normfactor( par );
  
  //return scale * (1.0/norm) * ((omega*v1)+(1.0-omega)*v2);
  return scale * (1.0/norm) * (((1.0-omega)*v1)+(omega*v2));
  
}


///////////////////////////////////////////
//Differential function definitions


double DpsiPDF( double x, void * params)
{
  
  struct pdf_params *p 
    = (struct pdf_params *) params;
  
  double R0     = p->R0;
  double Rt     = p->Rt;
  double wphase = p->phis;
  double Gamma  = p->gBar;
  double DGrate = p->Dgamma;
  double dms    = p->Dms;
  double omega  = p->omega;
  
  double G_H=Gamma*(1.0-DGrate/2.);
  double G_L=Gamma*(1.0+DGrate/2.);
  
  double K0    = R0;
  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
  double A0A0   = A1(K0,   G_L, G_H, Gamma, dms, wphase);
  double ApAp   = A2(Kp,   G_L, G_H, Gamma, dms, wphase);
  double AtAt   = A3(Kt,   G_L, G_H, Gamma, dms, wphase);
  
  double B0B0   = B1(K0,   G_L, G_H, Gamma, dms, wphase);
  double BpBp   = B2(Kp,   G_L, G_H, Gamma, dms, wphase);
  double BtBt   = B3(Kt,   G_L, G_H, Gamma, dms, wphase);

  const double pi = TMath::Pi();
  
  double sixtenothree = (16.0/3.0)*pi;
  
  double dv = 3.0*sixtenothree*(1.0-omega)*cos(x)*sin(x)*sin(x)*(ApAp+AtAt)
    + sixtenothree*omega*sin(x)*(-4.0*B0B0*cos(x)*sin(x)
                                 + 2.0*(BpBp + BtBt)*cos(x)*sin(x))
    + sixtenothree*omega*cos(x)*(2.0*B0B0*cos(x)*cos(x)
                                 + (BpBp + BtBt)*sin(x)*sin(x))
    + A0A0 * ( 2.0*sixtenothree*(1.0-omega)*pow(cos(x),3)
               - 4.0*sixtenothree*(1.0-omega)*cos(x)*sin(x)*sin(x));
  
  return dv;
  
}

double DphiPDF( double x, void * params)
{

  struct pdf_params *p 
    = (struct pdf_params *) params;
  
  double R0     = p->R0;
  double Rt     = p->Rt;
  double wphase = p->phis;
  double Gamma  = p->gBar;
  double DGrate = p->Dgamma;
  double dms    = p->Dms;
  double omega  = p->omega;
  
  double G_H=Gamma*(1.0-DGrate/2.);
  double G_L=Gamma*(1.0+DGrate/2.);

  double Kt    = Rt;
  double Kp    = 1.- R0 -Rt;
  
  double KpKt  = sqrt(Kp)*sqrt(Kt);
  
  double ApAp   = A2(Kp,   G_L, G_H, Gamma, dms, wphase);
  double AtAt   = A3(Kt,   G_L, G_H, Gamma, dms, wphase);
  
  double ImApAt = A5(KpKt, G_L, G_H, Gamma, dms, wphase);
  
  double BpBp   = B2(Kp,   G_L, G_H, Gamma, dms, wphase);
  double BtBt   = B3(Kt,   G_L, G_H, Gamma, dms, wphase);
  double ImBpBt = B5(KpKt, G_L, G_H, Gamma, dms, wphase);
  
  double sixteenonine = (16.0/9.0);
  
  double dv = sixteenonine*(1.0-omega)*(-4.0*ImApAt*cos(2.0*x)
                                        -2.0*(AtAt-ApAp)*sin(2.0*x))
    + sixteenonine*omega*(-4.0*ImBpBt*cos(2.0*x)
                          -2.0*(BtBt-BpBp)*sin(2.0*x));
  
  return dv;
  
}


