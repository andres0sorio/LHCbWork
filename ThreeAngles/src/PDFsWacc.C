// $Id: PDFsWacc.C,v 1.6 2006/12/09 18:17:55 aosorio Exp $
// Include files 
#include <Math/Integrator.h>
#include <Math/IntegrationTypes.h>
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/WrappedFunction.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

// local
#include "PDFsWacc.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PDFsWacc
//
// 2006-11-18 : Andres Felipe OSORIO OLIVEROS
//-----------------------------------------------------------------------------

double TimeAcceptance  ( double time , double *par )
{

  double KK  = 1.1;
  double arg = KK*time;
  return (pow( arg,3.0))/(1.0 + pow( arg,3.0));
  
}

double pTimeAcc     ( double * xx , double *par )
{
  
  double tt = xx[0];

  double value = TimeAcceptance(tt,par) * properTime (xx ,par);
  
  return value;
  
}

double pTimeAccPDF (  double * xx , double *par )
{
  
  double value = pTimeAcc( xx , par ) ;
  
  //Normalisation factor:  
  double nf = accnormfactor( par );
  
  return (1.0 / nf) * value;
  
}


double jpsiphiAccPDF   ( double *xx, double *par )
{
  
  double tt = xx[0];
  
  double value = TimeAcceptance(tt,par) * jpsiphi(xx ,par);
  
  //Normalisation factor:
  //
  //
  
  return value;
  
}

//--- theta

// double thetaAcc(double *x, double *par)
// {
  
//   //jpsiphi pdf three angles: this is the theta component
//   //time + two other angles are integrated
  
//   double fourothreepi = (4.0/3.0)*TMath::Pi();
  
//   double theta    = x[0];
  
//   double Gamma    = par[0];
//   double DGamma   = par[1];
//   double ft       = par[2];
//   double R0       = par[3];
//   //par[4] and par[5] free - they where tphase1 tphase2
//   double wphase   = par[6];
//   double xs       = par[7];
//   double omega    = par[8];
//   double scale    = par[10];
//   double accK1    = par[11];
//   double dms      = xs;
  
//   double DGrate=DGamma/Gamma;
//   double G_H=Gamma*(1.0-DGrate/2.);
//   double G_L=Gamma*(1.0+DGrate/2.);

//   double K0    = R0;
//   double Kt    = ft;
//   double Kp    = 1.- R0 -ft;
  
//   scale   = 1.0;
    
//   //Bs
//   double A0A0   = A1(K0,   G_L, G_H, Gamma, dms, wphase, accK1);
//   double ApAp   = A2(Kp,   G_L, G_H, Gamma, dms, wphase, accK1);
//   double AtAt   = A3(Kt,   G_L, G_H, Gamma, dms, wphase, accK1);
//   //Bar (Bs)
//   double B0B0   =  B1(K0,   G_L, G_H, Gamma, dms, wphase, accK1);
//   double BpBp   =  B2(Kp,   G_L, G_H, Gamma, dms, wphase, accK1);
//   double BtBt   =  B3(Kt,   G_L, G_H, Gamma, dms, wphase, accK1);
  
//   /////////////////////////////////////////
//   //W+
//   double v1 = fourothreepi * sin(theta) * ( (ApAp + AtAt)*(3.0+cos(2.0*theta))
//                                             + 4.0*A0A0*sin(theta)*sin(theta) );
//   //W-
//   double v2 = fourothreepi * sin(theta) * ( (BpBp + BtBt)*(3.0+cos(2.0*theta))
//                                             + 4.0 * B0B0*sin(theta)*sin(theta) );
  
//   return scale * (1.0/norm) * (((1.0-omega)*v1)+(omega*v2));
  
// }

// //--- psi

// double psiAcc(double *x, double *par)
// {
  
//   //jpsiphi pdf three angles: this is the psi component
//   //time + two other angles are integrated
  
//   double sixteenothreepi = (16.0/3.0)*TMath::Pi(); // (16/3)*pi
  
//   double psi      = x[0];
//   double Gamma    = par[0];
//   double DGamma   = par[1];
//   double ft       = par[2];
//   double R0       = par[3];
//   //par[4] and par[5] free - they where tphase1 tphase2
//   double wphase   = par[6];
//   double xs       = par[7];
//   double omega    = par[8];
//   double scale    = par[10];
//   double accK1    = par[11];
//   double dms      = xs;
  
//   double DGrate=DGamma/Gamma;
//   double G_H=Gamma*(1.0-DGrate/2.);
//   double G_L=Gamma*(1.0+DGrate/2.);

//   double K0    = R0;
//   double Kt    = ft;
//   double Kp    = 1.- R0 -Rt;
  
//   scale   = 1.0;
  
//   //Time depency integrated over
//   //Bs
//   double A0A0   = A1(K0,   G_L, G_H, Gamma, dms, wphase, accK1);
//   double ApAp   = A2(Kp,   G_L, G_H, Gamma, dms, wphase, accK1);
//   double AtAt   = A3(Kt,   G_L, G_H, Gamma, dms, wphase, accK1);
//   //Bar (Bs)
//   double B0B0   = B1(K0,   G_L, G_H, Gamma, dms, wphase, accK1);
//   double BpBp   = B2(Kp,   G_L, G_H, Gamma, dms, wphase, accK1);
//   double BtBt   = B3(Kt,   G_L, G_H, Gamma, dms, wphase, accK1);
  
//   //W+
//   double v1 = sixteenothreepi * sin(psi)* ( 2.0*A0A0*cos(psi)*cos(psi)
//                                             + (ApAp + AtAt)*sin(psi)*sin(psi));
//   //W-
//   double v2 =  sixteenothreepi * sin(psi)* ( 2.0*B0B0*cos(psi)*cos(psi)
//                                              + (BpBp + BtBt)*sin(psi)*sin(psi));
  
//   return scale * (((1.0-omega)*v1)+(omega*v2));
  
// }

// //--- phi
// double phiPDF(double *x, double *par)
// {
  
//   //jpsiphi pdf three angles: this is the phi component
//   //time + two other angles are integrated
//   double sixteenonine = (16.0/9.0);
  
//   double phi      = x[0];
//   double Gamma    = par[0];
//   double DGamma   = par[1];
//   double Rt       = par[2];
//   double R0       = par[3];
//   //par[4] and par[5] free - they where tphase1 tphase2
//   double wphase   = par[6];
//   double xs       = par[7];
//   double omega    = par[8];
//   double scale    = par[10];
//   double accK1    = par[11];
  
//   double dms      = xs;
  
//   double DGrate=DGamma/Gamma;
//   double G_H=Gamma*(1.0-DGrate/2.);
//   double G_L=Gamma*(1.0+DGrate/2.);

//   double K0    = R0;
//   double Kt    = Rt;
//   double Kp    = 1.- R0 -Rt;
  
//   double KpKt=sqrt(Kp)*sqrt(Kt);
  
//   scale   = 1.0;
  
//   //Time depency integrated over
//   //Bs
//   double A0A0   = A1(K0,   G_L, G_H, Gamma, dms, wphase, accK1);
//   double ApAp   = A2(Kp,   G_L, G_H, Gamma, dms, wphase, accK1);
//   double AtAt   = A3(Kt,   G_L, G_H, Gamma, dms, wphase, accK1);
//   double ImApAt = A5(KpKt, G_L, G_H, Gamma, dms, wphase, accK1);
//   //Bar (Bs)
//   double B0B0   = B1(K0,   G_L, G_H, Gamma, dms, wphase, accK1);
//   double BpBp   = B2(Kp,   G_L, G_H, Gamma, dms, wphase, accK1);
//   double BtBt   = B3(Kt,   G_L, G_H, Gamma, dms, wphase, accK1);
//   double ImBpBt = B5(KpKt, G_L, G_H, Gamma, dms, wphase, accK1);
  
//   //W+
//   double twophi = 2.0*phi;
  
//   double sin2phi = TMath::Sin(twophi);
//   double cos2phi = TMath::Cos(twophi);
  
//   double v1 = sixteenonine * (( AtAt - ApAp ) * cos2phi 
//                               + 2.0*(A0A0 + ApAp + AtAt 
//                                      - ( ImApAt * sin2phi )));
//   //W-
//   double v2 =  sixteenonine * (( BtBt - BpBp ) * cos2phi
//                                + 2.0*(B0B0 + BpBp + BtBt 
//                                       - ( ImBpBt * sin2phi )));
  
//   return scale * (1.0/norm) * (((1.0-omega)*v1)+(omega*v2));
  
// }


//--------------------------------------------------
/* Normalisation - auxiliary functions definition

   Using MathMore Library

*/
//...............  Numerical integration of pTimeAcc


double accnormfactor( double * pars )
{

  pdf_params * pms = new pdf_params();

  pms->gBar   = pars[0];
  pms->Dgamma = pars[1];
  pms->ft     = pars[2];
  pms->fp     = pars[3];
  pms->dt1    = 0.0;
  pms->dt2    = 0.0;
  pms->phis   = pars[6];
  pms->Dms    = pars[7];
  pms->omega  = pars[8];

  ROOT::Math::Integrator::GSLFuncPointer ff = &pTimeAcc;
  
  ROOT::Math::Integrator * nminteg = new ROOT::Math::Integrator(ff,ROOT::Math::Integration::ADAPTIVE, 1.E-9, 1E-6, 1000 );
  
  //--method for integration
  nminteg->SetIntegrationRule(ROOT::Math::Integration::GAUSS21);

  //--range [a,inf]
  double eval =  nminteg->IntegralUp(ff,pms,0.0);

  delete nminteg;
  delete pms;
  
  
  return eval;
  
}


double pTimeAcc (double x, void * params) {
  
  struct pdf_params *p 
    = (struct pdf_params *) params;
  
  double par[20];
  par[0]  = p->gBar;
  par[1]  = p->Dgamma;
  par[2]  = p->ft;
  par[3]  = p->fp;
  par[4]  = 0.0;
  par[5]  = 0.0;
  par[6]  = p->phis;
  par[7]  = p->Dms;
  par[8]  = p->omega;
  par[9]  = 1.0;
  par[10] = 1.0;
  par[11] = 1.0;
  
  double v1 = pTimeAcc( &x , par);
  
  return v1;
  
}
