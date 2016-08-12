// $Id: DiffPDFs.C,v 1.1 2006/12/09 18:19:18 aosorio Exp $
// Include files 


// local
#include "PDFs.h"
//-----------------------------------------------------------------------------
// Implementation file for class : DiffPDFs
//
// 2006-12-09 : Andres Osorio
//-----------------------------------------------------------------------------

///////////////////////////////////////////
//Differential function definitions

double DthetaWpPDF( double x, void * params)
{
  
  struct pdf_params *p 
    = (struct pdf_params *) params;
  
  const double pi     = TMath::Pi();
  double fourothreepi = (4.0/3.0) * pi;
  
  double fp     = p->fp;
  double ft     = p->ft;
  double wphase = p->phis;
  double Gamma  = p->gBar;
  double DGrate = p->Dgamma;
  double dms    = p->Dms;
    
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
  
  double dv     = fourothreepi *cos(x)*(6.0*A0A0 + ApAp + AtAt 
					+ 3.0*( AtAt + ApAp - 2.0*A0A0 )*cos(2*x));
  
  return dv;
  
}

double DthetaWmPDF( double x, void * params)
{
  
  struct pdf_params *p 
    = (struct pdf_params *) params;
  
  const double pi = TMath::Pi();
  double fourothreepi = (4.0/3.0) * pi;
  
  double fp     = p->fp;
  double ft     = p->ft;
  double wphase = p->phis;
  double Gamma  = p->gBar;
  double DGrate = p->Dgamma;
  double dms    = p->Dms;

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
  
  double dv =  fourothreepi *cos(x)*(6.0*B0B0 + BpBp + BtBt 
				     + 3.0*( BtBt + BpBp - 2.0*B0B0 )*cos(2*x));
				     
  return dv;
				     
}

double DpsiWpPDF( double x, void * params)
{
  
  struct pdf_params *p 
    = (struct pdf_params *) params;

  const double pi = TMath::Pi();
  double eightothree = (8.0/3.0);
  
  double fp     = p->fp;
  double ft     = p->ft;
  double wphase = p->phis;
  double Gamma  = p->gBar;
  double DGrate = p->Dgamma;
  double dms    = p->Dms;
  
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
  
  double fp     = p->fp;
  double ft     = p->ft;
  double wphase = p->phis;
  double Gamma  = p->gBar;
  double DGrate = p->Dgamma;
  double dms    = p->Dms;
  
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
  
  double fp     = p->fp;
  double ft     = p->ft;
  double wphase = p->phis;
  double Gamma  = p->gBar;
  double DGrate = p->Dgamma;
  double dms    = p->Dms;
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);

  double KtKt  = Kt(ft,fp);
  double KpKp  = Kp(ft,fp);
  //double K0K0  = K0(ft,fp);

  double KpKt  = sqrt(KpKp)*sqrt(KtKt);
  
  double ApAp   = A2def(KpKp,   tauL, tauH, tauBar, dms, wphase);
  double AtAt   = A3def(KtKt,   tauL, tauH, tauBar, dms, wphase);
  double ImApAt = A5def(KpKt, tauL, tauH, tauBar, dms, wphase);

  double dv = sixteenonine * ( 4.0*ImApAt*cos(2.0*x)
                              -2.0*(AtAt-ApAp)*sin(2.0*x));
  return dv;
  
}

double DphiWmPDF( double x, void * params)
{
  
  struct pdf_params *p 
    = (struct pdf_params *) params;
  
  double sixteenonine = (16.0/9.0);

  double fp     = p->fp;
  double ft     = p->ft;
  double wphase = p->phis;
  double Gamma  = p->gBar;
  double DGrate = p->Dgamma;
  double dms    = p->Dms;
  
  double G_H    = GHeavy(Gamma,DGrate);
  double G_L    = GLight(Gamma,DGrate);
  
  double tauH   = (1.0 / G_H);
  double tauL   = (1.0 / G_L);
  double tauBar = (1.0 / Gamma);
  
  double KtKt  = Kt(ft,fp);
  double KpKp  = Kp(ft,fp);
  //double K0K0  = K0(ft,fp);

  double KpKt  = sqrt(KpKp)*sqrt(KtKt);

  double BpBp   = B2def(KpKp,   tauL, tauH, tauBar, dms, wphase);
  double BtBt   = B3def(KtKt,   tauL, tauH, tauBar, dms, wphase);
  double ImBpBt = B5def(KpKt, tauL, tauH, tauBar, dms, wphase);
  
  double dv = sixteenonine * ( 4.0*ImBpBt*cos(2.0*x)
                              -2.0*(BtBt-BpBp)*sin(2.0*x));
  
  return dv;
  
}

