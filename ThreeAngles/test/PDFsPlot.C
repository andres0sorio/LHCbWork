// $Id: PDFsPlot.C,v 1.1 2006/11/16 21:25:41 aosorio Exp $
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
  
  double sqrtwo = 1.41421356237; // Sqrt[2]
  double intfactor[2];
  intfactor[0]=1.0;
  intfactor[1]=1.0;
  
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

  intfactor[0] = -1.*TMath::Sin(theta);
  intfactor[1] = -1.*TMath::Sin(psi);
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  
  //double dms      = xs / Gamma;
  double dms      = xs;
  
  double G_H=Gamma*(1.0-DGrate/2.);
  double G_L=Gamma*(1.0+DGrate/2.);
  
  double Kt  = Rt;
  double Kp  = Rp;
  double K0  = Rp;

  double K0Kp=sqrt(K0*Kp);
  double K0Kt=sqrt(K0*Kt);
  double KpKt=sqrt(Kp*Kt);
  
  //Fixed parameters - to be removed later - A.O.
  //tphase1 = 0.0;
  //tphase2 = TMath::Pi();
  
  //Bs
  double A0A0 = 0.5*K0*((1.0+cos(wphase))*exp(-G_L*time)+(1.0-cos(wphase))*exp(-G_H*time)
			      +2.0*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double ApAp = 0.5*Kp*((1.0+cos(wphase))*exp(-G_L*time)+(1.0-cos(wphase))*exp(-G_H*time)
			      +2.0*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double AtAt = 0.5*Kt*((1.0-cos(wphase))*exp(-G_L*time)+(1.0+cos(wphase))*exp(-G_H*time)
			      -2.0*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  double ReA0Ap = 0.5*K0Kp*cos(tphase2-tphase1)*((1.0+cos(wphase))*exp(-G_L*time)
						       +(1.0-cos(wphase))*exp(-G_H*time)
						       +2.0*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double ImApAt = KpKt * (exp(-Gamma*time)*(sin(tphase1)*cos(dms*time)
					      -cos(tphase1)*sin(dms*time)*cos(wphase))
			    - 0.5 *(exp(-G_H*time)-exp(-G_L*time))*cos(tphase1)*sin(wphase));
  double ImA0At = K0Kt * (exp(-Gamma*time)*(sin(tphase2)*cos(dms*time)
					      -cos(tphase2)*sin(dms*time)*cos(wphase))
			    - 0.5 *(exp(-G_H*time)-exp(-G_L*time))*cos(tphase2)*sin(wphase));
  
  //Bar (Bs)
  double B0B0 = 0.5*K0*((1.0+cos(wphase))*exp(-G_L*time)+(1.0-cos(wphase))*exp(-G_H*time)
			       -2.0*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double BpBp = 0.5 *Kp* ((1.0+cos(wphase))*exp(-G_L*time)+(1.0-cos(wphase))*exp(-G_H*time)
				-2.0*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double BtBt = 0.5*Kt*((1.0-cos(wphase))*exp(-G_L*time)+(1.0+cos(wphase))*exp(-G_H*time)
			       +2.0*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  double ReB0Bp = 0.5*K0Kp*cos(tphase2-tphase1)*((1.0+cos(wphase))*exp(-G_L*time)
						      +(1.0-cos(wphase))*exp(-G_H*time)
						      -2.0*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double ImBpBt = -1.0*KpKt*(exp(-Gamma*time)*(sin(tphase1)*cos(dms*time)
						  -cos(tphase1)*sin(dms*time)*cos(wphase))
				+ 0.5 *(exp(-G_H*time)-exp(-G_L*time))*cos(tphase1)*sin(wphase));
  double ImB0Bt = -1.0*K0Kt*(exp(-Gamma*time)*(sin(tphase2)*cos(dms*time)
						  - cos(tphase2)*sin(dms*time)*cos(wphase))
				+ 0.5 *(exp(-G_H*time)-exp(-G_L*time))*cos(tphase2)*sin(wphase));
  
  /////////////////////////////////////////
  //W+
  double v1 = 4.0*A0A0*sin(theta)*sin(theta)*cos(psi)*cos(psi)
    + ApAp* ( (1.0+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      - sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2.0*phi))
    + AtAt* ( (1.0+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      + sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2.0*phi))
    - 2.0*ImApAt * sin(theta)*sin(theta)*sin(psi)*sin(psi)*sin(2.0*phi)
    + sqrtwo*ReA0Ap * sin(2.0*theta) * sin(2.0*psi) * cos(phi)
    + sqrtwo*ImA0At * sin(2.0*theta) * sin(2.0*psi) * sin(phi);
  
  //W-
  double v2 = 4.0*B0B0*sin(theta)*sin(theta)*cos(psi)*cos(psi)
    + BpBp* ( (1.0+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      - sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2.0*phi))
    + BtBt* ( (1.0+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      + sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2.0*phi))
    - 2.0*ImBpBt * sin(theta)*sin(theta)*sin(psi)*sin(psi)*sin(2.0*phi)
    + sqrtwo*ReB0Bp * sin(2.0*theta) * sin(2.0*psi) * cos(phi)
    + sqrtwo*ImB0Bt * sin(2.0*theta) * sin(2.0*psi) * sin(phi);
  
  double norm = normfactor( par );
  
  
  return (1.0/norm) * ( (omega*v1) + (1.0-omega)*v2)*intfactor[0]*intfactor[1];
  
}

