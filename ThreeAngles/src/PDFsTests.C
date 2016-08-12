// $Id: PDFsTests.C,v 1.1 2006/10/18 09:45:18 aosorio Exp $
// Include files
// local
#include "PDFsTests.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PDFs
//
// 2006-06-12 : Andres Osorio Oliveros
//-----------------------------------------------------------------------------


double decayratepdfw(double *x, double *par)
{
  
  //Decay rate only as a fucntion of the wphase
  //Diff. decay rate - definition as a function of t and three angles* + parameters
  //from CERN-TH .... pg.42.
    
  double factor = 0.04476232774; // 9 / (64 pi)
  double sqrtwo = 1.41421356237; // Sqrt[2]
  
  double time      =  1.0;
  double theta     = -0.8;   // cosine of l+ polar angle in Jpsi rest frame [-1;1]
  double psi       =  0.8;   // cosine of K+ polar angle in Phi rest frame [-1;1]
  double phi       =  1.1;   // sum of azimuthal angles
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  double tphase1  = par[4];
  double tphase2  = par[5];
  double dms      = par[6];
  double omega    = par[8];
  
  double wphase   = x[0];
  
  Double_t G_H=Gamma*(1+DGrate/2.);
  Double_t G_L=Gamma*(1-DGrate/2.);
  
  Double_t KtKt=Rt;
  Double_t KpKp=(1-Rt)*Rp;
  Double_t K0K0=(1-Rt)*(1-Rp);
  Double_t K0Kp=sqrt(K0K0*KpKp);
  Double_t K0Kt=sqrt(K0K0*KtKt);
  Double_t KpKt=sqrt(KpKp*KtKt);
  
  Double_t A0A0 = 0.5*K0K0*((1+cos(wphase))*exp(-G_L*time)
			  +(1-cos(wphase))*exp(-G_H*time)
			  +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  Double_t ApAp = KpKp/2*((1+cos(wphase))*exp(-G_L*time)
			  +(1-cos(wphase))*exp(-G_H*time)
			  +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  Double_t AtAt = 0.5*KtKt*((1-cos(wphase))*exp(-G_L*time)
			  +(1+cos(wphase))*exp(-G_H*time)
			  -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  Double_t ReA0Ap = 1/2*K0Kp*cos(tphase2-tphase1)*((1+cos(wphase))*exp(-G_L*time)
						   +(1-cos(wphase))*exp(-G_H*time)
						   +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  Double_t ImApAt = KpKt*(exp(-Gamma*time)*(sin(tphase1)*cos(dms*time)
					    -cos(tphase1)*sin(dms*time)*cos(wphase))
			  - 1/2*(exp(-G_H*time)-exp(-G_L*time))*cos(tphase1)*sin(wphase));
  Double_t ImA0At = K0Kt*(exp(-Gamma*time)*(sin(tphase2)*cos(dms*time)
					    -cos(tphase2)*sin(dms*time)*cos(wphase))
			  - 1/2*(exp(-G_H*time)-exp(-G_L*time))*cos(tphase2)*sin(wphase));
  
  Double_t B0B0 = 0.5*K0K0*((1+cos(wphase))*exp(-G_L*time)
			  +(1-cos(wphase))*exp(-G_H*time)
			  -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  Double_t BpBp = KpKp/2*((1+cos(wphase))*exp(-G_L*time)
			  +(1-cos(wphase))*exp(-G_H*time)
			  -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  Double_t BtBt = 0.5*KtKt*((1-cos(wphase))*exp(-G_L*time)
			  +(1+cos(wphase))*exp(-G_H*time)
			  +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  
  Double_t ReB0Bp = 1/2*K0Kp*cos(tphase2-tphase1)*((1+cos(wphase))*exp(-G_L*time)
						   +(1-cos(wphase))*exp(-G_H*time)
						   -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  Double_t ImBpBt = -KpKt*(exp(-Gamma*time)*(sin(tphase1)*cos(dms*time)
					     -cos(tphase1)*sin(dms*time)*cos(wphase))
			   +1/2*(exp(-G_H*time)-exp(-G_L*time))*cos(tphase1)*sin(wphase));
  Double_t ImB0Bt = -K0Kt*(exp(-Gamma*time)*(sin(tphase2)*cos(dms*time)
					     -cos(tphase2)*sin(dms*time)*cos(wphase))
			   +1/2*(exp(-G_H*time)-exp(-G_L*time))*cos(tphase2)*sin(wphase));
  
  /////////////////////////////////////////
  //W+
  Double_t v1 = 4*A0A0*sin(theta)*sin(theta)*cos(psi)*cos(psi)
    + ApAp* ( (1+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      - sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2*phi))
    + AtAt* ( (1+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      + sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2*phi))
    + 2*ImApAt * sin(theta)*sin(theta)*sin(psi)*sin(psi)*sin(2*phi)
    - sqrtwo*ReA0Ap * sin(2*theta) * sin(2*psi) * cos(phi)
    + sqrtwo*ImA0At * sin(2*theta) * sin(2*psi) * sin(phi);
  
  //W-
  Double_t v2 = 4*B0B0*sin(theta)*sin(theta)*cos(psi)*cos(psi)
    + BpBp* ( (1+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      - sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2*phi))
    + BtBt* ( (1+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      + sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2*phi))
    + 2*ImBpBt * sin(theta)*sin(theta)*sin(psi)*sin(psi)*sin(2*phi)
    - sqrtwo*ReB0Bp * sin(2*theta) * sin(2*psi) * cos(phi)
    + sqrtwo*ImB0Bt * sin(2*theta) * sin(2*psi) * sin(phi);
  
  return factor* (omega*v1 + (1.0 - omega)*v2);
  
}



double jpsiphipdfA(double *x, double *par)
{
  
  /*
    Decay rate definition - only returns Part A 
  */
  
  //double factor = 0.04476232774; // 9 / (64 pi)
  double sqrtwo = 1.41421356237; // Sqrt[2]

  double time      = x[0];
  double theta     = x[1];   // cosine of l+ polar angle in Jpsi rest frame [-1;1]
  double psi       = x[2];   // cosine of K+ polar angle in Phi rest frame [-1;1]
  double phi       = x[3];   // sum of azimuthal angles
  
  double Gamma     = par[0];
  double DGrate    = par[1];
  double Rt        = par[2];
  double Rp        = par[3];
  double tphase1   = par[4];
  double tphase2   = par[5];
  double wphase    = par[6];
  double dms       = par[7];
  
  Double_t G_H     = Gamma*(1+DGrate/2.);
  Double_t G_L     = Gamma*(1-DGrate/2.);
  
  Double_t KtKt    = Rt;
  Double_t KpKp    = (1-Rt)*Rp;
  Double_t K0K0    = (1-Rt)*(1-Rp);
  Double_t K0Kp    = sqrt(K0K0*KpKp);
  Double_t K0Kt    = sqrt(K0K0*KtKt);
  Double_t KpKt    = sqrt(KpKp*KtKt);
  
  Double_t A0A0 = 0.5*K0K0*((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			  +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  Double_t ApAp = 0.5*KpKp*((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			  +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  Double_t AtAt = 0.5*KtKt*((1-cos(wphase))*exp(-G_L*time)+(1+cos(wphase))*exp(-G_H*time)
			  -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  Double_t ReA0Ap = 1/2*K0Kp*cos(tphase2-tphase1)*((1+cos(wphase))*exp(-G_L*time)
						   +(1-cos(wphase))*exp(-G_H*time)
						   +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  Double_t ImApAt = KpKt*(exp(-Gamma*time)*(sin(tphase1)*cos(dms*time)-cos(tphase1)*sin(dms*time)*cos(wphase))-
			  1/2*(exp(-G_H*time)-exp(-G_L*time))*cos(tphase1)*sin(wphase));
  Double_t ImA0At = K0Kt*(exp(-Gamma*time)*(sin(tphase2)*cos(dms*time)-cos(tphase2)*sin(dms*time)*cos(wphase))-
			  1/2*(exp(-G_H*time)-exp(-G_L*time))*cos(tphase2)*sin(wphase));
  
  
    /////////////////////////////////////////
    //W+
  Double_t v1 = 4*A0A0*sin(theta)*sin(theta)*cos(psi)*cos(psi)
    + ApAp* ( (1+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      - sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2*phi))
    + AtAt* ( (1+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      + sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2*phi))
    + 2*ImApAt * sin(theta)*sin(theta)*sin(psi)*sin(psi)*sin(2*phi)
    - sqrtwo*ReA0Ap * sin(2*theta) * sin(2*psi) * cos(phi)
    + sqrtwo*ImA0At * sin(2*theta) * sin(2*psi) * sin(phi);
    
  return v1;
  
}

double jpsiphipdfB(double *x, double *par)
{

  //Decay rate definition - only returns Part B
  //double factor = 0.04476232774; // 9 / (64 pi)
  double sqrtwo = 1.41421356237; // Sqrt[2]

  double time      = x[0];
  double theta     = x[1];   // cosine of l+ polar angle in Jpsi rest frame [-1;1]
  double psi       = x[2];   // cosine of K+ polar angle in Phi rest frame [-1;1]
  double phi       = x[3];   // sum of azimuthal angles
  
  double Gamma = par[0];
  double DGrate = par[1];
  double Rt = par[2];
  double Rp = par[3];
  double tphase1 = par[4];
  double tphase2 = par[5];
  double wphase = par[6];
  double dms = par[7];
  
  Double_t G_H=Gamma*(1+DGrate/2.);
  Double_t G_L=Gamma*(1-DGrate/2.);
  
  Double_t KtKt=Rt;
  Double_t KpKp=(1-Rt)*Rp;
  Double_t K0K0=(1-Rt)*(1-Rp);


  
  Double_t K0Kp=sqrt(K0K0*KpKp);
  Double_t K0Kt=sqrt(K0K0*KtKt);
  Double_t KpKt=sqrt(KpKp*KtKt);

  Double_t B0B0 = 0.5*K0K0*((1+cos(wphase))*exp(-G_L*time)
			  +(1-cos(wphase))*exp(-G_H*time)
			  -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  Double_t BpBp = 0.5*KpKp*((1+cos(wphase))*exp(-G_L*time)
			  +(1-cos(wphase))*exp(-G_H*time)
			  -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  Double_t BtBt = 0.5*KtKt*((1-cos(wphase))*exp(-G_L*time)
			  +(1+cos(wphase))*exp(-G_H*time)
			  +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  
  Double_t ReB0Bp = 0.5*K0Kp*cos(tphase2-tphase1)*((1+cos(wphase))*exp(-G_L*time)
						   +(1-cos(wphase))*exp(-G_H*time)
						   -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  Double_t ImBpBt = -KpKt*(exp(-Gamma*time)*(sin(tphase1)*cos(dms*time)
					     -cos(tphase1)*sin(dms*time)*cos(wphase))
			   +1/2*(exp(-G_H*time)-exp(-G_L*time))*cos(tphase1)*sin(wphase));
  Double_t ImB0Bt = -K0Kt*(exp(-Gamma*time)*(sin(tphase2)*cos(dms*time)
					     -cos(tphase2)*sin(dms*time)*cos(wphase))
			   +1/2*(exp(-G_H*time)-exp(-G_L*time))*cos(tphase2)*sin(wphase));
  
  /////////////////////////////////////////
  //W-
  Double_t v2 = 4*B0B0*sin(theta)*sin(theta)*cos(psi)*cos(psi)
    + BpBp* ( (1+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      - sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2*phi))
    + BtBt* ( (1+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      + sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2*phi))
    + 2*ImBpBt * sin(theta)*sin(theta)*sin(psi)*sin(psi)*sin(2*phi)
    - sqrtwo*ReB0Bp * sin(2*theta) * sin(2*psi) * cos(phi)
    + sqrtwo*ImB0Bt * sin(2*theta) * sin(2*psi) * sin(phi);
  
  
  return v2;
  
}

