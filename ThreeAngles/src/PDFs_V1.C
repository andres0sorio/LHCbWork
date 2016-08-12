// $Id: PDFs_V1.C,v 1.2 2006/11/06 19:22:15 aosorio Exp $
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
  //* the definition is in terms of dcos(theta') dcos(theta")
  
  double sqrtwo = 1.41421356237; // Sqrt[2]
  
  double time      = x[0];
  double theta     = x[1];   // l+ polar angle in Jpsi rest frame
  double psi       = x[2];   // K+ polar angle in Phi rest frame
  double phi       = x[3];   // sum of azimuthal angles
  
  double Gamma    = par[0];
  //double DGrate   = par[1];
  double DGamma   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  double tphase1  = par[4];
  double tphase2  = par[5];
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  //double tmin     = par[9];
  
  //double dms      = xs / Gamma;
  double dms      = xs;
  
  double DGrate=DGamma/Gamma;
  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  //double K0K0=1.0 / ( 1.0 + Rt + Rp);
  //double KtKt=Rt*K0K0;
  //double KpKp=Rp*K0K0;
  
  //double KtKt=Rt;
  //double KpKp=(1-Rt)*Rp;
  //double K0K0=(1-Rt)*(1-Rp);
  
  double KtKt=Rt;
  double KpKp=Rp;
  double K0K0=1 - Rt - Rp;
  
  double K0Kp=sqrt(K0K0*KpKp);
  double K0Kt=sqrt(K0K0*KtKt);
  double KpKt=sqrt(KpKp*KtKt);
  
  //Fixed parameters - to be removed later
  tphase1 = 0.0;
  tphase2 = TMath::Pi();
  
  //Bs
  double A0A0 = 0.5*K0K0*((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			    +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double ApAp = 0.5*KpKp*((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			    +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double AtAt = 0.5*KtKt*((1-cos(wphase))*exp(-G_L*time)+(1+cos(wphase))*exp(-G_H*time)
			    -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  double ReA0Ap = 0.5*K0Kp*cos(tphase2-tphase1)*((1+cos(wphase))*exp(-G_L*time)
						   +(1-cos(wphase))*exp(-G_H*time)
						   +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double ImApAt = KpKt*(exp(-Gamma*time)*(sin(tphase1)*cos(dms*time)
					    -cos(tphase1)*sin(dms*time)*cos(wphase))
			  -0.5*(exp(-G_H*time)-exp(-G_L*time))*cos(tphase1)*sin(wphase));
  double ImA0At = K0Kt*(exp(-Gamma*time)*(sin(tphase2)*cos(dms*time)-cos(tphase2)*sin(dms*time)*cos(wphase))
			  -0.5*(exp(-G_H*time)-exp(-G_L*time))*cos(tphase2)*sin(wphase));
  
  //Bar (Bs)
  double B0B0 = 0.5*K0K0*((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			    -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double BpBp = 0.5*KpKp*((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			    -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double BtBt = 0.5*KtKt*((1-cos(wphase))*exp(-G_L*time)+(1+cos(wphase))*exp(-G_H*time)
			    +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  double ReB0Bp = 0.5*K0Kp*cos(tphase2-tphase1)*((1+cos(wphase))*exp(-G_L*time)
						   +(1-cos(wphase))*exp(-G_H*time)
						   -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double ImBpBt = -KpKt*(exp(-Gamma*time)*(sin(tphase1)*cos(dms*time)
					     -cos(tphase1)*sin(dms*time)*cos(wphase))
			   +0.5*(exp(-G_H*time)-exp(-G_L*time))*cos(tphase1)*sin(wphase));
  double ImB0Bt = -K0Kt*(exp(-Gamma*time)*(sin(tphase2)*cos(dms*time)
					     -cos(tphase2)*sin(dms*time)*cos(wphase))
			   +0.5*(exp(-G_H*time)-exp(-G_L*time))*cos(tphase2)*sin(wphase));
  
  
  /////////////////////////////////////////
  //W+
  double v1 = 4*A0A0*sin(theta)*sin(theta)*cos(psi)*cos(psi)
    + ApAp* ( (1+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      - sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2*phi))
    + AtAt* ( (1+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      + sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2*phi))
    + 2*ImApAt * sin(theta)*sin(theta)*sin(psi)*sin(psi)*sin(2*phi)
    - sqrtwo*ReA0Ap * sin(2*theta) * sin(2*psi) * cos(phi)
    + sqrtwo*ImA0At * sin(2*theta) * sin(2*psi) * sin(phi);
  
  //W-
  double v2 = 4*B0B0*sin(theta)*sin(theta)*cos(psi)*cos(psi)
    + BpBp* ( (1+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      - sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2*phi))
    + BtBt* ( (1+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      + sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2*phi))
    + 2*ImBpBt * sin(theta)*sin(theta)*sin(psi)*sin(psi)*sin(2*phi)
    - sqrtwo*ReB0Bp * sin(2*theta) * sin(2*psi) * cos(phi)
    + sqrtwo*ImB0Bt * sin(2*theta) * sin(2*psi) * sin(phi);
  
  double norm = normfactor( par );
  
  return (1.0 / norm ) * ( (omega*v1) + (1.0-omega)*v2);
  
}


double normfactor(double *par) 
{
  
  double factor =  22.3402144255274173;
  
  double Gamma    = par[0];
  //double DGrate   = par[1];
  double DGamma   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  double tphase1  = par[4];
  double tphase2  = par[5];
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double tmin     = par[9];
  //double scale    = par[10];
  
  //double dms      = xs / Gamma;
  double dms      = xs;
  double DGrate=DGamma/Gamma;
  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  //double K0K0=1.0 / ( 1.0 + Rt + Rp);
  //double KtKt=Rt*K0K0;
  //double KpKp=Rp*K0K0;

  //double KtKt=Rt;
  //double KpKp=(1-Rt)*Rp;
  //double K0K0=(1-Rt)*(1-Rp);

  double KtKt=Rt;
  double KpKp=Rp;
  double K0K0=1 - Rt - Rp;

  //Fixed parameters - to be removed later - A.O.
  tphase1 = 0.0;
  tphase2 = TMath::Pi();
  tmin    = 0.0;
  
  // Bs
  double gammadms = (1.0 / ((Gamma*Gamma) + (dms*dms)) );
  
  double A0A0 = 0.5* K0K0*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			   + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  double ApAp = 0.5* KpKp*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			   + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  double AtAt = 0.5* KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			   - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  // Bar (Bs)
  double B0B0 = 0.5*K0K0*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  double BpBp = 0.5*KpKp*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  double BtBt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  double v1 = A0A0 + ApAp + AtAt;
  
  double v2 = B0B0 + BpBp + BtBt;
  
  return factor * ( (omega*v1) + (1.0-omega)*v2);
  
}


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
  double tphase1  = par[4];
  double tphase2  = par[5];
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];

  //double dms      = xs / Gamma;
  double dms      = xs;

  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);

  //double K0K0=1.0 / ( 1.0 + Rt + Rp);
  //double KtKt=Rt * K0K0;
  //double KpKp=Rp * K0K0;
  
  //double KtKt=Rt;
  //double KpKp=(1-Rt)*Rp;
  //double K0K0=(1-Rt)*(1-Rp);

  double KtKt=Rt;
  double KpKp=Rp;
  double K0K0=1 - Rt - Rp;


  double K0Kp=sqrt(K0K0*KpKp);
  double K0Kt=sqrt(K0K0*KtKt);
  double KpKt=sqrt(KpKp*KtKt);

  //Fixed parameters - to be removed later - A.O.
  tphase1 = 0.0;
  tphase2 = TMath::Pi();
  
  //Bs
  double A0A0 = 0.5* K0K0 *((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			      +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double ApAp = 0.5* KpKp *((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			      +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double AtAt = 0.5* KtKt *((1-cos(wphase))*exp(-G_L*time)+(1+cos(wphase))*exp(-G_H*time)
			      -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  double ReA0Ap = 0.5 * K0Kp * cos(tphase2-tphase1)*((1+cos(wphase))*exp(-G_L*time)
						       +(1-cos(wphase))*exp(-G_H*time)
						       +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double ImApAt = KpKt * (exp(-Gamma*time)*(sin(tphase1)*cos(dms*time)
					      -cos(tphase1)*sin(dms*time)*cos(wphase))
			    - 0.5 *(exp(-G_H*time)-exp(-G_L*time))*cos(tphase1)*sin(wphase));
  double ImA0At = K0Kt * (exp(-Gamma*time)*(sin(tphase2)*cos(dms*time)
					      -cos(tphase2)*sin(dms*time)*cos(wphase))
			    - 0.5 *(exp(-G_H*time)-exp(-G_L*time))*cos(tphase2)*sin(wphase));
  
  //Bar (Bs)
  double B0B0 = 0.5 * K0K0 *((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			       -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double BpBp = 0.5 * KpKp * ((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
				-2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double BtBt = 0.5 * KtKt *((1-cos(wphase))*exp(-G_L*time)+(1+cos(wphase))*exp(-G_H*time)
			       +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  double ReB0Bp = 0.5 * K0Kp *cos(tphase2-tphase1)*((1+cos(wphase))*exp(-G_L*time)
						      +(1-cos(wphase))*exp(-G_H*time)
						      -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  double ImBpBt = -1.0 *KpKt*(exp(-Gamma*time)*(sin(tphase1)*cos(dms*time)
						  -cos(tphase1)*sin(dms*time)*cos(wphase))
				+ 0.5 *(exp(-G_H*time)-exp(-G_L*time))*cos(tphase1)*sin(wphase));
  double ImB0Bt = -1.0 *K0Kt*(exp(-Gamma*time)*(sin(tphase2)*cos(dms*time)
						  - cos(tphase2)*sin(dms*time)*cos(wphase))
				+ 0.5 *(exp(-G_H*time)-exp(-G_L*time))*cos(tphase2)*sin(wphase));
  
  /////////////////////////////////////////
  //W+
  double v1 = 4*A0A0*sin(theta)*sin(theta)*cos(psi)*cos(psi)
    + ApAp* ( (1+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      - sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2*phi))
    + AtAt* ( (1+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      + sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2*phi))
    + 2*ImApAt * sin(theta)*sin(theta)*sin(psi)*sin(psi)*sin(2*phi)
    - sqrtwo*ReA0Ap * sin(2*theta) * sin(2*psi) * cos(phi)
    + sqrtwo*ImA0At * sin(2*theta) * sin(2*psi) * sin(phi);
  
  //W-
  double v2 = 4*B0B0*sin(theta)*sin(theta)*cos(psi)*cos(psi)
    + BpBp* ( (1+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      - sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2*phi))
    + BtBt* ( (1+cos(theta)*cos(theta))*sin(psi)*sin(psi) 
	      + sin(theta)*sin(theta)*sin(psi)*sin(psi)*cos(2*phi))
    + 2*ImBpBt * sin(theta)*sin(theta)*sin(psi)*sin(psi)*sin(2*phi)
    - sqrtwo*ReB0Bp * sin(2*theta) * sin(2*psi) * cos(phi)
    + sqrtwo*ImB0Bt * sin(2*theta) * sin(2*psi) * sin(phi);
  
  double norm = normfactor( par );
  
  
  return (1.0/norm) * ( (omega*v1) + (1-omega)*v2)*intfactor[0]*intfactor[1];
  
}


/*

  PDF describing the single observables
  Mathematica assisted the integration
  Andres Osorio - 14 Aug 2006
  
*/

double properTimePDF(double *x, double *par)
{
  
  //jpsiphi pdf three angles: this is the propertime component
  //the three angles are integrated
  //The MC integral matches this function
  
  //double factor = 0.04476232774; // 9 / (64 pi)
  //double sqrtwo = 1.41421356237; // Sqrt[2]
  
  double time      = x[0];

  if( time < 0 ) return 0.0;
  
  double Gamma    = par[0];
  //double DGrate   = par[1];
  double DGamma   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double scale    = par[10];

  //double dms      = xs / Gamma;
  double dms      = xs;
  double DGrate   = DGamma/Gamma;
  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  //double K0K0= 1.0 / ( 1.0 + Rt + Rp);
  //double KtKt=Rt * K0K0;
  //double KpKp=Rp * K0K0;

  //double KtKt=Rt;
  //double KpKp=(1-Rt)*Rp;
  //double K0K0=(1-Rt)*(1-Rp);
  
  double KtKt=Rt;
  double KpKp=Rp;
  double K0K0=1 - Rt - Rp;
  
  scale   = 1.0;
  
  //Bs
  double A0A0 = 0.5* K0K0 *((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			    +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  double ApAp = 0.5* KpKp *((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			    +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  double AtAt = 0.5* KtKt*((1-cos(wphase))*exp(-G_L*time)+(1+cos(wphase))*exp(-G_H*time)
			   -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  
  //Bar (Bs)
  double B0B0 = 0.5* K0K0*((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			   -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  double BpBp = 0.5* KpKp*((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			   -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  double BtBt = 0.5* KtKt*((1-cos(wphase))*exp(-G_L*time)+(1+cos(wphase))*exp(-G_H*time)
			   +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  
  /////////////////////////////////////////
  //W+
  double v1 = A0A0 + ApAp + AtAt;
  
  //W-
  double v2 = B0B0 + BpBp + BtBt;
  
  double norm = normfactor( par );
  
  double factor =  22.3402144255274173; //Integral of the pdf w.r.t. three angles
  
  return (1.0/norm) * factor * scale * (omega*v1+(1-omega)*v2);
  
  
}


double thetaPDF(double *x, double *par)
{
  
  //jpsiphi pdf three angles: this is the theta component
  //time + two other angles are integrated
  
  double fourothreepi = 4.18879020000; // (4/3)*pi
  
  double theta    = x[0];
  
  double Gamma    = par[0];
  //double DGrate   = par[1];
  double DGamma   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double tmin     = par[9];
  double scale    = par[10];
  
  //double dms = xs / Gamma;
  double dms      = xs;
  double DGrate=DGamma/Gamma;
  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  double KtKt=Rt;
  double KpKp=(1-Rt)*Rp;
  double K0K0=(1-Rt)*(1-Rp);
  
  //double KtKt=Rt;
  //double KpKp=Rp;
  //double K0K0=1 - Rt - Rp;
  
  //double K0K0= 1.0 / ( 1.0 + Rt + Rp);
  //double KtKt=Rt * K0K0;
  //double KpKp=Rp * K0K0;
  
  tmin    = 0.0;
  scale   = 1.0;
  
  //Time depency integrated over
  double gammadms = (1.0 / (Gamma*Gamma + dms*dms) );
  
  double A0A0 = 0.5* K0K0*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			   + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  double ApAp = 0.5* KpKp*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			   + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  double AtAt = 0.5* KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			   - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  
  //Bar (Bs)
  double B0B0 = 0.5*K0K0*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  double BpBp = 0.5*KpKp*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  double BtBt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  /////////////////////////////////////////
  //W+
  double v1 = fourothreepi * ( - 2.0 * A0A0 
				 - 3.0 * ApAp 
			       - 3.0 * AtAt 
			       + 2.0 * A0A0 * cos(2*theta) 
			       - ApAp * cos(2*theta)
				 - AtAt * cos (2*theta) ) * sin(theta);
  
  //W-
  double v2 = fourothreepi * ( - 2.0 * B0B0 
			       - 3.0 * BpBp 
			       - 3.0 * BtBt 
			       + 2.0 * B0B0 * cos(2*theta) 
			       - BpBp * cos(2*theta)
			       - BtBt * cos (2*theta) ) * sin(theta);
  
  double norm = normfactor( par );

  return -1.0  * (1.0/norm) * scale * ((omega*v1)+(1.0-omega)*v2);
  
}

double psiPDF(double *x, double *par)
{
  
  //jpsiphi pdf three angles: this is the psi component
  //time + two other angles are integrated

  double eightothreepi = 8.37758040000; // (8/3)*pi
  
  double psi      = x[0];
  
  double Gamma    = par[0];
  //double DGrate   = par[1];
  double DGamma   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double tmin     = par[9];
  double scale    = par[10];
  
  //double dms = xs / Gamma;
  double dms      = xs;
  double DGrate=DGamma/Gamma;
  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  //double K0K0= 1.0 / ( 1.0 + Rt + Rp);
  //double KtKt=Rt * K0K0;
  //double KpKp=Rp * K0K0;

  //double KtKt=Rt;
  //double KpKp=(1-Rt)*Rp;
  //double K0K0=(1-Rt)*(1-Rp);

  double KtKt=Rt;
  double KpKp=Rp;
  double K0K0=1 - Rt - Rp;

  //Fixed parameters - to be removed later - A.O.
  tmin    = 0.0;
  scale   = 1.0;

  //Time depency integrated over
  double gammadms = (1.0 / (Gamma*Gamma + dms*dms) );
  
  double A0A0 = 0.5*K0K0*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  double ApAp = 0.5*KpKp*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  double AtAt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  //Bar (Bs)
  double B0B0 = 0.5*K0K0*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  double BpBp = 0.5*KpKp*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  double BtBt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  //W+
  double v1 =  - eightothreepi * (  2.0 * A0A0 
				    + ApAp 
				    + AtAt 
				    + 2.0 * A0A0 * cos(2*psi) 
				    - ApAp * cos(2*psi)
				    - AtAt * cos(2*psi) ) * sin(psi);
  
  //W-
  double v2 =  - eightothreepi * ( 2.0 * B0B0 
				   + BpBp 
				   + BtBt 
				   + 2.0 * B0B0 * cos(2*psi) 
				   - BpBp * cos(2*psi)
				   - BtBt * cos(2*psi) ) * sin(psi);
  
  double norm = normfactor( par );
  
  return -1.0 * (1.0/norm) * scale * ((omega*v1)+(1.0-omega)*v2);
  
}

double phiPDF(double *x, double *par)
{
  

  //jpsiphi pdf three angles: this is the phi component
  //time + two other angles are integrated
  double factor2       = 0.02083333333; // (1/48)
  double factor3       = 170.666666666; // 512/3.0;
  double factor4       =  85.333333333; // 256/3.0;
    
  double phi      = x[0];
  
  double Gamma    = par[0];
  //double DGrate   = par[1];
  double DGamma   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  double tphase1  = par[4];
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double tmin     = par[9];
  double scale    = par[10];
  
  //double dms    = xs / Gamma;
  double dms      = xs;
  
  double DGrate=DGamma/Gamma;
  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  //double K0K0= 1.0 / ( 1.0 + Rt + Rp);
  //double KtKt=Rt * K0K0;
  //double KpKp=Rp * K0K0;
  
  //double KtKt=Rt;
  //double KpKp=(1-Rt)*Rp;
  //double K0K0=(1-Rt)*(1-Rp);
  
  double KtKt=Rt;
  double KpKp=Rp;
  double K0K0=1 - Rt - Rp;
  
  double KpKt=sqrt(KpKp*KtKt);
  
  //Fixed parameters - to be removed later - A.O.
  tphase1 = 0.0;
  tmin    = 0.0;
  scale   = 1.0;

  //Time depency integrated over
  double gammadms = (1.0 / (Gamma*Gamma + dms*dms) );
  
  double A0A0 = 0.5*K0K0*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			    +(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
									+dms*cos(dms*tmin))));
  
  double ApAp = 0.5*KpKp*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			    +(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
									+dms*cos(dms*tmin))));
  
  double AtAt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			    +(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
									+dms*cos(dms*tmin))));
  
  
  
  //Bar (Bs)
  double B0B0 = 0.5*K0K0*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			    +(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
									+dms*cos(dms*tmin))));
  
  double BpBp = 0.5*KpKp*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			    +(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
									+dms*cos(dms*tmin))));
  
  double BtBt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			    +(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
									+dms*cos(dms*tmin))));
  
  double ImApAt = KpKt*(gammadms*sin(tphase1)*(exp(-Gamma*tmin)*(Gamma*cos(dms*tmin)
								   -dms*sin(dms*tmin)))
			  - gammadms*cos(tphase1)*cos(wphase)*(exp(-Gamma*tmin)
							       *(Gamma*sin(dms*tmin)
								 +dms*cos(dms*tmin)))
			  - 0.5*((1.0/G_H)*exp(-G_H*tmin)
				 -(1.0/G_L)*exp(-G_L*tmin))*cos(tphase1)*sin(wphase));
  
  double ImBpBt = - KpKt*(gammadms*sin(tphase1)*(exp(-Gamma*tmin)*(Gamma*cos(dms*tmin)
								     -dms*sin(dms*tmin)))
			    - gammadms*cos(tphase1)*cos(wphase)
			    *(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin)))
			    + 0.5*((1.0/G_H)*exp(-G_H*tmin)-(1.0/G_L)
				   *exp(-G_L*tmin))*cos(tphase1)*sin(wphase));

  ///////////////////////////////////////// 
  /* 
     double factor2       = 0.02083333333; // (1/48)
     double factor3       = 170.666666666; // 512/3.0;
     double factor4       =  85.333333333; // 256/3.0; 
  */
  /////////////////////////////////////////
  //W+
  double v1 =  factor2 * (  factor3 * A0A0 
			      + factor3 * ApAp 
			      + factor3 * AtAt 
			      - factor4 * ApAp * cos(2*phi)
			      + factor4 * AtAt * cos(2*phi) 
			      + factor3 * ImApAt * sin(2*phi) );
  
  //W-
  double v2 =  factor2 * (  factor3 * B0B0 
			      + factor3 * BpBp 
			      + factor3 * BtBt 
			      - factor4 * BpBp * cos(2*phi)
			      + factor4 * BtBt * cos(2*phi) 
			      + factor3 * ImBpBt * sin(2*phi) );

  double norm = normfactor( par );
  
  return scale * (1.0/norm) * ((omega*v1)+(1.0-omega)*v2);
  
  
}
