// $Id: EvenOdd_PDFs.C,v 1.2 2006/11/06 19:22:15 aosorio Exp $
// Include files
// local
#include "EvenOdd_PDFs.h"

//-----------------------------------------------------------------------------
// Implementation file for class : EvenOdd_PDFs
//
// 2006-10-18 : Andres Osorio Oliveros
//-----------------------------------------------------------------------------

double properTimePDFOdd(double *x, double *par)
{
  
  //jpsiphi pdf three angles: this is the propertime component
  //the three angles are integrated
  //The MC integral matches this function
  
  //double factor = 0.04476232774; // 9 / (64 pi)
  //double sqrtwo = 1.41421356237; // Sqrt[2]
  
  double time      = x[0];

  if( time < 0 ) return 0.0;
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  double wphase   = par[6];
  double xs      = par[7];
  double omega    = par[8];
  //double tmin     = par[9];
  double scale    = par[10];
  //double scale    = 1.0;

  double dms = xs / Gamma;

  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  double K0K0=1.0 / (1 + Rt + Rp);
  double KtKt=Rt * K0K0;
  
  //Bs
  double AtAt = 0.5* KtKt*((1-cos(wphase))*exp(-G_L*time)+(1+cos(wphase))*exp(-G_H*time)
			     -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  //Bar (Bs)
  double BtBt = 0.5* KtKt*((1-cos(wphase))*exp(-G_L*time)+(1+cos(wphase))*exp(-G_H*time)
			     +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  /////////////////////////////////////////
  //W+
  double v1 = AtAt;
  
  //W-
  double v2 = BtBt;
  
  
  return scale * ((omega*v1)+(1-omega)*v2);
  
}

double properTimePDFEven(double *x, double *par)
{
  
  //jpsiphi pdf three angles: this is the propertime component
  //the three angles are integrated
  //The MC integral matches this function
  
  //double factor = 0.04476232774; // 9 / (64 pi)
  //double sqrtwo = 1.41421356237; // Sqrt[2]
  
  double time      = x[0];

  if( time < 0 ) return 0.0;
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  double wphase   = par[6];
  double xs      = par[7];
  double omega    = par[8];
  //double tmin     = par[9];
  double scale    = par[10];
  //double scale    = 1.0;

  double dms = xs / Gamma;

  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  double K0K0= 1.0 / ( 1.0 + Rt + Rp);
  double KpKp=Rp * K0K0;
  
  //Bs
  double A0A0 = 0.5* K0K0 *((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			      +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  double ApAp = 0.5* KpKp *((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			      +2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  //Bar (Bs)
  double B0B0 = 0.5* K0K0*((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			     -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  double BpBp = 0.5* KpKp*((1+cos(wphase))*exp(-G_L*time)+(1-cos(wphase))*exp(-G_H*time)
			     -2*exp(-Gamma*time)*sin(dms*time)*sin(wphase));
  
  /////////////////////////////////////////
  //W+
  double v1 = A0A0 + ApAp;
  
  //W-
  double v2 = B0B0 + BpBp;
  
  
  return scale * ((omega*v1)+(1-omega)*v2);
  
}

double thetaEven(double *x, double *par)
{
  
  //jpsiphi pdf three angles: this is the theta component
  //time + two other angles are integrated
  
  //double sqrtwo       = 1.41421356237; // Sqrt[2]
  double factor       = 0.04476232774; // 9 / (64 pi)
  double fourothreepi = 4.18879020000; // (4/3)*pi
  
  double theta    = x[0];
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  double wphase   = par[6];
  double xs      = par[7];
  double omega    = par[8];
  double tmin     = par[9];
  double scale    = par[10];

  double dms = xs / Gamma;
  
  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  //double KtKt=Rt;
  double KpKp=(1-Rt)*Rp;
  double K0K0=(1-Rt)*(1-Rp);

  //double K0Kp=sqrt(K0K0*KpKp);
  //double K0Kt=sqrt(K0K0*KtKt);
  //double KpKt=sqrt(KpKp*KtKt);
  
  //Time depency integrated over
  double gammadms = ( 1.0 / (Gamma*Gamma + dms*dms) );
  
  double A0A0 = 0.5*K0K0*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			  +(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
								      +dms*cos(dms*tmin))));
  
  double ApAp = 0.5*KpKp*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			  +(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
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
  
  /////////////////////////////////////////
  //W+
  double v1 = fourothreepi * ( - 2.0 * A0A0 
				 - 3.0 * ApAp 
				 + 2.0 * A0A0 * cos(2*theta) 
				 - ApAp * cos(2*theta) ) * sin(theta);
  
  //W-
  double v2 = fourothreepi * ( - 2.0 * B0B0 
				 - 3.0 * BpBp 
				 + 2.0 * B0B0 * cos(2*theta) 
				 - BpBp * cos(2*theta) ) * sin(theta);
  
  
  return -1.0*scale*factor * ((omega*v1)+(1.0-omega)*v2);
  
}

double thetaOdd(double *x, double *par)
{
  
  //jpsiphi pdf three angles: this is the theta component
  //time + two other angles are integrated

  //double sqrtwo       = 1.41421356237; // Sqrt[2]
  double factor       = 0.04476232774; // 9 / (64 pi)
  double fourothreepi = 4.18879020000; // (4/3)*pi
  
  double theta    = x[0];
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  //double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  double wphase   = par[6];
  double xs      = par[7];
  double omega    = par[8];
  double tmin     = par[9];
  double scale    = par[10];
  
  double dms = xs / Gamma;

  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  double KtKt=Rt;
  //double KpKp=(1-Rt)*Rp;
  //double K0K0=(1-Rt)*(1-Rp);

  //double K0Kp=sqrt(K0K0*KpKp);
  //double K0Kt=sqrt(K0K0*KtKt);
  //double KpKt=sqrt(KpKp*KtKt);

  //Time depency integrated over
  double gammadms = (1.0 / (Gamma*Gamma + dms*dms) );
  
  double AtAt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			  +(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
								      +dms*cos(dms*tmin))));
  
  
  //Bar (Bs)
  double BtBt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			  +(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
								      +dms*cos(dms*tmin))));
  
  
  /////////////////////////////////////////
  //W+
  double v1 = fourothreepi * (- 3.0 * AtAt 
				- AtAt * cos (2*theta) ) * sin(theta);
  
  //W-
  double v2 = fourothreepi * ( - 3.0 * BtBt 
				 - BtBt * cos (2*theta) ) * sin(theta);
  
  
  return -1.0*scale*factor * ((omega*v1)+(1.0-omega)*v2);
  
}

////////////////////////////////
//

double psiPDFEven(double *x, double *par)
{
  
  //jpsiphi pdf three angles: this is the psi component
  //time + two other angles are integrated
  double factor        = 0.04476232774; // 9 / (64 pi)
  double eightothreepi = 8.37758040000; // (8/3)*pi
  
  double psi      = x[0];
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  double wphase   = par[6];
  double xs      = par[7];
  double omega    = par[8];
  double tmin     = par[9];
  double scale    = par[10];
  
  double dms = xs / Gamma;

  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  double KpKp=(1-Rt)*Rp;
  double K0K0=(1-Rt)*(1-Rp);
  
  //Time depency integrated over
  double gammadms = (1.0 / (Gamma*Gamma + dms*dms) );
  
  double A0A0 = 0.5*K0K0*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  double ApAp = 0.5*KpKp*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  //Bar (Bs)
  double B0B0 = 0.5*K0K0*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  double BpBp = 0.5*KpKp*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  //W+
  double v1 =  - eightothreepi * (  2.0 * A0A0 
				      + ApAp 
				      + 2.0 * A0A0 * cos(2*psi) 
				      - ApAp * cos(2*psi) ) * sin(psi);
  
  //W-
  double v2 =  - eightothreepi * ( 2.0 * B0B0 
				     + BpBp 
				     + 2.0 * B0B0 * cos(2*psi) 
				     - BpBp * cos(2*psi) ) * sin(psi);
  
  return -1.0 * scale * factor * ((omega*v1)+(1.0-omega)*v2);
  
}

double psiPDFOdd(double *x, double *par)
{
  
  //jpsiphi pdf three angles: this is the psi component
  //time + two other angles are integrated

  //double sqrtwo        = 1.41421356237; // Sqrt[2]
  double factor        = 0.04476232774; // 9 / (64 pi)
  double eightothreepi = 8.37758040000; // (8/3)*pi
  
  double psi      = x[0];
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  //double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  double wphase   = par[6];
  double xs      = par[7];
  double omega    = par[8];
  double tmin     = par[9];
  double scale    = par[10];
  
  double dms = xs / Gamma;

  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  double KtKt=Rt;
   
  //Time depency integrated over
  double gammadms = (1.0 / (Gamma*Gamma + dms*dms) );
  
  
  double AtAt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  //Bar (Bs)
  						  
  double BtBt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  //W+
  double v1 =  - eightothreepi * (  AtAt 
				    - AtAt * cos(2*psi) ) * sin(psi);
  
  //W-
  double v2 =  - eightothreepi * ( BtBt 
				   - BtBt * cos(2*psi) ) * sin(psi);
  
  
  return -1.0 * scale * factor * ((omega*v1)+(1.0-omega)*v2);
  
}


double SinpsiEven(double *x, double *par) {

  //jpsiphi pdf three angles: this is the theta component
  //time + two other angles are integrated
  
  double factor        = 0.04476232774;  // 9 / (64 pi)

  double eightothreepi = 8.37758040000;  // (8/3)*pi
  double sintheta    = x[0];
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  double wphase   = par[6];
  double xs      = par[7];
  double omega    = par[8];
  double tmin     = par[9];
  //double scale    = par[10];
  
  double dms = xs / Gamma;

  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  double KpKp=(1-Rt)*Rp;
  double K0K0=(1-Rt)*(1-Rp);

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
  
  //Bar (Bs)
  double B0B0 = 0.5*K0K0*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			  +(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
								      +dms*cos(dms*tmin))));
  
  double BpBp = 0.5*KpKp*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			  +(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
								      +dms*cos(dms*tmin))));
  
  /////////////////////////////////////////
  //W+
  double v1 = - eightothreepi * ( (2.0 * A0A0) + ApAp + 2.0 * A0A0 * ( 1 - 2.0*sintheta*sintheta )
				  - ApAp * ( 1 - 2.0*sintheta*sintheta ) );
  
  //W-
  double v2 = - eightothreepi * ( (2.0 * B0B0) + BpBp + 2.0 * B0B0 * ( 1 - 2.0*sintheta*sintheta )
				  - BpBp * ( 1 - 2.0*sintheta*sintheta ) );
  
  
  
  return  factor * ((omega*v1)+(1.0-omega)*v2);
  
}

double SinpsiOdd(double *x, double *par)
{
  
  //jpsiphi pdf three angles: this is the theta component
  //time + two other angles are integrated
  
  double factor       = 0.04476232774;  // 9 / (64 pi)
  
  double eightothreepi = 8.37758040000; // (8/3)*pi
  
  double sintheta    = x[0];
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  //double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double tmin     = par[9];
  //double scale    = par[10];
  
  double dms = xs / Gamma;

  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  double KtKt=Rt;
  
  //Time depency integrated over
  double gammadms = (1.0 / (Gamma*Gamma + dms*dms) );
  
  double AtAt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			    +(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
									+dms*cos(dms*tmin))));
  
  
  //Bar (Bs)
  double BtBt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			    +(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
									+dms*cos(dms*tmin))));
  
  
  /////////////////////////////////////////
  //W+
  double v1 = - eightothreepi * ( AtAt - AtAt * ( 1.0 - 2.0*sintheta * sintheta) );
  
  //W-
  double v2 = - eightothreepi * ( BtBt - BtBt * ( 1.0 - 2.0*sintheta * sintheta) );
  
  
  return  factor * ((omega*v1)+(1.0-omega)*v2);
  
}


double CosthetaTrEven(double *x, double *par)
{
  
  //Cosine of theta_Tr 
  // angle relations used to transformed for helicity basis to transversity basis
  double scale  = par[10];

  double result = scale * SinpsiEven(x,par) * SinphiEven(x,par);
  
  return result;
  
}

double CosthetaTrOdd(double *x, double *par)
{
  
  //Cosine of theta_Tr 
  // angle relations used to transformed for helicity basis to transversity basis
  double scale  = par[10];

  double result = scale * SinpsiOdd(x,par) * SinphiOdd(x,par);
  
  return result;
  
}

double CosthetaTrTotal(double *x, double *par)
{
  
  double result = CosthetaTrEven(x,par) + CosthetaTrOdd(x,par);
  
  return result;
  
}

double SinphiEven(double *x, double *par)
{
  
  //jpsiphi pdf three angles: this is the phi component
  //time + two other angles are integrated
  //double sqrtwo      = 1.41421356237; // Sqrt[2]
  double factor        = 0.04476232774; // 9 / (64 pi)
  double factor2       = 0.02083333333; // (1/48)
  double factor3       = 170.666666666; // 512/3.0;
  double factor4       =  85.333333333; // 256/3.0;
    
  double sinphi      = x[0];
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  double wphase   = par[6];
  double xs      = par[7];
  double omega    = par[8];
  double tmin     = par[9];
  //double scale    = par[10];
  
  double dms = xs / Gamma;

  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  double KpKp=(1-Rt)*Rp;
  double K0K0=(1-Rt)*(1-Rp);
  
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
  
  
  //Bar (Bs)
  double B0B0 = 0.5*K0K0*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			  +(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
								      +dms*cos(dms*tmin))));
  
  double BpBp = 0.5*KpKp*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			  +(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
								      +dms*cos(dms*tmin))));
						  
  ///////////////////////////////////////// 
  /* factor               = 0.04476232774; // 9 / (64 pi)
     double sqrtwo        = 1.41421356237; // Sqrt[2]
     double factor        = 0.04476232774; // 9 / (64 pi)
     double factor2       = 0.02083333333; // (1/48)
     double factor3       = 170.666666666; // 512/3.0;
     double factor4       =  85.333333333; // 256/3.0; 
  */

  /////////////////////////////////////////
  //W+
  double v1 =  factor2 * (  factor3 * A0A0 
			      + factor3 * ApAp 
			      - factor4 * ApAp * (1 - 2.0*sinphi*sinphi) );
  
  //W-
  double v2 =  factor2 * (  factor3 * B0B0 
			      + factor3 * BpBp 
			      - factor4 * BpBp * (1 - 2.0*sinphi*sinphi) );
  
  
  return factor * ((omega*v1)+(1.0-omega)*v2);
  
}

double SinphiOdd(double *x, double *par)
{
  
  
  //jpsiphi pdf three angles: this is the phi component
  //time + two other angles are integrated

  double factor        = 0.04476232774; // 9 / (64 pi)
  double factor2       = 0.02083333333; // (1/48)
  double factor3       = 170.666666666; // 512/3.0;
  double factor4       =  85.333333333; // 256/3.0;
  
  double sinphi      = x[0];
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  //double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double tmin     = par[9];
  //double scale    = par[10];
  
  double dms = xs / Gamma;

  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  double KtKt=Rt;
  
  //Time depency integrated over
  double gammadms = (1.0 / (Gamma*Gamma + dms*dms) );
  
  double AtAt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			  +(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
								      +dms*cos(dms*tmin))));
  //Bar (Bs)
  
  double BtBt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			  +(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
								      +dms*cos(dms*tmin))));
  
  ///////////////////////////////////////// 
  /* factor            = 0.04476232774; // 9 / (64 pi)
     double sqrtwo     = 1.41421356237; // Sqrt[2]
     double factor        = 0.04476232774; // 9 / (64 pi)
     double factor2       = 0.02083333333; // (1/48)
     double factor3       = 170.666666666; // 512/3.0;
     double factor4       =  85.333333333; // 256/3.0; 
  */

  /////////////////////////////////////////
  //W+
  double v1 =  factor2 * (  factor3 * AtAt 
			    + factor4 * AtAt * (1.0 - 2.0*sinphi*sinphi) );
  
  //W-
  double v2 =  factor2 * (  factor3 * BtBt 
			    + factor4 * BtBt * (1.0 - 2.0*sinphi*sinphi) );
  
  
  
  return factor * ((omega*v1)+(1.0-omega)*v2);
  
  
}


double phiPDFEven(double *x, double *par)
{
  
  
  //jpsiphi pdf three angles: this is the phi component
  //time + two other angles are integrated

  double factor        = 0.04476232774; // 9 / (64 pi)
  double factor2       = 0.02083333333; // (1/48)
  double factor3       = 170.666666666; // 512/3.0;
  double factor4       =  85.333333333; // 256/3.0;
    
  double phi      = x[0];
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double tmin     = par[9];
  double scale    = par[10];
  
  double dms = xs / Gamma;
  
  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  double KpKp=(1-Rt)*Rp;
  double K0K0=(1-Rt)*(1-Rp);
  
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

  //Bar (Bs)
  double B0B0 = 0.5*K0K0*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			  +(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
								      +dms*cos(dms*tmin))));
  
  double BpBp = 0.5*KpKp*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			  +(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
								      +dms*cos(dms*tmin))));
						  
  ///////////////////////////////////////// 
  /* factor               = 0.04476232774; // 9 / (64 pi)
     double sqrtwo        = 1.41421356237; // Sqrt[2]
     double factor        = 0.04476232774; // 9 / (64 pi)
     double factor2       = 0.02083333333; // (1/48)
     double factor3       = 170.666666666; // 512/3.0;
     double factor4       =  85.333333333; // 256/3.0; 
  */
  
  /////////////////////////////////////////
  //W+
  double v1 =  factor2 * (  factor3 * A0A0 
			      + factor3 * ApAp 
			      - factor4 * ApAp * cos(2*phi) );
  
  //W-
  double v2 =  factor2 * (  factor3 * B0B0 
			      + factor3 * BpBp 
			      - factor4 * BpBp * cos(2*phi));
  
  
  
  return scale * factor * ((omega*v1)+(1.0-omega)*v2);
  
  
}


double phiPDFOdd(double *x, double *par)
{
  
  
  //jpsiphi pdf three angles: this is the phi component
  //time + two other angles are integrated
  //double sqrtwo        = 1.41421356237; // Sqrt[2]
  double factor        = 0.04476232774; // 9 / (64 pi)
  double factor2       = 0.02083333333; // (1/48)
  double factor3       = 170.666666666; // 512/3.0;
  double factor4       =  85.333333333; // 256/3.0;
    
  double phi      = x[0];
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  //double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  double wphase   = par[6];
  double xs       = par[7];
  double omega    = par[8];
  double tmin     = par[9];
  double scale    = par[10];
  
  double dms = xs / Gamma;
  
  double G_H=Gamma*(1-DGrate/2.);
  double G_L=Gamma*(1+DGrate/2.);
  
  double KtKt=Rt;
  
  //Time depency integrated over
  double gammadms = (1.0 / (Gamma*Gamma + dms*dms) );
  
 

  double AtAt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			  +(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
								      +dms*cos(dms*tmin))));
  
  //Bar (Bs)
  
  
  double BtBt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
			  +(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			  + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)
								      +dms*cos(dms*tmin))));
  ///////////////////////////////////////// 
  /* factor               = 0.04476232774; // 9 / (64 pi)
     double sqrtwo        = 1.41421356237; // Sqrt[2]
     double factor        = 0.04476232774; // 9 / (64 pi)
     double factor2       = 0.02083333333; // (1/48)
     double factor3       = 170.666666666; // 512/3.0;
     double factor4       =  85.333333333; // 256/3.0; 
  */
  
  /////////////////////////////////////////
  //W+
  double v1 =  factor2 * (  factor3 * AtAt 
			      + factor4 * AtAt * cos(2*phi) );
  
  //W-
  double v2 =  factor2 * (  factor3 * BtBt 		      
			      + factor4 * BtBt * cos(2*phi) );
  
  
  
  return scale * factor * ((omega*v1) + (1.0-omega)*v2);
  
  
}
