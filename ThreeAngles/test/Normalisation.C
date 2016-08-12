
double normfactor(double *x, double *par) 
{

  //double factor = 0.04476232774; // 9 / (64 pi)
  //double factor =  22.3402144255274173;
  double factor = 1.0;
  double sqrtwo = 1.41421356237; // Sqrt[2]
  
  double time      = x[0];
  double theta     = x[1];   // l+ polar angle in Jpsi rest frame
  double psi       = x[2];   // K+ polar angle in Phi rest frame
  double phi       = x[3];   // sum of azimuthal angles
  
  double Gamma    = par[0];
  double DGrate   = par[1];
  double Rt       = par[2];
  double Rp       = par[3];
  double tphase1  = par[4];
  double tphase2  = par[5];
  double wphase   = par[6];
  double dms      = par[7];
  double omega    = par[8];
  double tmin     = par[9];
  double scale    = par[10];
  

  Double_t G_H=Gamma*(1+DGrate/2.);
  Double_t G_L=Gamma*(1-DGrate/2.);
  
  Double_t KtKt=Rt;
  Double_t KpKp=(1-Rt)*Rp;
  Double_t K0K0=(1-Rt)*(1-Rp);
  
  Double_t K0Kp=sqrt(K0K0*KpKp);
  Double_t K0Kt=sqrt(K0K0*KtKt);
  Double_t KpKt=sqrt(KpKp*KtKt);
  
  //Bs
  //Time depency integrated over

  tmin = 0.0;

  Double_t gammadms = (1.0 / (Gamma*Gamma + dms*dms) );
  
  Double_t A0A0 = 0.5* K0K0*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			     + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  Double_t ApAp = 0.5* KpKp*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			     + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  Double_t AtAt = 0.5* KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			     - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  Double_t ReA0Ap = 0.5*K0Kp*cos(tphase2-tphase1)*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
						   +(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
						   +2*sin(wphase)*gammadms
						   *(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  Double_t ImApAt = KpKt*((sin(tphase1)*gammadms*exp(-Gamma*tmin)
			   *(Gamma*cos(dms*tmin) - dms*sin(dms*tmin)))
			  - cos(tphase1)*cos(wphase)*gammadms*exp(-Gamma*tmin)
			  *(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))
			  -0.5*((1.0/G_H)*exp(-G_H*tmin)-(1.0/G_L)*exp(-G_L*tmin))*cos(tphase1)*sin(wphase));
  
  Double_t ImA0At = K0Kt*((sin(tphase2)*gammadms*exp(-Gamma*tmin)
			   *(Gamma*cos(dms*tmin) - dms*sin(dms*tmin)))
			  - cos(tphase2)*cos(wphase)*gammadms*exp(-Gamma*tmin)
			  *(Gamma*sin(dms*tmin) + dms*cos(dms*tmin))
			  -0.5*((1.0/G_H)*exp(-G_H*tmin)-(1.0/G_L)*exp(-G_L*tmin))*cos(tphase2)*sin(wphase));
  
  
  //Bar (Bs)
  Double_t B0B0 = 0.5*K0K0*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  Double_t BpBp = 0.5*KpKp*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    - 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  Double_t BtBt = 0.5*KtKt*((1-cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)+(1+cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
			    + 2*gammadms*sin(wphase)*(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  
  
  Double_t ReB0Bp = 0.5*K0Kp*cos(tphase2-tphase1)*((1+cos(wphase))*(1.0/G_L)*exp(-G_L*tmin)
						   +(1-cos(wphase))*(1.0/G_H)*exp(-G_H*tmin)
						   -2*sin(wphase)*gammadms
						   *(exp(-Gamma*tmin)*(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))));
  
  Double_t ImBpBt = -KpKt*((sin(tphase1)*gammadms*exp(-Gamma*tmin)
			    *(Gamma*cos(dms*tmin) - dms*sin(dms*tmin)))
			   - cos(tphase1)*cos(wphase)*gammadms*exp(-Gamma*tmin)
			   *(Gamma*sin(dms*tmin)+dms*cos(dms*tmin))
			   +0.5*((1.0/G_H)*exp(-G_H*tmin)-(1.0/G_L)*exp(-G_L*tmin))*cos(tphase1)*sin(wphase));
  
  Double_t ImB0Bt = -K0Kt*((sin(tphase2)*gammadms*exp(-Gamma*tmin)
			    *(Gamma*cos(dms*tmin) - dms*sin(dms*tmin)))
			   - cos(tphase2)*cos(wphase)*gammadms*exp(-Gamma*tmin)
			   *(Gamma*sin(dms*tmin) + dms*cos(dms*tmin))
			   +0.5*((1.0/G_H)*exp(-G_H*tmin)-(1.0/G_L)*exp(-G_L*tmin))*cos(tphase2)*sin(wphase));
  
  
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
  
  return factor* (omega*v1+(1.0-omega)*v2);

}
