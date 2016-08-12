// $Id: GetLogPlots.C,v 1.1 2007/04/09 11:49:08 aosorio Exp $
{

  double par[25];
  
  par[0]=0.68;
  par[1]=0.0;
  par[2]=0.0;
  par[3]=0.0;
  par[4]=0.0;
  par[5]=3.1415;
  par[6]=0.0;
  par[7]=20.00;
  par[8]=0.50;
  par[9]=0.04;
  par[10]=0.80;
  par[11]=0.0022;
  par[12]=0.0584;
  par[13]=-0.111;
  par[14]=0.471;
  par[15]=0.0;
  par[16]=0.50;
  par[17]=0.555;

  double xx[8];
  xx[0] = -0.431357;
  xx[1] = 0.0;
  xx[2] = 0.0;
  xx[3] = 0.0;
  xx[4] = 0.0;
  xx[5] = 0.0;
  xx[6] = nfactorjpsiphi( par );
  xx[7] = 0.0448;
  
  double val = totalpdfWRes( xx , par );
  
  std::cout << val << std::endl;
  
  
}
