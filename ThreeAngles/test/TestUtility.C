// $Id: TestUtility.C,v 1.3 2006/11/24 22:14:57 aosorio Exp $
// Include files

#include "TestUtility.h"
//-----------------------------------------------------------------------------
// Implementation file  : TestUtility.C
// The method of unbinned maximum likelihood
// Testing CP measurements - multiangular fit
// 2006-08-24 : Andres Osorio 
//---------------------------------------------------------------------------


void fcn(int &npar, double *grad, double &fval, double *par, int iflag) 
{
  
  double sum  = 0.0;
  double lnfi = 0.0;
  
  double xx[4];

  for(int i = 0; i < MAXEVTS; ++i) {
    
    xx[0]=xi[0][i]; //propertime
    xx[1]=xi[1][i];
    xx[2]=xi[2][i];
    xx[3]=xi[3][i];
  
    //double num = angpdfwres(xx,par); //with convolution
    double num = jpsiphipdf ( xx, par ); //no convolution
    lnfi = TMath::Log (num);
    
    sum += lnfi;
    
  }
  
  fval = -1.0*sum;
  
}


double llh( double * x, double * par) {
  
  int npar = MAXPARS;
  int iflag = 1;
  double grad(1.0);
  double fval(0.);  
  
  int rp  = (int)par[15];
  par[rp] = x[0];

  fcn(npar, &grad, fval , par, iflag);
  
  return fval;
  
}


int TestOne(const char *infile, int nevts)
{
  
  gROOT->Reset();
  
  for(int k = 0; k < 6; ++ k) xi[k].clear();
  
  double fitpar[20];
  std::vector<double> v1;
  
  MAXEVTS = nevts;
  
  for(int i = 0; i< 20; ++i) {
    fitpar[i] = 0.0;
  }
  
  getData( infile ); //Open event file
  
  std::cout << "Will run the scan for: " << MAXEVTS << " events" << std::endl;
  
  //read initial parameter values
  char initpars[30];
  std::string file = std::string( infile );
  int pos1 = 11;
  int pos2 = file.find(".",pos1);
  std::string opt = file.substr(pos1,(pos2-pos1));
  
  sprintf(initpars,"inputParameters_%s.dat",opt.c_str());
  
  int SC = readData(initpars, v1);
  
  if ( !SC ) return 0;
  
  if( (int)v1.size() < MAXPARS ) {std::cout << "stepTwo> Error reading Data" << '\n'; return 0;}
  
  for(int i =0; i < MAXPARS; ++i) {
    fitpar[i] = v1[i];
    std::cout <<  i << '\t' 
	      << fitpar[i] << '\n';
  }
  
  int  np = 10;
  double x[1];
  double xp[100];
  double yp[100];
  double fval(0.);
  
  double xmin[15];
  double xmax[15];
  double step[15];

  //xmin xmax definitions
  xmin[0] = 0.50; //Gamma Bar
  xmax[0] = 0.90;
  xmin[1] = -0.10; //DeltaGamma
  xmax[1] = 0.40;
  xmin[2] = 0.10; //Rt
  xmax[2] = 0.39;
  xmin[3] = 0.40; //R0
  xmax[3] = 0.80;
  xmin[4] = 0.00; //delta1
  xmax[4] = 0.30;
  xmin[5] = 3.10; //delta2
  xmax[5] = 3.20;
  xmin[6] = -0.15;//phi_s
  xmax[6] = 0.05;
  xmin[7] = 15.0; //DMs
  xmax[7] = 30.0;

  for(int i=0; i < 8; ++i) step[i] = (xmax[i] - xmin[i])/ (double)np;

  TCanvas *cv[1];
  cv[0]= new TCanvas("cv1","Plots",200,10,1000,700);
  cv[0]->Divide(4,2);
  
  TGraph *gt[15];

  int itest (0);

  itest = 2;

  for( itest = 0; itest < 8; ++itest) {

    fitpar[15] = itest;
    
    for(int i = 0; i < np; ++i) {
      x[0] = xmin[itest] + i*step[itest];
      fval = llh(x,fitpar);
      xp[i] = x[0];
      yp[i] = fval;
    }  
    
    //reset the fit parameters
    for(int k=0; k < MAXPARS; ++k) fitpar[k] = v1[k];
    
    if ( (itest != 5) && (itest != 7) ) {
      gt[itest] = new TGraph(np,xp,yp);
      cv[0]->cd(itest+1);
      gt[itest]->Draw("AC");
      cv[0]->Update();
    }
    else {std::cout << "test: " << itest << " ignored" << std::endl;}

    }
  
  
  return 1;

}

void getData(const char *infile) {
  
  TFile *f1 = new  TFile (infile, "READ");
  
  if( f1->IsOpen() ) {
    std::cout << "getData> File: " << std::string(infile) 
	      << " opened." << std::endl;
  }
  else {
    std::cout << "getData> Error: File: " << std::string(infile) 
	      << " not accessible" << std::endl;
    exit(1);
  }
  
  f1->cd();

  TTree *T = new TTree();

  T = (TTree*)gROOT->FindObject("data");
  
  Double_t var[7];
  
  T->SetBranchAddress("ptnres"   ,&var[0]);
  T->SetBranchAddress("theta"    ,&var[1]);
  T->SetBranchAddress("psi"      ,&var[2]);
  T->SetBranchAddress("phi"      ,&var[3]);
  T->SetBranchAddress("ptwres"   ,&var[4]);
  T->SetBranchAddress("thwres"   ,&var[5]);
  T->SetBranchAddress("mass"     ,&var[6]);
  
  
  int n_events = T->GetEntriesFast();
  
  for ( int i = 0; i< n_events; ++i) {
 
    T->GetEntry(i);
    xi[0].push_back(var[0]);
    xi[1].push_back(var[1]);
    xi[2].push_back(var[2]);
    xi[3].push_back(var[3]);
    xi[4].push_back(var[4]);
    xi[5].push_back(var[5]);

    if ( i < 5 ) {
      std::cout << i << '\t' 
		<< xi[0][i] << '\t' 
		<< xi[1][i] << '\t' 
		<< xi[2][i] << '\t' 
		<< xi[3][i] << '\t' 
		<< xi[4][i] << '\t' 
		<< xi[5][i] << std::endl;
    }
    
  }
  
  std::cout << "getData> Total events read:" << std::endl;
  std::cout << xi[0].size() << '\t' 
	    << xi[1].size() << '\t' 
	    << xi[2].size() << '\t' 
	    << xi[3].size() << '\t' 
	    << xi[4].size() << '\t' 
	    << xi[5].size() << std::endl;

  delete T;

  f1->Close();
  
  delete f1;


}

