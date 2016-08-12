// $Id: UMLlh_Fit.C,v 1.11 2006/11/29 08:46:35 aosorio Exp $
// Include files

#include "UMLlh_Fit.h"
//-----------------------------------------------------------------------------
// Implementation file  : UMLlh_Fit.C
// The method of unbinned maximum likelihood
// Testing CP measurements - multiangular fit
// 2006-07-23 : Andres Osorio Oliveros
//---------------------------------------------------------------------------

///////////////////
//Globals

std::vector<double> xi[6];
int MAXEVTS=100;
int MAXPARS=15;
/////////////////////

void getData(const char *infile) {
  
  TFile *f1 = new  TFile (infile, "READ");
  
  if( f1 ) {
    std::cout << "getData> File: " << std::string(infile) 
	      << " opened." << std::endl;
  }
  else {
    std::cout << "getData> Error: File: " << std::string(infile) 
	      << " not accessible" << std::endl;
    exit(1);
  }
  
  f1->cd();
  TTree *T = (TTree*)gROOT->FindObject("data");
  
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
    
  }
  
  std::cout << "getData> Total events read:" << std::endl;
  std::cout << xi[0].size() << '\t' 
	    << xi[1].size() << '\t' 
	    << xi[2].size() << '\t' 
	    << xi[3].size() << '\t' 
	    << xi[4].size() << '\t' 
	    << xi[5].size() << std::endl;
  
}

void stepOne(const char *infile, int nevts)
{
  
  //using TMinuit to find the Maximum of Loglikelihood
  //gROOT->Reset();
  
  double par[MAXPARS];
  double errpar[MAXPARS];
  double fitpar[MAXPARS];
  double fixpar[MAXPARS];
  std::vector<double> v1;
  std::vector<double> v2;
  
  MAXEVTS = nevts;
  
  for(int i = 0; i< MAXPARS; ++i) {
    par[i]    = 0.0;
    errpar[i] = 0.0;
    fitpar[i] = 0.0;
    fixpar[i] = 0.0;
  }
  
  getData( infile ); //Open event file
  
  std::cout << "Will run the fit for: " << MAXEVTS << " events" << std::endl;
  
  //read initial parameter values
  char initpars[30];
  std::string file = std::string( infile );
  int pos1 = 11;
  int pos2 = file.find(".",pos1);
  std::string opt = file.substr(pos1,(pos2-pos1));
  
  sprintf(initpars,"initialSO_%s.dat",opt.c_str());
  
  readData(initpars, v1, v2);
  
  if( (int)v1.size() < MAXPARS ) {std::cout << "stepTwo> Error readingData" << '\n'; exit(1);}
  
  for(int i =0; i < MAXPARS; ++i) {
    fitpar[i] = v1[i];
    fixpar[i] = (int) v2[i];
  }
  
  ///Start a new Fit
  
  TMinuit *gMinuit = new TMinuit( MAXPARS ); 
  gMinuit->SetFCN(fcn);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //Set starting values and step sizes for parameters
  Double_t vstart[MAXPARS];
  Double_t step[MAXPARS];
  
  for(int i=0; i < MAXPARS; ++i) {
    vstart[i] = fitpar[i];
    step[i]   = 0.001;
  }
  
  step[7]=0.1;

  gMinuit->mnparm(0, "Gamma",     vstart[0], step[0], 0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(1, "DGamma",    vstart[1], step[1], 0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(2, "Rt",        vstart[2], step[2], 0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(3, "R0",        vstart[3], step[3], 0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(4, "tphase1",   vstart[4], step[4], -0.1   , 0.1    ,ierflg);
  gMinuit->mnparm(5, "tphase2",   vstart[5], step[5],  3.0   , 3.5    ,ierflg);
  gMinuit->mnparm(6, "wphase",    vstart[6], step[6], -0.1   , 0.1    ,ierflg);
  gMinuit->mnparm(7, "dmas",      vstart[7], step[7], 15.0   , 25.0   ,ierflg);
  gMinuit->mnparm(8, "omega",     vstart[8], step[8],  0.0   , 0.0    ,ierflg);
  gMinuit->mnparm(9, "sPT",       vstart[9], step[9],  0.0   , 0.0    ,ierflg);
  gMinuit->mnparm(10,"f1",        vstart[10], step[10], 0.0   , 0.0    ,ierflg);
  gMinuit->mnparm(11,"mu1",       vstart[11], step[11], 0.0   , 0.0    ,ierflg);
  gMinuit->mnparm(12,"s1",        vstart[12], step[12], 0.0   , 0.0    ,ierflg);
  gMinuit->mnparm(13,"mu2",       vstart[13], step[13], 0.0   , 0.0    ,ierflg);
  gMinuit->mnparm(14,"s2",        vstart[14], step[14], 0.0   , 0.0    ,ierflg);
 
  //Fix parameters
  int nfixed = 0;
  int nfree  = 0;
  for ( int k = 0; k < MAXPARS; ++k) {
    if( fixpar[k] ) {
      gMinuit->FixParameter(k);
      ++nfixed;
    }
  }
  
  nfree = MAXPARS - nfixed;
  
  std::cout << "Number of free parameters: " << nfree << std::endl;
  

  // Now ready for minimization step
  arglist[0] = 800;
  arglist[1] = 1.;
  
  gMinuit->mnexcm("MIGRAD", arglist , 2, ierflg);
  int status = ierflg;
  
  /* 
     Always run HESSE:
     "As a general rule, anyone seriously interested in the parameter errors 
     should always put at least a HESSE command after each MIGRAD 
     (or MINIMIZE) command." ( from Minuit Manual )
  */
  
  //gMinuit->mnexcm("HESSE", arglist , 1, ierflg);
  gMinuit->mnhess();
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(3,amin);
  
  ///Store output from fit
  char results[30];
  sprintf(results,"resultsStepOne_%s.dat",opt.c_str());
  
  std::ofstream *os = new std::ofstream(results,ofstream::out);
  
  (*os) << "//Results from the log-likelihood fit" << std::endl;
  
  (*os) << status << '\t' << MAXPARS << std::endl;
  
  if ( status == 0 ) {
    
    for(int i=0;i<MAXPARS;++i) {
      fitpar[i] =0.0;
      gMinuit->GetParameter(i, fitpar[i], errpar[i]);
      (*os) << fitpar[i] << '\t' << errpar[i] << '\n';
    }
    
    (*os) << "Amin:" << '\t' << amin << '\n';
  }
  
  (*os) << status << '\t' << nfree << '\n';
  
  double errmat[100];
  int matdim = nfree;
  int k = 1;
  
  if ( nfree > 1 && status == 0 ) {
    gMinuit->mnemat(errmat, matdim);
    for(int i=0; i < matdim*matdim; ++i){
      (*os) << errmat[i] << " ";
      if ( k == matdim ) {
        (*os) << '\n';
        k = 0;
      }
      ++k;
    }
  }
  
  os->close();
  delete os;
  
}

void stepTwo(const char *infile, int nevts) {
  
  //using TMinuit explicitly to get the min of Loglikelihood
  //gROOT->Reset();
  
  double errpar[10];
  double fitpar[10];
  std::vector<double> v1;
  std::vector<double> v2;
  
  MAXEVTS = nevts;
  
  for(int i = 0; i< 10; ++i) {
    errpar[i] = 0.0;
    fitpar[i] = 0.0;
  }
  
  getData( infile );  
  std::cout << "Will run the fit for: " << MAXEVTS << " events" << std::endl;

  //read results from previous step

  char prevres[30];
  
  std::string file = std::string( infile );
  int pos1 = 11;
  int pos2 = file.find(".",pos1);
  std::string opt = file.substr(pos1,(pos2-pos1));

  sprintf(prevres,"resultsSO_%s.out",opt.c_str());
    
  readData(prevres, v1, v2);
  
  if( (int)v1.size() <8 ) {std::cout << "stepTwo> Error readingData" << '\n'; exit(1);}
  
  for(int i = 0; i< 10; ++i) {
    fitpar[i] = v1[i];
    errpar[i] = v2[i];
  }
  
  ///Start a new Fit
  
  TMinuit * gMinuit = new TMinuit( 8 ); //7 parameters describe the flight of the particle
  gMinuit->SetFCN(fcn);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //Set starting values and step sizes for parameters
  Double_t vstartp[8] = { fitpar[0],  
			  fitpar[1],
			  fitpar[2],
			  fitpar[3],
			  fitpar[4],
			  fitpar[5],
			  fitpar[6],
			  fitpar[7]};
  
  static Double_t step[8]   = {0.001, 0.001 , 0.001 , 0.001, 0.001, 0.01 , 0.001 , 0.1};
  
  gMinuit->mnparm(0, "Gamma",     vstartp[0], step[0], 0.0    , 2.0    ,ierflg);
  gMinuit->mnparm(1, "DGrate",    vstartp[1], step[1], 0.0    , 0.2    ,ierflg);
  gMinuit->mnparm(2, "Rt",        vstartp[2], step[2], 0.1    , 0.3    ,ierflg);
  gMinuit->mnparm(3, "Rp",        vstartp[3], step[3], 0.0    , 1.0    ,ierflg);
  gMinuit->mnparm(4, "tphase1",   vstartp[4], step[4], -0.1   , 0.1    ,ierflg);
  gMinuit->mnparm(5, "tphase2",   vstartp[5], step[5], 3.0    , 3.5    ,ierflg);
  gMinuit->mnparm(6, "wphase",    vstartp[6], step[6], -0.1   , 0.1    ,ierflg);
  gMinuit->mnparm(7, "dmas",      vstartp[7], step[7], 18.0   , 22.0   ,ierflg);
  
  //Fix parameters
  gMinuit->FixParameter(0);
  gMinuit->FixParameter(1);
  gMinuit->FixParameter(2);
  gMinuit->FixParameter(3);
  gMinuit->FixParameter(4);
  gMinuit->FixParameter(5);
  gMinuit->FixParameter(7);
  
  // Now ready for minimization step
  arglist[0] = 800;
  arglist[1] = 1.;
  
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(3,amin);
  
  std::ofstream *os = new std::ofstream("resultsStepTwo.out",ofstream::out);
  
  (*os) << "//Results from fit" << std::endl;
  
  for(int i=0;i<8;++i) {
    fitpar[i] =0.0;
    gMinuit->GetParameter(i, fitpar[i], errpar[i]);
    (*os) << fitpar[i] << '\t' << errpar[i] << '\n';
  }
  
  os->close();
  delete os;

  
}           

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
    //double num = properTimePDF ( xx , par );
    //double num = thetaPDF( xx , par );
    //double num = psiPDF ( xx , par );
    //double num = phiPDF ( xx , par );
    
    //double den = integrate( jpsiphipdf, xx, par, 0.0, 10.0 );
    //lnfi = TMath::Log (num) - TMath::Log (den);
    
    lnfi = TMath::Log (num);
    
    sum += lnfi;
    
  }
  
  fval = -1.0*sum;
  
}

