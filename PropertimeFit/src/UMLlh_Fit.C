// $Id: UMLlh_Fit.C,v 1.22 2007/04/03 00:31:09 aosorio Exp $
// Include files

#include "UMLlh_Fit.h"
//-----------------------------------------------------------------------------
// Implementation file  : UMLlh_Fit.C
// Using the method of unbinned maximum likelihood
// 2006-07-23 : Andres Osorio Oliveros
//---------------------------------------------------------------------------

//...............................
// Globals

std::vector<double> xi[9];
int MAXEVTS=100;
int MAXPARS=18;
//...............................

int getData(const char *infile) {
  
  //Clear containers:
  for ( int i = 0; i < 9; ++i) xi[i].clear();
  
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
  
  Double_t var[10];
  
  T->SetBranchAddress("ptnres"   ,&var[0]);
  T->SetBranchAddress("theta"    ,&var[1]);
  T->SetBranchAddress("psi"      ,&var[2]);
  T->SetBranchAddress("phi"      ,&var[3]);
  T->SetBranchAddress("ptwres"   ,&var[4]);
  T->SetBranchAddress("thwres"   ,&var[5]);
  T->SetBranchAddress("mass"     ,&var[6]);
  T->SetBranchAddress("tag"      ,&var[7]);
  T->SetBranchAddress("ptsigma"  ,&var[8]);
  
  int n_events = T->GetEntriesFast();
  
  for ( int i = 0; i< n_events; ++i) {
    
    T->GetEntry(i);
    xi[0].push_back(var[0]);
    xi[1].push_back(var[1]);
    xi[2].push_back(var[2]);
    xi[3].push_back(var[3]);
    xi[4].push_back(var[4]);
    xi[5].push_back(var[5]);
    xi[6].push_back(var[6]);
    xi[7].push_back(var[7]); // <- tag
    xi[8].push_back(var[8]); // <- tag
    
  }
  
  std::cout << "getData> Total events read:" << std::endl;
  std::cout << xi[0].size() << '\t' 
            << xi[1].size() << '\t' 
            << xi[2].size() << '\t' 
            << xi[3].size() << '\t' 
            << xi[4].size() << '\t'
            << xi[5].size() << '\t'
            << xi[6].size() << '\t'
            << xi[7].size() << '\t'
            << xi[8].size() << std::endl;
  
  return n_events;
    
}

void stepOne(const char *infile, int nevts, int fitopt)
{
  
  //using TMinuit to find the Maximum of Loglikelihood
  //gROOT->Reset();
  
  double par[MAXPARS];
  double errpar[MAXPARS];
  double fitpar[MAXPARS];
  double fixpar[MAXPARS];
  std::vector<double> v1;
  std::vector<double> v2;
  
  for(int i = 0; i< MAXPARS; ++i) {
    par[i]    = 0.0;
    errpar[i] = 0.0;
    fitpar[i] = 0.0;
    fixpar[i] = 0.0;
  }
  
  int maxentries = getData( infile ); //Open event file
  
  if ( maxentries < nevts ) {
    std::cout << "stepOne> nevts higher than entries" << std::endl;
    MAXEVTS = maxentries;
  } else MAXEVTS = nevts;
  
  std::cout << "Will run the fit for: " << MAXEVTS << " events" << std::endl;
  
  //read initial parameter values
  char initpars[30];
  std::string file = std::string( infile );
  int pos1 = 11;
  int pos2 = file.find(".",pos1);
  std::string opt = file.substr(pos1,(pos2-pos1));
  
  sprintf(initpars,"initialSO_%s.dat",opt.c_str());
  
  readData(initpars, v1, v2);
  
  if( (int)v1.size() < MAXPARS ) {std::cout << "stepOne> Error readingData" << '\n'; exit(1);}
  
  for(int i =0; i < MAXPARS; ++i) {
    fitpar[i] = v1[i];
    fixpar[i] = (int) v2[i];
  }
  
  ///Start a new Fit
  TMinuit *gMinuit = new TMinuit( MAXPARS ); 
  
  std::cout << "stepOne> \n"
            << "Fit option set to " << fitopt << ": using resolution model. \n";
  gMinuit->SetFCN(fcnwresol);
  
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
  
  gMinuit->mnparm(0, "Gamma",     vstart[0], step[0], 0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(1, "DGamma",    vstart[1], step[1], 0.0    , 0.00    ,ierflg);
  gMinuit->mnparm(2, "Rt",        vstart[2], step[2], 0.0    , 0.00    ,ierflg);
  gMinuit->mnparm(3, "R0",        vstart[3], step[3], 0.0    , 0.00    ,ierflg);
  gMinuit->mnparm(4, "tphase1",   vstart[4], step[4], 0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(5, "tphase2",   vstart[5], step[5], 0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(6, "wphase",    vstart[6], step[6], 0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(7, "dms",       vstart[7], step[7], 16.0   , 22.0    ,ierflg);
  gMinuit->mnparm(8, "omega",     vstart[8], step[8],  0.0   , 0.0    ,ierflg);
  gMinuit->mnparm(9, "sPT",       vstart[9], step[9],  0.0   , 0.0    ,ierflg);
  gMinuit->mnparm(10,"f1",        vstart[10], step[10], 0.0   , 0.0   ,ierflg);
  gMinuit->mnparm(11,"mu1",       vstart[11], step[11], 0.0   , 0.0   ,ierflg);
  gMinuit->mnparm(12,"s1",        vstart[12], step[12], 0.0   , 0.0   ,ierflg);
  gMinuit->mnparm(13,"mu2",       vstart[13], step[13], 0.0   , 0.0   ,ierflg);
  gMinuit->mnparm(14,"s2",        vstart[14], step[14], 0.0   , 0.0   ,ierflg);
  gMinuit->mnparm(15,"tag",       vstart[15], step[15], 0.0   , 0.0   ,ierflg);
  gMinuit->mnparm(16,"fsig",      vstart[16], step[16], 0.0   , 0.0   ,ierflg);
  gMinuit->mnparm(17,"tbkg",      vstart[17], step[17], 0.0   , 0.0   ,ierflg);
  
  
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
  
  gMinuit->mnhess();
  
  // Print results
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(3,amin);
  
  if ( status == 0 ) status = 0;
  else status = 4;
  
  
    
  ///Store output from fit
  char results[30];
  sprintf(results,"resultsStepOne_%s.dat",opt.c_str());
  
  std::ofstream *os = new std::ofstream(results,ofstream::out);
  
  //increase the output precision
  (*os) << setprecision(14) << "//Results from the log-likelihood fit" << std::endl;
  
  //For backwards compatibility add 10 units to status
  
  (*os) << (status+10) << '\t' << MAXPARS << std::endl;
  
  if ( status == 0 ) {
    
    for(int i=0;i<MAXPARS;++i) {
      fitpar[i] =0.0;
      gMinuit->GetParameter(i, fitpar[i], errpar[i]);
      (*os) << fitpar[i] << '\t' << errpar[i] << '\n';
    }
    
    (*os) << "Amin:" << '\t' << amin << '\t' << edm  << '\t'
          << errdef << '\n';
  }
  
  (*os) << icstat << '\t' << nfree << '\n'; // <- output the status of errmat calculation
  
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

void fcnwresol(int &npar, double *grad, double &fval, double *par, int iflag) 
{
  
  
  /*
    fcn function considering resolution effects
    Normalisation factor is numerically obtained using MathMore routines 
    2007-02-08 : Andres Osorio Oliveros
  */
  
  int q             = 0;
  double xx[8];
  double normfact[3];
  
  double sumoflogs  = 0.0;
  double lnfi       = 0.0;
  double func       = 0.0;
  
  par[15]=-1.0;
  normfact[0] = nfactorjpsiphi( par );
    
  par[15]= 0.0;
  normfact[1] = nfactorjpsiphi( par );
  
  par[15]= 1.0;
  normfact[2] = nfactorjpsiphi( par );
  
  for(int i = 0; i < MAXEVTS; ++i) {
    
    xx[0]=xi[4][i]; //propertime with resolution
    xx[5]=xi[7][i]; //tag
    xx[7]=xi[8][i]; //evt-by-evt resolution
    
    q= (int)xx[5];
    xx[6]  = normfact[q+1];
    
    func   = totalpdfWRes( xx , par );

    if ( func <= 0.0 ) {
      std::cout << "fcnwres> totalpdfWRes returns 0! at evt: " << i << std::endl;
      continue;
    }
          
    lnfi   = TMath::Log(func);
    sumoflogs = lnfi + sumoflogs;
    
  }
  
  fval = -2.0 * ( sumoflogs );

}


