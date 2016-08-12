// Include files

#include "OneAngleFit.h"
//-----------------------------------------------------------------------------
// 2006-11-27 : Andres Osorio 
//---------------------------------------------------------------------------

///////////////////
//Globals

std::vector<double> xi[6];
int MAXEVTS=100;
int MAXPARS=15;

/////////////////////


void fcn(int &npar, double *grad, double &fval, double *par, int iflag) 
{
  
  double sum  = 0.0;
  double lnfi = 0.0;
  double xx[4];
  
  for(int i = 0; i < MAXEVTS; ++i) {
    
    xx[0]=xi[0][i]; //propertime
    xx[1]=xi[1][i];
    
    double num = TotalPDF ( xx, par ); 
    lnfi = TMath::Log (num);
    sum += lnfi;
  }
  
  fval = -1.0*sum;
  
}

int getData(const char *infile) {
  
  ifstream * in = new ifstream(infile);
  
  if(!in->is_open()) {
    print_message("could not open input file.");
    print_message(infile);
    exit(1);
  }
  
  print_message("opening files:");
  
  double theta(0.);
  double time (0.);
  int    qfact(0);
  int n_events(0);
  int i(0);
  
  while (1) {
    
    (*in) >> theta >> time >> qfact;
    
    xi[0].push_back(time);
    xi[1].push_back(theta);
    xi[2].push_back(0.0);
    xi[3].push_back(0.0);
    xi[4].push_back(0.0);
    xi[5].push_back(0.0);
    
    if ( i < 5 ) {
      std::cout << i << '\t' 
                << xi[0][i] << '\t' 
                << xi[1][i] << '\t' 
                << xi[2][i] << '\t' 
                << xi[3][i] << '\t' 
                << xi[4][i] << '\t' 
                << xi[5][i] << std::endl;
    }
    
    ++i;
    
    ++n_events;
    
    if (!in->good()) break;
    
    
  }    
  
  in->close();
  
  delete in;
  
  return n_events;
  
  std::cout << "getData> Total events read:" << std::endl;
  std::cout << xi[0].size() << '\t' 
            << xi[1].size() << '\t' 
            << xi[2].size() << '\t' 
            << xi[3].size() << '\t' 
            << xi[4].size() << '\t' 
            << xi[5].size() << std::endl;
  
  return n_events;
  
}

void fitData( const char * infile, const char * paramfile , int nevts)
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
  
  int maxentries = getData( infile ); //Open event file
  
  if ( maxentries < nevts ) {
    std::cout << "fitData> nevts higher than entries" << std::endl;
    MAXEVTS = maxentries;
  } else MAXEVTS = nevts;
  
  std::cout << "Will run the fit for: " << MAXEVTS << " events" << std::endl;
  
  readData( paramfile, v1, v2);
  
  if( (int)v1.size() < MAXPARS ) {std::cout << "fitData> Error readingData" << '\n'; exit(1);}
  
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
    
  gMinuit->mnparm(0, "w",         vstart[0], step[0],   0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(1, "q",         vstart[1], step[1],   0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(2, "Amix",      vstart[2], step[2],   0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(3, "DGamma",    vstart[3], step[3],   0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(4, "Gamma",     vstart[4], step[4],   0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(5, "dm",        vstart[5], step[5],   0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(6, "mean",      vstart[6], step[6],   0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(7, "width",     vstart[7], step[7],   0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(8, "sPT",       vstart[8], step[8],   0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(9, "Rt",        vstart[9], step[9],   0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(10,"f1",        vstart[10], step[10], 0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(11,"s1",        vstart[11], step[11], 0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(12,"s2",        vstart[12], step[12], 0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(13,"mu1",       vstart[13], step[13], 0.0    , 0.0    ,ierflg);
  gMinuit->mnparm(14,"mu2",       vstart[14], step[14], 0.0    , 0.0    ,ierflg);
   
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

  std::cout << "MIGRAD: status: " << status << std::endl;
  
  /* 
     Always run HESSE:
     "As a general rule, anyone seriously interested in the parameter errors 
     should always put at least a HESSE command after each MIGRAD 
     (or MINIMIZE) command." ( from Minuit Manual )
  */
  
  gMinuit->mnexcm("HESSE", arglist , 1, ierflg);
  
  //gMinuit->mnhess();

  std::cout << "HESSE: status: " << ierflg << std::endl;
  
  // Print results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(3,amin);
  
  ///Store output from fit
  char results[30];
  std::string file = std::string( infile );
  int pos1 = 7;
  int pos2 = file.find(".",pos1);
  std::string opt = file.substr(pos1,(pos2-pos1));
  
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
