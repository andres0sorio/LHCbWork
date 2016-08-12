// $Id: RandomDists.C,v 1.3 2006/11/24 22:14:57 aosorio Exp $
// Include files 

// local
#include "RandomDists.h"

//
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/RootFinder.h"
#include "Math/RootFinderAlgorithms.h"
#include "Math/WrappedFunction.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

//-----------------------------------------------------------------------------
// Implementation file for class : RandomDists
//
// 2006-07-04 : Andres Osorio Oliveros
//-----------------------------------------------------------------------------

double myfunc(double x, double *par);

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
RandomDists::RandomDists(  ) {
  
  dist[0]  = new TH1D("h1","This is the total distribution",100,0.0,20.0);
  dist[0]->Sumw2();
  dist[1]  = new TH1D("h2","This is the total distribution",100,0.0,3.15);
  dist[1]->Sumw2();
  dist[2]  = new TH1D("h3","This is the total distribution",100,0.0,3.15);
  dist[2]->Sumw2();
  dist[3]  = new TH1D("h4","This is the total distribution",100,0.0,6.5);
  dist[3]->Sumw2();
  dist[4]  = new TH1D("h5","This is the total distribution",100,-5.0,20.0);
  dist[4]->Sumw2();
  dist[5]  = new TH1D("h6","This is the total distribution",100,-5.0,20.0);
  dist[5]->Sumw2();
  dist[6]  = new TH1D("h7","This is the total distribution",100,-5.0,20.0);
  dist[6]->Sumw2();
  dist[7]  = new TH1D("h8","This is the total distribution",100,-5.0,20.0);
  dist[7]->Sumw2();
  
  max_events = 1000;
  
  std::cout << "RandomDists> the dafault MaxEvents is: " << max_events << std::endl;
  
  t1 = new TTree("data","MC toy events");
  
  outfile = new TFile("mctoydata.root","RECREATE");

  var[0] = 0.0;
  var[1] = 0.0;
  var[2] = 0.0;
  var[3] = 0.0;
  var[4] = 0.0;
  var[5] = 0.0;
  var[6] = 0.0;
  var[7] = 0.0;
  var[8] = 0.0;
  
  outfile = new TFile("mctoydata.root","RECREATE");
  
  t1 = new TTree("data","MC toy events");
  
  t1->Branch("ptnres"    ,&var[0] ,"ptnres/D");
  t1->Branch("theta"     ,&var[1] ,"theta/D");
  t1->Branch("psi"       ,&var[2] ,"psi/D");
  t1->Branch("phi"       ,&var[3] ,"phi/D");
  t1->Branch("ptwres"    ,&var[4] ,"ptwres/D");
  t1->Branch("thwres"    ,&var[5] ,"thwres/D");
  t1->Branch("ptaccnres" ,&var[6] ,"ptaccnres/D");
  t1->Branch("ptaccwres" ,&var[7] ,"ptaccwres/D");
  t1->Branch("mass"      ,&var[8] ,"mass/D");

  pi = TMath::Pi();

  params = new pdf_params();
    

}

RandomDists::RandomDists( const char *filename , int max ) {
  
  dist[0]  = new TH1D("h1","This is the total distribution",100,0.0,20.0);
  dist[0]->Sumw2();
  dist[1]  = new TH1D("h2","This is the total distribution",100,0.0,3.15);
  dist[1]->Sumw2();
  dist[2]  = new TH1D("h3","This is the total distribution",100,0.0,3.15);
  dist[2]->Sumw2();
  dist[3]  = new TH1D("h4","This is the total distribution",100,0.0,6.5);
  dist[3]->Sumw2();
  dist[4]  = new TH1D("h5","This is the total distribution",100,0.0,20.0);
  dist[4]->Sumw2();
  dist[5]  = new TH1D("h6","This is the total distribution",100,-5.0,20.0);
  dist[5]->Sumw2();
  dist[6]  = new TH1D("h7","This is the total distribution",100,-5.0,20.0);
  dist[6]->Sumw2();
  dist[7]  = new TH1D("h8","This is the total distribution",100,-5.0,20.0);
  dist[7]->Sumw2();

  max_events = max;
  
  var[0] = 0.0;
  var[1] = 0.0;
  var[2] = 0.0;
  var[3] = 0.0;
  var[4] = 0.0;
  var[5] = 0.0;
  var[6] = 0.0;
  var[7] = 0.0;
  var[8] = 0.0;  

  outfile = new TFile(filename , "RECREATE");
  
  t1 = new TTree("data","MC toy events");
  
  t1->Branch("ptnres"    ,&var[0] ,"ptnres/D");
  t1->Branch("theta"     ,&var[1] ,"theta/D");
  t1->Branch("psi"       ,&var[2] ,"psi/D");
  t1->Branch("phi"       ,&var[3] ,"phi/D");
  t1->Branch("ptwres"    ,&var[4] ,"ptwres/D");
  t1->Branch("thwres"    ,&var[5] ,"thwres/D");
  t1->Branch("ptaccnres" ,&var[6] ,"ptaccnres/D");
  t1->Branch("ptaccwres" ,&var[7] ,"ptaccwres/D");
  t1->Branch("mass"      ,&var[8] ,"mass/D");
  
  pi = TMath::Pi();
  
  params = new pdf_params();

}


//=============================================================================
// Destructor
//=============================================================================
RandomDists::~RandomDists() 
{

  for(int i=0; i <8; ++i){
    if(dist[i]) {
      dist[i]->Delete();
    }
  }
  
  if(outfile) {
    outfile->Close();
    delete outfile;
  }
  
  delete params;
    
} 

//=============================================================================
void RandomDists::RunExperiment(double *pm, int npar) 
{
  
  //Invoking MathMore functions
  //double test = ROOT::Math::cyl_bessel_i(2.0,2.0);
  
  double xx[1];
  double x = 0.0;
  double u = 0.0;
  double f = 0.0;
  double fmax[4];

  for( int i = 0; i < npar; ++i) {
    par[i] = pm[i];
    std::cout << par[i] << std::endl;
  }

  findMaxima(par);
  
  fmax[0] =ptMax;
  fmax[1] =thetaMax;
  fmax[2] =psiMax;
  fmax[3] =phiMax;
  
  Ptr2pdf astPDF[4];
  astPDF[0] = &properTimePDF;
  astPDF[1] = &thetaPDF;
  astPDF[2] = &psiPDF;
  astPDF[3] = &phiPDF;
  
  //Generate random numbers
  TRandom3 *rdn = new TRandom3();

  rdn->SetSeed();
  
  xx[0]=0.0;
  
  xmin[0]=0.0;  //ps
  xmax[0]=100.00;
  xmin[1]=0.0; //(theta)
  xmax[1]= pi;
  xmin[2]=0.0; //(psi)
  xmax[2]= pi;
  xmin[3]=0.0; //phi
  xmax[3]=2.0*pi;
  
  int counter = 0;
  Double_t r1(0.0);
  Double_t r2(0.0);
  
  for( int i=0; i < 4; ++i) {

    counter = 0;
    
    while (counter < max_events) {
      
      //r1 = gRandom->Uniform(1);
      r1 = rdn->Rndm();
      x = xmin[i] + r1*(xmax[i]-xmin[i]);
      xx[0] = x;
      
      //r2 = gRandom->Uniform(1);
      r2 = rdn->Rndm();

      u = r2 * fmax[i];
      f = (*astPDF[i])(xx,par);
      
      if ( u < f ) {
        dist[i]->Fill(xx[0]);	
        data[i].push_back(xx[0]);
        if( ((counter+1)%10000)==0 ) std::cout << "RunExperiment> " << (counter+1) << std::endl;
        ++counter;
      }
      
    }
    
  }
  
  if ( data[0].size() != 0 ) { 
    
    std::cout << "RunExperiment> event=0 " 
	      << data[0][0] << '\t'
	      << data[1][0] << '\t'
	      << data[2][0] << '\t'
	      << data[3][0] << std::endl;
    
  }
  
  //Detector resolution effects (use a convolution to model the effects)
  
  double convpars[6];
  convpars[0] = 0.030;
  convpars[1] = 0.80;     //f1
  convpars[2] = 0.0022;   //mu1
  convpars[3] = 0.0584;   //s1
  convpars[4] = -0.1114;  //mu2
  convpars[5] = 0.471;    //s2
  
  
  double ptmin(-2.0);
  double ptmax(10.0);
  
  counter = 0;
  //Propertime convolution
  while ( counter < max_events) {
    
    // r1 = gRandom->Uniform(1);
    r1 = rdn->Rndm();
    x = ptmin + r1*(ptmax-ptmin);
    xx[0] = x;
    
    // r2 = gRandom->Uniform(1);
    r2 = rdn->Rndm();
    
    u = r2 * 1.0;
    //f = convolve( properTimePDF, withGaussian, par, xx , convpars);
    f = convolve( properTimePDF, with2Gaussians, par, xx , convpars);
    
    if ( u < f ) {
      dist[4]->Fill(xx[0]);	
      data[4].push_back(xx[0]);
      if( ((counter+1)%10000)==0 ) std::cout << "RunExperiment> " << (counter+1) << std::endl;
      ++counter;
      
    }
    
  }
  
  double smearpars[1];
  smearpars[0] = 0.020;
  
  //will be theta smearing
  std::vector<double>::iterator itr;
  for (itr = data[1].begin(); itr != data[1].end(); ++itr) {
    double thetasmear = smearValue( withGaussian, (*itr), smearpars);
    data[5].push_back(thetasmear);
  }

  
  //Only Acceptance
  counter = 0;
  ptmin   = 0.0;
  ptmax   = 10.0;
  
  while ( counter < max_events) {
    
    r1 = rdn->Rndm();
    x = ptmin + r1*(ptmax-ptmin);
    xx[0] = x;
    
    r2 = rdn->Rndm();
    
    u = r2 * 1.0;
    
    //f = lfTimePDFAcc( xx , par );
    
    if ( u < f ) {
      dist[6]->Fill(xx[0]);	
      data[6].push_back(xx[0]);
      if( ((counter+1)%10000)==0 ) std::cout << "RunExperiment> " << (counter+1) << std::endl;
      ++counter;
      
    }
    
  }
  
  
  //Fill tree now
  for (int k = 0; k < max_events; ++k) {
    
    var[0] = data[0][k];
    var[1] = data[1][k];
    var[2] = data[2][k];
    var[3] = data[3][k];
    var[4] = data[4][k];
    var[5] = data[5][k];
    var[6] = data[6][k];
    var[7] = 0.0;
    var[8] = 0.0;
    t1->Fill();
    
  }
  
  t1->Write();
  
}

void RandomDists::findMaxima( double *par )
{
  
  //get the maximum for each pdf

  //For properTime:assume maximum is when it intersects the y-axis
  double xmax = 0.0;
  ptMax = properTimePDF(&xmax,par);
  
  //Theta:
  xmax = TMath::Pi() / 2.0;
  thetaMax    = thetaPDF(&xmax,par);
  
  setParameters( par );
  
  //psi:
  ///////////////////////////////
  //Using gsl
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  T = gsl_root_fsolver_brent;
  gsl_function F;
  F.function = &DpsiPDF;
  F.params = params;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);

  gsl_root_fsolver_set (s, &F, 0.01, 1.00);
  
  int iter = 0;
  int max_iter = 100;
  double x_lo  = 0.01;
  double x_hi  = 1.00;
  double froot = 0.9;
  
  int status = GSL_CONTINUE;
  
  while (status == GSL_CONTINUE && iter < max_iter)
  {
    status = gsl_root_fsolver_iterate (s);
    froot = gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
  
    status = gsl_root_test_interval (x_lo, x_hi,0, 0.0001);
    
    if (status == GSL_SUCCESS) 
      std::cout << " GSL: succeded after: " 
                << iter << " iterations" << std::endl;
    ++iter;
    
  }
  
  psiMax = psiPDF( &froot , par );
  
  ////////////////////////////////
  F.function = &DphiPDF;
  
  iter = 0;
  x_lo  = 1.50;
  x_hi  = 3.00;
  froot = 1.5;
    
  std::vector<double> minmax;
  
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  status = GSL_CONTINUE;
  
  while (status == GSL_CONTINUE && iter < max_iter)
  {
    status = gsl_root_fsolver_iterate (s);
    froot = gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    
    status = gsl_root_test_interval (x_lo, x_hi,0, 0.0001);
    
    if (status == GSL_SUCCESS) {
      std::cout << " GSL: succeded after: " 
                << iter << " iterations" << std::endl;
      minmax.push_back(phiPDF( &froot , par ));
    }
    
    ++iter;
  }
  
  
  phiMax = phiPDF( &froot , par );
  
  std::cout << "ptmax: "    << ptMax << '\t' 
            << "thetaMax: " << thetaMax << '\t' 
            << "psiMax: "   << psiMax << '\t'
            << "phiMax: "   << phiMax << std::endl;
  
  
}


void RandomDists::CopyDataTo(int v, double *out)
{
  
  std::vector<double>::iterator itr;
  int index = 0;
  for(itr = data[v].begin(); itr != data[v].end(); ++itr){
    out[index]=(*itr);
    ++index;
  }
  
}

void RandomDists::setParameters( double * par )
{
  
  
  params->gBar   = par[0];
  params->Dgamma = par[1];
  params->Rt     = par[2];
  params->R0     = par[3];
  params->dt1    = 0.0;
  params->dt2    = 0.0;
  params->phis   = par[6];
  params->Dms    = par[7];
  params->omega  = par[8];
  
  
}

