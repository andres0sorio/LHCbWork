// $Id: RandomDists.C,v 1.18 2006/12/09 18:17:56 aosorio Exp $
// Include files 

// local
#include "RandomDists.h"




//-----------------------------------------------------------------------------
// Implementation file for class : RandomDists
//
// 2006-07-04 : Andres Osorio Oliveros
//-----------------------------------------------------------------------------

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
  
  rdn[0] = new TRandom3(0);
  rdn[1] = new TRandom3(0);
  rdn[2] = new TRandom3(0);
  rdn[3] = new TRandom3(0);
  reg    = new TRandom3(0);

}

RandomDists::RandomDists( const char *filename , int max ) {

  max_events = max;
  
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
  
  params = new pdf_params();
  
  rdn[0] = new TRandom3(0);
  rdn[1] = new TRandom3(0);
  rdn[2] = new TRandom3(0);
  rdn[3] = new TRandom3(0);
  reg    = new TRandom3(0);
  tag    = new TRandom3(0);
  
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
    data[i].clear();
  }
  
  if(outfile) {
    outfile->Close();
    delete outfile;
  }
  
  delete params;
  for (int i=0; i < 4; ++i) delete rdn[i];
  delete reg;
  delete tag;

} 

void RandomDists::initialise(double *pm, int npar)
{
  
  pi = TMath::Pi();
  
  for( int i = 0; i < npar; ++i) {
    par[i] = pm[i];
    std::cout << par[i] << std::endl;
  }

  setParameters( pm );
  
  xmin[0]=0.0;  //ps
  xmax[0]=20.00;
  xmin[1]=0.0; //(theta)
  xmax[1]=pi;
  xmin[2]=0.0; //(psi)
  xmax[2]=pi;
  xmin[3]=-pi; //phi
  xmax[3]=pi;
  
  for ( int i=0; i < 8; ++i) {
    var[i] = 0.0;
  }
  
  for (int i=0; i< 4; ++i) {
    tempoval[i] = 0.0;
    foundval[i] = false;
  }
  
  pdfmax  = 0.0;
  
  double omega  = 0.30; // <- mistag fraction
  
  tagfrac[0][0] = (1.0-omega);
  tagfrac[0][1] = omega;
  tagfrac[1][0] = omega;
  tagfrac[1][1] = (1.0-omega);
  
  generatingBs = true;

  findMaxima();

  
  
}


//=============================================================================
void RandomDists::RunExperiment() 
{
  
  double u = 0.0;
  double f = 0.0;

  int counter = 0;
  double r1(0.0);
  double r2(0.0);
  
  counter = 0;
  
  while (counter < max_events) {
    
    SelectTagging(); // <- select whether is going to be treated as a Bs / Bs_bar

    updateMaxima();
    
    for ( int i = 0; i < 4; ++i ) {
      r1 =  rdn[i]->Rndm();
      tempoval[i] = xmin[i] + r1*(xmax[i]-xmin[i]);
    }
    
    r2 = reg->Rndm();
    u  = r2 * pdfmax;
    f  = (*astPDF)(tempoval , par); // <- evaluate the pdf
    
    if ( u < f ) {
      for ( int i = 0; i < 4; ++i ) foundval[i] = true;
    }
    
    if ( foundval[0] && foundval[1] && foundval[2] && foundval[3] ) {
      
      for( int k=0; k < 4; ++k) {
        dist[k]->Fill(tempoval[k]);
        if ( counter < 5 )  data[k].push_back(tempoval[k]);
        foundval[k] = false; // reset the found flag
      }
      
      WriteToTree();
      
      if( ((counter+1)%10000)==0 ) std::cout << "RunExperiment> " << (counter+1) << std::endl;
      ++counter;
    }
    
    
  }
  
  if ( data[0].size() != 0 ) { 
    
    std::cout << "RunExperiment> event=0 " 
              << data[0][0] << '\t'
              << data[1][0] << '\t'
              << data[2][0] << '\t'
              << data[3][0] << std::endl;
  }
  
  t1->Write();
  
}

void RandomDists::SelectTagging( )
{
  
  int    q     = 1;
  double bstag = tagfrac[1-q][0];
  double r3    = tag->Rndm();
  
  if ( r3 < bstag ) {
    //Bs
    generatingBs = true;
    astPDF = &jpsiphiWpPDF;
  } else {
    //Bs_bar
    generatingBs = false;
    astPDF = &jpsiphiWmPDF;
  }
  
}

void RandomDists::WriteToTree() 
{
  
  //Write to tree on a event by event basis - for large data sets
  var[0] = tempoval[0];
  var[1] = tempoval[1];
  var[2] = tempoval[2];
  var[3] = tempoval[3];
  var[4] = 0.0; // <- not yet implemented
  var[5] = 0.0;
  var[6] = 0.0;
  var[7] = 0.0;
  var[8] = 0.0;
  t1->Fill();
  
}


void RandomDists::findMaxima()
{
  
  //For properTime:assume maximum is when it intersects the y-axis
  double tmax  (0.);
  double xx[4];
  
  xx[0]      = properTimeWpPDF(&tmax,par);
  xx[1]      = findMaximum(1, thetaWpPDF, DthetaWpPDF);
  xx[2]      = findMaximum(2, psiWpPDF, DpsiWpPDF);
  xx[3]      = findMaximum(3, phiWpPDF, DphiWpPDF);
    
  fmax[0]    = jpsiphiWpPDF( xx , par );
  std::cout << "max{W+} = " << fmax[0] << std::endl;
  
  xx[0]      = properTimeWmPDF(&tmax,par);
  xx[1]      = findMaximum(1, thetaWmPDF, DthetaWmPDF);
  xx[2]      = findMaximum(2, psiWmPDF  , DpsiWmPDF);
  xx[3]      = findMaximum(3, phiWmPDF  , DphiWmPDF);
  
  fmax[1]    = jpsiphiWmPDF( xx , par );
  std::cout << "max{W-} = " << fmax[1] << std::endl;
  
}

double RandomDists::findMaximum( int vindex, Ptr2pdf pdf_func, Ptr2Dpdf Dpdf_func)
{

  double xsi   = 0.;
  int iter     = 0;
  int max_iter = 20;
  double x_lo  = 0.00;
  double x_hi  = 1.00;
  double froot = 0.9;
  
  //Run a quick scan on this function to determine a maximum
  int np = 20;
  double dx = (xmax[vindex] - xmin[vindex]) / (double)np;
  double xatmax=(0.);
  double yatmax=(0.);
  for( int i = 0; i < np; ++i) {
    xsi = 0.0+i*dx;
    double yy = pdf_func( &xsi , par);
    if ( yy > yatmax ) {
      yatmax = yy;
      xatmax = xsi;
    }
  }
  
  x_lo  = xatmax - dx;
  x_hi  = xatmax + dx;
  
  gsl_set_error_handler_off();
  
  T = gsl_root_fsolver_bisection;
  s = gsl_root_fsolver_alloc (T);
  
  F.function = Dpdf_func;
  F.params = params;
  
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  
  int status = GSL_CONTINUE;
  
  while (status == GSL_CONTINUE && iter < max_iter)
  {
    status = gsl_root_fsolver_iterate (s);
    froot = gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    
    status = gsl_root_test_interval (x_lo, x_hi,0, 0.00001);
    
    if (status) {
      if (status == -1 ) {
        fprintf (stderr, "failed \n");
      } 
    }
    
    if (status == GSL_SUCCESS) std::cout << " GSL: succeded after: " 
                                         << iter 
                                         << " iterations" 
                                         << std::endl;
    ++iter;
    
  }

  gsl_root_fsolver_free ( s );

  return froot;
  
}


void RandomDists::updateMaxima()
{  
  
  if (generatingBs) {
    pdfmax = fmax  [0];
  }
  else { 
    pdfmax = fmax  [1];
  }
  
}


void RandomDists::setParameters( double * par )
{
  
  params->gBar   = par[0];
  params->Dgamma = par[1];
  params->ft     = par[2];
  params->fp     = par[3];
  params->dt1    = 0.0;
  params->dt2    = 0.0;
  params->phis   = par[6];
  params->Dms    = par[7];
  params->omega  = par[8];
  params->time   = 0.0;
  
}


