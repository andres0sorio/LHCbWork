// $Id: RandomDists.C,v 1.6 2007/02/24 15:19:56 aosorio Exp $
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
  
  //dist[0]  = new TH1D("h1","This is the total distribution",100,0.0,20.0);
  //dist[0]->Sumw2();
    
  max_events = 1000;
  
  std::cout << "RandomDists> the dafault MaxEvents is: " 
            << max_events << std::endl;
  
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
  var[9] = 0.0;
  
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
  t1->Branch("tag"       ,&var[9] ,"tag/D");
  
  pi = TMath::Pi();
  
  rdn[0] = new TRandom3(0);
  rdn[1] = new TRandom3(0);
  rdn[2] = new TRandom3(0);
  rdn[3] = new TRandom3(0);
  reg    = new TRandom3(0);
  tag[0] = new TRandom3(0);
  tag[1] = new TRandom3(0);
  tag[2] = new TRandom3(0);
  tag[3] = new TRandom3(0);

}

RandomDists::RandomDists( const char *filename , int max ) {

  max_events = max;
  
  //dist[0]  = new TH1D("h1","This is the total distribution",100,0.0,20.0);
  //dist[0]->Sumw2();
  
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
  t1->Branch("tag"       ,&var[9] ,"tag/D");
  
  rdn[0] = new TRandom3(0);
  rdn[1] = new TRandom3(0);
  rdn[2] = new TRandom3(0);
  rdn[3] = new TRandom3(0);
  reg    = new TRandom3(0);
  tag[0] = new TRandom3(0);
  tag[1] = new TRandom3(0);
  tag[2] = new TRandom3(0);
  tag[3] = new TRandom3(0);
  
}


//=============================================================================
// Destructor
//=============================================================================
RandomDists::~RandomDists() 
{
  
  //for(int i=0; i <8; ++i){
  //  if(dist[i]) {
  //    dist[i]->Delete();
  //  }
  //  data[i].clear();
  //}
  
  if(outfile) {
    outfile->Close();
    delete outfile;
  }
  
  for (int i=0; i < 4; ++i) delete rdn[i];
  delete reg;
  for (int i=0; i < 4; ++i) delete tag[i];

    
} 

void RandomDists::initialise(double *pm, int npar)
{
  
  pi = TMath::Pi();
  
  for( int i = 0; i < npar; ++i) {
    par[i] = pm[i];
    std::cout << par[i] << std::endl;
  }

  xmin[0]=0.0;  //ps
  xmax[0]=20.00;
  xmin[1]=-1.0; // cos(theta)
  xmax[1]= 1.0;
  xmin[2]=-1.0; // cos(psi)
  xmax[2]= 1.0;
  xmin[3]=-pi;  // phi
  xmax[3]= pi;
  
  for ( int i=0; i < 8; ++i) {
    var[i] = 0.0;
  }
  
  for (int i=0; i< 4; ++i) {
    tempoval[i] = 0.0;
  }
  
  pdfmax  = 0.05;
  
  bkg_ov_sig = 1.00; //<- B/S ratio
    
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
    
    SelectSignalType();
    tempoval[5] = qfactor;
    
    for ( int i = 0; i < 4; ++i ) {
      r1 =  rdn[i]->Rndm();
      tempoval[i] = xmin[i] + r1*(xmax[i]-xmin[i]);
    }
    
    tempoval[1] = acos( tempoval[1] ); // cos -> theta
    tempoval[2] = acos( tempoval[2] ); // cos -> psi
    
    r2 = reg->Rndm();
    u  = r2 * pdfmax;
    f  = (*astPDF)(tempoval , par); // <- evaluate the pdf
    
    if (f > pdfmax) {
      updateMaxima( f );
      u  = r2 * pdfmax;
    }
    
    if ( u < f ) {
      
      for( int k=0; k < 4; ++k) {
        if ( counter < 5 )  data[k].push_back(tempoval[k]);
      }
      
      WriteToTree();
      
      if( ((counter+1)%10000)==0 ) std::cout << "RunExperiment> " 
                                             << (counter+1) 
                                             << std::endl;
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

void RandomDists::SelectSignalType( )
{

  astPDF   = &bkgPDF;
  qfactor  = 0;
      
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
  var[9] = qfactor;
  t1->Fill();
  
}

void RandomDists::updateMaxima( double nmax )
{  
  
  fmax  [0] = nmax;
  std::cout << "updateMaxima: new maximum:"
            << fmax[0]
            << std::endl;
  pdfmax = fmax  [0];
    
}

