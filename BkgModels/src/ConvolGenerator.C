// $Id: ConvolGenerator.C,v 1.9 2007/02/24 15:19:55 aosorio Exp $
// Include files

// local
#include "ConvolGenerator.h"

//-----------------------------------------------------------------------------
// Implementation file for class : ConvolGenerator
//
// 2007-01-17 : Andres Osorio
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
ConvolGenerator::ConvolGenerator(  ) {}

ConvolGenerator::ConvolGenerator(  const char *filename , int max ) : RandomDists(filename, max) {
  
}

//=============================================================================
// Destructor
//=============================================================================
ConvolGenerator::~ConvolGenerator() {} 

//=============================================================================
void ConvolGenerator::initialise(double *pm, int npar)
{
  
  pi = TMath::Pi();
  
  for( int i = 0; i < npar; ++i) {
    par[i] = pm[i];
    std::cout << par[i] << std::endl;
  }

  xmin[0]= -2.0;  //ps
  xmax[0]=20.00;
  xmin[1]=-1.0; // cos(theta)
  xmax[1]= 1.0;
  xmin[2]=-1.0; // cos(psi)
  xmax[2]= 1.0;
  xmin[3]=-pi;  // phi
  xmax[3]= pi;
  
  for ( int i=0; i < 10; ++i) {
    var[i] = 0.0;
    tempoval[i] = 0.0;
  }
  
  pdfmax  = 0.07;
  
  //Common convolution parameters - default
  convpars[0] = 0.030;    //sigma_{propertime}
  //Double gaussian parameters
  convpars[1] = 0.801;     //f1
  convpars[2] = 0.0022;    //mu1
  convpars[3] = 0.0584;    //s1
  convpars[4] = -0.1114;   //mu2
  convpars[5] = 0.471;     //s2
  convpars[6] = 0.20;      //sigma for theta smearing
  
  
}

void ConvolGenerator::RunExperiment() 
{

  std::cout << "ConvolGenerator> RunExperiment starts now. \n";
  
  double u = 0.0;
  double f = 0.0;

  int counter = 0;
  double r1(0.0);
  double r2(0.0);
  
  counter = 0;

  astPDFwConv = &convolveBkg;
  
  while (counter < max_events) {
    
    tempoval[5] = qfactor;
    
    for ( int i = 0; i < 4; ++i ) {
      r1 =  rdn[i]->Rndm();
      tempoval[i] = xmin[i] + r1*(xmax[i]-xmin[i]);
    }
    
    tempoval[1] = acos( tempoval[1] ); // cos -> theta
    tempoval[2] = acos( tempoval[2] ); // cos -> psi
    
    r2 = reg->Rndm();
    u  = r2 * pdfmax;
    f  = (*astPDFwConv)( with2Gaussians, par, tempoval, convpars);
        
    if (f > pdfmax) {
      updateMaxima( f );
      u  = r2 * pdfmax;
    }
    
    if ( u < f ) {
      
      if ( counter < 5 )  data[4].push_back(tempoval[0]);
      
      for( int k=1; k < 4; ++k) {
        if ( counter < 5 )  data[k].push_back(tempoval[k]);
      } 
      
      WriteToTree();
      
      if( ((counter+1)%10000)==0 ) std::cout << "RunExperiment> " 
                                             << (counter+1) << std::endl;
      ++counter;
      
    }
        
  }
  
  if ( data[1].size() != 0 ) { 
    
    std::cout << "RunExperiment> event=0 " 
              << data[4][0] << '\t'
              << data[1][0] << '\t'
              << data[2][0] << '\t'
              << data[3][0] << std::endl;
  }
  
  t1->Write();
  
}

void ConvolGenerator::WriteToTree() 
{
  
  //Write to tree on a event by event basis - for large data sets
  var[0] = 0.0; //this is time without convolution
  var[1] = tempoval[1];//angle 1
  var[2] = tempoval[2];//angle 2
  var[3] = tempoval[3];//angle 3
  var[4] = tempoval[0];//lifetime w. convolution//
  var[5] = 0.0;//angle 1 w. smearing
  var[6] = 0.0;
  var[7] = 0.0;
  var[8] = 0.0;
  var[9] = qfactor;
  t1->Fill();
  
}
