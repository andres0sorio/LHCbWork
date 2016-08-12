// $Id: TestToyMC.C,v 1.1 2006/11/24 22:15:17 aosorio Exp $
// Include files 



// local
#include "TestToyMC.h"

//-----------------------------------------------------------------------------
// Implementation file for class : TestToyMC
//
// 2006-11-23 : Andres Felipe OSORIO OLIVEROS
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TestToyMC::TestToyMC( const char *infile ) {


  hh[0] = new TH1D("Transversity","Distributions",100,-1.0,1.0);
  hh[1] = new TH1D("phi","Distributions",100,-3.15,3.15);
  hh[2] = new TH1D("Cos(psi)","Distributions",100,-1.0,1.0);
  
  GetData( infile );
  
  acanv = new TCanvas("acanc","Results from fits",200,10,500,500);

  maxhistos = 3;
  
}
//=============================================================================
// Destructor
//=============================================================================
TestToyMC::~TestToyMC() {

  for ( int i=0; i < maxhistos; ++i) delete hh[i];
  delete T;
  
  f1->Close();
    
  delete f1;

  delete acanv;
  

}

void TestToyMC::DrawDistributions()
{
  
  int maxevts = T->GetEntriesFast();
    
  for ( int i = 0; i < (int) maxevts; ++i) {
    
    T->GetEntry(i);
    
    double theta = var[1];
    double psi   = var[2];
    double phi   = var[3];
    
    hh[0]->Fill( TMath::Sin(psi)*TMath::Sin(phi));
    hh[1]->Fill( phi );
    hh[2]->Fill( TMath::Cos( psi ) );
  
  }
  

  
  acanv->cd();
  hh[1]->SetMarkerStyle(8);
  hh[1]->SetMarkerSize(0.7);
  hh[1]->Draw("P");
  
}


//=============================================================================
void TestToyMC::GetData(const char *infile) {
  
  f1 = new  TFile (infile, "READ");
  
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
  
  T = new TTree();
  
  T = (TTree*)gROOT->FindObject("data");
  
 
  
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
  
  
}
