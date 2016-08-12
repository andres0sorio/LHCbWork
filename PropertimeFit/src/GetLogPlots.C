// $Id: GetLogPlots.C,v 1.1 2007/04/09 11:49:08 aosorio Exp $
{

  TFile *f1 = new  TFile ("phis_study_dc04_0.root", "READ");
  
  if( f1 ) {
    std::cout << "getData> File: " << std::string("phis_study_dc04_0.root") 
              << " opened." << std::endl;
  }
  else {
    std::cout << "getData> Error: File: " << std::string("phis_study_dc04_0.root") 
              << " not accessible" << std::endl;
    exit(1);
  }
  
  TStyle *st1 = new TStyle();
  setStyleOptions(st1);
  st1->cd();
  
  f1->cd();
  
  TTree *T = new TTree();
  T = (TTree*)gROOT->FindObject("data");
  
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

  TFile * f2 = new TFile( "projections.root" , "RECREATE" );
  
  f2->cd();
  
  Projections * pr[6];
  pr[0] = new Projections("Propertime", 50, -2.0, 15.0 );
  //pr[1] = new Projections("Propertime_zoom", 50, -5.0, 20.0 );
  
  double time  = 0.0;
  double qfac  = 0.0;
  
  f1->cd();
  
  for(int i=0; i <  n_events; ++i){
    T->GetEntry(i);
    time= var[4];
    qfac = var[7];
    if ( qfac == 0.0) pr[0]->Initialise( time , false );
    //if ( qfac == 0.0) pr[1]->Initialise( time , false );
  }  
  delete T;
  f1->Close();
  delete f1;
  
  TF1 * fun[18];
  
  pr[0]->Initialise( 1.0 , true );
  //pr[1]->Initialise( 1.0 , true );

  fun[0] = new TF1("pdf1", properTimeWRes2,  -2.0, 15.0, 20);
  
  for( int k=0; k < 1; ++k)
    for ( int j=0; j < 1; ++j ) pr[k]->AddFunction( fun[k*3 + j], 20);
  
  for( int k = 0; k < 1; ++k ) {
    pr[k]->SetParameters( "initialSO_dc04_0.dat" ,  1.0 );
    pr[k]->Draw( );
  }
  
  
}
