// $Id: GetResults.C,v 1.8 2006/12/17 01:24:30 aosorio Exp $
// Include files 

// local
#include "GetResults.h"

//-----------------------------------------------------------------------------
// Implementation file for class : GetResults
//
// 2006-11-21 : Andres Felipe OSORIO OLIVEROS
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
GetResults::GetResults(  ) {

  fitres = new FitResults();

  npars = 0;
  
}

//=============================================================================
// Destructor
//=============================================================================
GetResults::~GetResults() {

  std::vector<TObject*>::iterator itr;
  itr = histograms.begin();
  while ( itr != histograms.end() ) {
    delete (*itr);
    ++itr;
  }

  std::vector<TF1*>::iterator itrf;
  itrf = fitfunctions.begin();
  while ( itrf != fitfunctions.end() ) {
    delete (*itrf);
    ++itrf;
  }
  
  for ( int k=0; k < npars; ++k) {
    delete plotarea[k];
  }
  
  if( fitres ) delete fitres;
  
} 

//=============================================================================

void GetResults::addDataFile( const char * infile ) //< Results from the fit
{
  
  std::string file = std::string(infile);
  fitresults.push_back(file);
  
}


void GetResults::addDataList( const char * infile ) //< Results from the fit
{
  
  ifstream in;
  in.open(infile);
  
  if(!in) {
    print_message("could not open input file.");
    print_message(infile);
    exit(1);
  }
  
  print_message("opening files:");
  
  //char comment[256];
  //in.getline (comment,256);
  //print_message(comment);
  
  while (1) {
    char str[60];
    in >> str;
    if (!in.good()) break;
    std::string file = std::string(str);
    int pos = file.find("results/", 0);
    if ( pos == 0 ) file = std::string("./") + file;
    else file = std::string("./results/") + file;
    fitresults.push_back(file);
  }    
  
  in.close();
  
}


void GetResults::addDataSpec( const char * infile ) 
{
  
  //Read input paremeters;
  readData(infile,tmcinputs);
  std::cout << "Total input parameters imported: " 
            << tmcinputs.size()
            << std::endl;
  
}

void GetResults::addHistogram( const char *name , int nbins, double xlo, double xhi)
{
  
  h1d = new TH1D( name , "Results from ToyMC experiments", nbins, xlo, xhi);
  h1d->Sumw2();
  histograms.push_back(h1d);
  
}


void GetResults::addFitFunction( const char *name, double xlo, double xhi)
{
  
  ff = new TF1( name , "gaus", xlo, xhi);
  fitfunctions.push_back(ff);
  
}


void GetResults::getFitResults()
{
  
  std::vector<std::string>::iterator itr;
  itr = fitresults.begin();
  
  while ( itr != fitresults.end() ) {
    
    fitres->AddResults( (*itr).c_str() );
    ++itr;
    
  }

  std::cout << "Total ToyMC experiments imported: " 
            << fitres->Size()
            << std::endl;
  
  std::cout << "Total number of fits which succeded: "
            << fitres->NSucceded()
            << std::endl;

  GetFailedExps();
  
}

void GetResults::GetFailedExps()
{
  
  std::vector<Parameters*>::iterator itpm;
  itpm = fitres->FirstParameters();
  int counter = 1;

  std::cout << "The following experiments failed to converge: ";
  
  while ( itpm != fitres->LastParameters() ) {
    
    if ( ! (*itpm)->has_succeded ) {
      std::cout << counter << '\t';
    }
    
    ++itpm;
    ++counter;
  }
  
  std::cout << '\n';
  
}

void GetResults::plotFittedValue( int param , int histo ) 
{

  gStyle->SetOptFit();
  
  std::vector<Parameters*>::iterator itpm;
  itpm = fitres->FirstParameters();
  while ( itpm != fitres->LastParameters() ) {
    
    if ( ! (*itpm)->has_succeded ) {
      ++itpm;
      continue;
    }
    
    double val = (*itpm)->param[param];
    h1d = (TH1D*) histograms[histo];
    h1d->Fill( val );
    ++itpm;
  }
  
  h1d->SetMarkerStyle(21);
  h1d->SetMarkerSize(0.7);
  h1d->Draw("e1p");
  
  ff = (TF1*) fitfunctions[histo];
  ff->SetLineWidth(2);
  ff->SetLineColor(1);
  h1d->Fit(ff,"0");
  ff->Draw("same");
  
}

void GetResults::plotFittedError( int param , int histo ) 
{
  
  std::vector<Parameters*>::iterator itpm;
  itpm = fitres->FirstParameters();
  while ( itpm != fitres->LastParameters() ) {

    if ( ! (*itpm)->has_succeded ) {
      ++itpm;
      continue;
    }

    double val = (*itpm)->errparam[param];
    h1d = (TH1D*) histograms[histo];
    h1d->Fill( val );
    ++itpm;
  }
  
  h1d->SetMarkerStyle(21);
  h1d->SetMarkerSize(0.7);
  h1d->Draw("e1p");

}

void GetResults::plotPull( int param , int histo ) 
{
  
  std::vector<Parameters*>::iterator itpm;
  itpm = fitres->FirstParameters();
  
  double x_true =  tmcinputs[param];
    
  while ( itpm != fitres->LastParameters() ) {

    if ( ! (*itpm)->has_succeded ) {
      ++itpm;
      continue;
    }
    
    double val    = (*itpm)->param[param];
    double errval = (*itpm)->errparam[param];
    double pull   = ( val - x_true ) / errval;
    h1d = (TH1D*) histograms[histo];
    h1d->Fill( pull );
    ++itpm;
  }
  
  h1d->SetMarkerStyle(21);
  h1d->SetMarkerSize(0.7);
  h1d->Draw("e1p");
  
  ff = (TF1*) fitfunctions[histo];
  ff->SetLineWidth(2);
  ff->SetLineColor(2);
  h1d->Fit(ff,"0");
  ff->Draw("same");
  
}

void GetResults::addSelection( int ipar, const char *name, double xmin, double xmax) 
{
  
  selpar.push_back( ipar );
  selparname.push_back( std::string(name) );
  selparmin.push_back( xmin );
  selparmax.push_back( xmax );
  
  char p1[20];
  sprintf( p1 ,"%s_error",name);
  char p2[20];
  sprintf( p2 ,"%s_pull",name);
  
  addHistogram(name,50 , xmin  , xmax);
  addFitFunction(name, xmin, xmax);
  addHistogram(p2  ,100, -8.0 , 8.0);
  addFitFunction(p2, xmin, xmax);
  
  //addHistogram(p1  ,50 , 0.001 , 0.010);

  ++npars;
  
}

void GetResults::drawHistograms() 
{

  std::cout << "GetResults> drawing results for "
            << npars
            << " parameters.\n";

  int hindex = 0;
  
  char canvasname[20];
  
  for ( int i = 0; i < npars; ++i ){
    
    int par = selpar[i];
    int xoffset = i*10;
    int yoffset = i*10;
  
    sprintf(canvasname,"plotare%d",i);
      
    plotarea[i] = new  TCanvas(canvasname,"Results from fits", 200+xoffset, 10+yoffset, 1200, 500);
    plotarea[i]->Divide(3,1);
    plotarea[i]->cd();
    
    plotarea[i]->cd(1);
    plotFittedValue(par, i+hindex);
    ++hindex;
    //plotarea[i]->cd(2);
    //plotPull       (par, i+hindex);
   
    //++hindex;
    //plotarea[i]->cd(3);
    //plotFittedError(par, i+hindex);

    plotarea[i]->Update();
    plotarea[i]->cd();
    
  }

}

void GetResults::evalGCorrel( int ipar )
{
  
  std::vector<ErrMatrix*>::iterator itpm;
  double rho(0.);
  int co(0);
  
  itpm = fitres->FirstErrMatrix();
  while ( itpm != fitres->LastErrMatrix() ) {
    
    if ( ! (*itpm)->has_errmatrix ) {
      ++itpm;
      continue;
    }
    
    //get the covariance/error matrix
    
    int maxparams = (*itpm)->maxparams;
    TMatrixD V(maxparams, maxparams);
    
    for(int i=0; i < maxparams; ++i){
      for(int j=0; j < maxparams; ++j){
        V[i][j] = sqrt( (*itpm)->matrix[i][j] );
      }
    }
    
    V.Print();
    
    std::cout << "Inverting matrix: " << std::endl;
    Double_t det       = V.Determinant();
    std::cout << det << std::endl;
    
    TMatrixD invV      = V.Invert();
    
    std::cout <<  V[ipar][ipar] << " " << invV[ipar][ipar] << std::endl;
    
    rho =  V[ipar][ipar]*invV[ipar][ipar];
    std::cout << "glob corr: " << rho << std::endl;
    
    ++co;
    ++itpm;
  
    break;
    
  }
  
  

}

//void GetResults::plotCorrel()
//{
//  
//
//
// 
//}
