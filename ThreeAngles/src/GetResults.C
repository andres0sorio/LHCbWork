// $Id: GetResults.C,v 1.6 2006/11/27 14:26:02 aosorio Exp $
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
  acanv = new TCanvas("acanc","Results from fits",200,10,1200,500);
  acanv->Divide(3,1);
  hpos = 1;
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
  
  if( fitres ) delete fitres;
  
  if( acanv ) delete acanv;
  
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
  
  acanv->cd(hpos);
  h1d->SetMarkerStyle(21);
  h1d->SetMarkerSize(0.7);
  h1d->Draw("e1p");
  acanv->Update();
  ++hpos;
  if ( hpos > 4 ) hpos = 1;
  acanv->cd();

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
  
  acanv->cd(hpos);
  h1d->SetMarkerStyle(21);
  h1d->SetMarkerSize(0.7);
  h1d->Draw("e1p");
  acanv->Update();
  ++hpos;
  if ( hpos > 4 ) hpos = 1;
  acanv->cd();
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
  
  acanv->cd(hpos);
  h1d->SetMarkerStyle(21);
  h1d->SetMarkerSize(0.7);
  h1d->Draw("e1p");
  acanv->Update();
  ++hpos;
  if ( hpos > 4 ) hpos = 1;
  acanv->cd();
}

