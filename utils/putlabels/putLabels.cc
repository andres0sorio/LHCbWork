// $Id: putLabels.cc,v 1.2 2006/03/26 16:10:20 aooliver Exp $
// Include files 



// local
#include "putLabels.h"

//-----------------------------------------------------------------------------
// Implementation file for class : putLabels
// little class to read histograms from Gaudi and add labels

// 2006-03-26 : Andres Osorio Oliveros
//-----------------------------------------------------------------------------
std::fstream& operator>>(std::fstream &istr, Options &rhs) {
  istr >> rhs.id;
  istr >> rhs.xLabel;
  istr >> rhs.yLabel; 
  return istr;
}

std::ostream& operator<<(std::ostream &ostr, Options &rhs) {
  ostr << rhs.id << '\t';
  ostr << rhs.xLabel << '\t';
  ostr << rhs.yLabel << '\n';
  ostr << std::endl;
  return ostr;
}

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
putLabels::putLabels( const char *inopts, const char *infile ) {
  
  opts = new std::fstream(inopts,std::ios_base::in);
    
  if(! opts->is_open()) {
    std::cout << "putLabels> could not open input file> " 
              << inopts
              << std::endl;
    exit(1);
  }
  
  while(1) {
          
    Options *temp = new Options();
    (*opts) >> *temp;
    if ( !opts->eof() ) {
      vopts.push_back(temp);}
    else { 
      delete temp;
      break;
    }
    
  }

  opts->close();
  delete opts;
  
  input = new TFile(infile,"READ");
  
  if (!input) { 
    std::cout << "putLabels> could not open input file " 
              << infile
              << std::endl;
    exit(1);
  }
  
  std::cout << "putLabels> New file will be output.root "
	    << std::endl;

  output = new TFile("output.root","RECREATE");
  
}

//=============================================================================
// Destructor
//=============================================================================
putLabels::~putLabels() {

  if (input) { input->Close(); delete input; }

  if (output) { output->Close(); delete output; }

} 

//=============================================================================
void putLabels::processFile()
{

  loopOverHistograms(output,input);

}

inline void putLabels::setOptions(TH1D *histo){

  std::vector<Options*>::iterator itr;
  itr = vopts.begin();
  
  while (itr != vopts.end()){

    int hid = (*itr)->id;
    char hname[10];
    sprintf(hname,"%d",hid); //temporary - histograms have no name but a number
    
    TString histo_name(histo->GetName());
    
    if( histo_name == TString(hname)) {
      TAxis::TAxis *axis1;
      axis1 = histo->GetXaxis();
      axis1->SetTitle( (*itr)->xLabel.c_str() );
      axis1 = histo->GetYaxis();
      axis1->SetTitle( (*itr)->yLabel.c_str() );
      
    }
    
    ++itr;
    
  }
  
  
}

inline void putLabels::setOptions(TH2D *histo){



}


void putLabels::loopOverHistograms(TDirectory *target, TDirectory *source) 
{
  
  TString path( (char*)strstr( target->GetPath(), ":"));
  path.Remove(0,2);
  
  source->cd(path);
  TDirectory *current_sourcedir = gDirectory;
  
  // Loop over all keys in this directory
  TIter nextkey (current_sourcedir->GetListOfKeys());
  TKey *key;
  
  while ((key = (TKey*)nextkey())) {
    
    source->cd(path);
    
    TObject *obj = key->ReadObj();
    
    if( obj->IsA()->InheritsFrom("TH1D")) {
      
      TH1D *h1 = (TH1D*)obj;
      setOptions(h1);
    }
    
    else if( obj->IsA()->InheritsFrom("TH2D")) {
      
      TH2D *h1 = (TH2D*)obj;
      setOptions(h1);
    }
    
    else if ( obj->IsA()->InheritsFrom("TDirectory") ) {
      
      target->cd();
      TDirectory *newdir = target->mkdir (obj->GetName(), obj->GetTitle() );
      loopOverHistograms(newdir, source);
      
    } 
    
    else {
    }
    
    if(obj) {
      target->cd();
      obj->Write( key->GetName() );
    }
    
  }
  
  target->Write();
  
}
 
