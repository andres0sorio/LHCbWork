// $Id: printHistograms.cc,v 1.2 2006/07/28 13:16:04 aooliver Exp $
// Include files 



// local
#include "root_utilities.h"

//-----------------------------------------------------------------------------
// Implementation file for class : printHistograms
//
// 2006-03-03 : Andres Osorio Oliveros
//-----------------------------------------------------------------------------

//=============================================================================

// Print histograms 

void printHistograms(TDirectory *source) {
  
  TCanvas *plotarea;
  TPad *onepad;
  TString histoname;
  
  TStyle *st1 = new TStyle("st1", "General output style");
  setStyleOptions(st1);
  st1->cd();
  
  TDirectory *current_sourcedir = gDirectory;
  
  // Loop over all keys in this directory
  
  TIter nextkey (current_sourcedir->GetListOfKeys());
  TKey *key;
  
  while((key = (TKey*)nextkey())) {
    
    TObject *obj = key->ReadObj();
    
    if( obj->IsA()->InheritsFrom("TH1")) {

      TH1D *histo = (TH1D*)obj;
            
      print_debug_message(obj->GetName());
      
      plotarea = new TCanvas("plotarea","General output",0,0,100,100);
      onepad = new TPad("pad","General output",0.01,0.04,0.99,0.99,10);
      plotarea->cd();
      onepad->Draw();
      plotarea->Update();
      
      setHistogramsOptions(histo);
            
      histoname = TString(obj->GetName()) + TString(".eps");
      onepad->cd();
      histo->Draw("hist");
      plotarea->Update();
      plotarea->Print(histoname);
      
      delete plotarea;
      
    }    
    
    else if( obj->IsA()->InheritsFrom("TH2")) {   
      
      print_debug_message(obj->GetName());
      
      plotarea = new TCanvas("plotarea","General output",0,0,100,80);
      onepad = new TPad("pad","General output",0.01,0.04,0.99,0.99,10);
      
      plotarea->cd();
      onepad->Draw();
      plotarea->Update();
      
      histoname = TString(obj->GetName()) + TString(".eps");
      onepad->cd();
      obj->Draw();
      plotarea->Update();
      gErrorIgnoreLevel = 1;
      plotarea->Print(histoname);
      
      delete plotarea;
      
    } 
    
    else if ( obj->IsA()->InheritsFrom("TDirectory") ) {
      
      print_debug_message(obj->GetName());
      
      mkdir(obj->GetName(),S_IRWXU);
      chdir(obj->GetName());
      
      source->cd(obj->GetName());      
      TDirectory *newdir = gDirectory;
      printHistograms(newdir);
      
      chdir("../");
      
    }
    
    else print_message("printHistograms> What is this object?");
    
    if(obj) source->cd();
    
  }
  
  delete st1;
  
  print_debug_message("printHistograms : Done.");

}
