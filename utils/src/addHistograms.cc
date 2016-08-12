// $Id: addHistograms.cc,v 1.2 2006/07/28 13:16:04 aooliver Exp $
// Include files 



// local
#include "root_utilities.h"

//-----------------------------------------------------------------------------
// Implementation file for class : addHistograms
//
// 2006-03-03 : Andres Osorio Oliveros
//-----------------------------------------------------------------------------

//////////////////////////////////////////////////////////
//// Add 1D histograms

void addRootHistos(TDirectory *target, TList *sourcelist) {
  
  TString path( (char*)strstr( target->GetPath(), ":"));
  path.Remove(0,2);
  
  TFile *first_source = (TFile*)sourcelist->First();
  first_source->cd(path);
  TDirectory *current_sourcedir = gDirectory;
  
  // Loop over all keys in this directory
  TChain *globChain = 0;
  TIter nextkey (current_sourcedir->GetListOfKeys());
  TKey *key;
  
  while ((key = (TKey*)nextkey())) {
    
    first_source->cd(path);
    
    TObject *obj = key->ReadObj();
    
    if( obj->IsA()->InheritsFrom("TH1D")) {
      
      TH1D *h1 = (TH1D*)obj;
      // loop over all source files
      TFile *nextsource = (TFile*)sourcelist->After(first_source);
      
      while (nextsource) {
        
        if(nextsource->cd(path)) {
          TH1D *h2 = (TH1D*)gDirectory->Get(h1->GetName());
          if( h2 ) {
            h1->Add ( h2 );
            delete h2;
          }
        }
        
        else print_message("This path doesn't exist: ", path.Data());
        
        nextsource = (TFile*)sourcelist->After(nextsource);
      }
    }
    
    else if( obj->IsA()->InheritsFrom("TH2D")) {
      
      TH2D *h1 = (TH2D*)obj;
      // loop over all source files
      TFile *nextsource = (TFile*)sourcelist->After(first_source);
      
      while (nextsource) {
        if(nextsource->cd(path)) {
          TH2D *h2 = (TH2D*)gDirectory->Get(h1->GetName());
          if( h2 ) {
            h1->Add ( h2 );
            delete h2;
          }
        }
        
        else print_message("This path doesn't exist: ", path.Data());
	
        nextsource = (TFile*)sourcelist->After(nextsource);
      }
    }
    
    else if ( obj->IsA()->InheritsFrom("TDirectory") ) {
      
      target->cd();
      TDirectory *newdir = target->mkdir (obj->GetName(), obj->GetTitle() );
      
      addRootHistos(newdir, sourcelist);
      
    } 

    else if(obj->IsA()->InheritsFrom( "TTree" )) {
      std::cout << "Tree-Key name..." <<  key->GetName() << "...." << std::endl;
      globChain->Write( key->GetName() );
      
    } 

    else {
      print_message("Unknown object type, name: ");
      print_message(obj->GetName(),obj->GetTitle());
    }
    
    // now write merged histograms
    
    if (obj) {
      
      target->cd();
      
      //!!if the object is a tree, it is stored in globChain...
      //if(obj->IsA()->InheritsFrom( "TTree" ))
      //globChain->Write( key->GetName() );
      //else obj->Write( key->GetName() );

      obj->Write( key->GetName() );
      
    }
    
  }
  
  target->Write();
  
}
