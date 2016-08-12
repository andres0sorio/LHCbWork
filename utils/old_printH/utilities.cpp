#include "utilities.hpp"

std::string getBaseName(const char *fileName) {
  std::string fileNameString = std::string(fileName);
  int iStart = fileNameString.rfind('/')+1;  
  if(iStart < 0) iStart = 0;
  int iEnd = fileNameString.find('.',iStart);
  return fileNameString.substr(iStart,(iEnd-iStart));
}


/////////////////////////////////////////////////////////
// Print histograms in the same way as OutputOrganizer

void printHistograms(TDirectory *source) {
  
  TCanvas *plotarea;
  TPad *onepad;
  TString histoname;
  
  TStyle *st1 = new TStyle("st1", "Farah' style");
  setStyleOptions(st1);
  st1->cd();
  
  TDirectory *current_sourcedir = gDirectory;
  
  // Loop over all keys in this directory
  
  TIter nextkey (current_sourcedir->GetListOfKeys());
  TKey *key;
  
  while((key = (TKey*)nextkey())) {
    
    TObject *obj = key->ReadObj();
    
    if( obj->IsA()->InheritsFrom("TH1D")) {

#if _DEBUB
      std::cout << obj->GetName() << std::endl;
#endif
      
      plotarea = new TCanvas("plotarea","lhcAnalysis",0,0,100,80);
      onepad = new TPad("pad","lhcAnalysis",0.01,0.04,0.99,0.99,10);
      plotarea->cd();
      onepad->Draw();
      plotarea->Update();
      
      TH1D *h1 = (TH1D*)obj;

      //set histograms options: please edit
      h1->SetMarkerStyle(21);
      h1->SetMarkerSize(0.4);
      h1->SetLineColor(1);
      TAxis::TAxis *axis1;
      TAxis::TAxis *axis2;
      axis1 = h1->GetXaxis();
      axis2 = h1->GetYaxis();
      axis1->SetTitle(h1->GetName());
      axis2->SetTitle("Events");
      //this function is defined in plot_util.cpp
      setAxeOptions(axis1);
      setAxeOptions(axis2);
      
      histoname = TString(obj->GetName()) + TString(".eps");
      onepad->cd();
      h1->Draw();
      plotarea->Update();
      plotarea->Print(histoname);
      
      //only if you want to create gif files
      //      createGIF(TString(obj->GetName()));
      
      delete h1;
      delete plotarea;
      
    }    
    
    else if( obj->IsA()->InheritsFrom("TH2D")) {   
      
#if _DEBUB
      std::cout << obj->GetName() << std::endl;
#endif      
      plotarea = new TCanvas("plotarea","wwsAnalysis",0,0,100,80);
      onepad = new TPad("pad","wwsAnalysis",0.01,0.04,0.99,0.99,10);
      
      plotarea->cd();
      onepad->Draw();
      plotarea->Update();

      TH2D *h2 = (TH2D*)obj;
      
      //set histograms options: please edit
      h2->SetMarkerStyle(21);
      h2->SetMarkerSize(0.4);
      TAxis::TAxis *axis1;
      TAxis::TAxis *axis2;
      axis1 = h2->GetXaxis();
      axis2 = h2->GetYaxis();
      axis1->SetTitle(h2->GetName());
      axis2->SetTitle("Events");
      //this function is defined in plot_util.cpp
      setAxeOptions(axis1);
      setAxeOptions(axis2);
      
      histoname = TString(obj->GetName()) + TString(".eps");
      onepad->cd();
      h2->Draw();
      plotarea->Update();
      plotarea->Print(histoname);
      
      //only if you want to create gif files
      //createGIF(TString(obj->GetName()));
      
      delete h2;
      delete plotarea;
      
    } 
    
    else if ( obj->IsA()->InheritsFrom("TDirectory") ) {
      
#if _DEBUB
      std::cout << obj->GetName() << std::endl;
#endif     
      mkdir(obj->GetName(),S_IRWXU);
      chdir(obj->GetName());
      
      source->cd(obj->GetName());      
      TDirectory *newdir = gDirectory;
      printHistograms(newdir);
      
      chdir("../");
      
    }
    
    else std::cout << "printHistograms> What is this object?" << std::endl;
    
    if(obj) source->cd();
    
  }
  
  delete st1;
  
#if _DEBUG
  std::cout << "printHistograms : Done." << std::endl;
#endif 
  
}

