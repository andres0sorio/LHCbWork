// $Id: combHistograms.cc,v 1.1 2006/03/03 17:20:20 aooliver Exp $
// Include files 



// local
#include "root_utilities.h"
//-----------------------------------------------------------------------------
// Implementation file for class : combHistograms
//
// 2006-03-03 : Andres Osorio Oliveros
//-----------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////
// Combine histograms

void combineHistograms( TDirectory *basesource, TList *sourcelist) {
  
  TString path( (char*)strstr( basesource->GetPath(), ":"));
  path.Remove(0,2);
  
  TDirectory *current_sourcedir = gDirectory;
  
  //--------------------------------------------------------------------//
  TCanvas *plotarea;
  TPad *onepad;
  TString histoname;
  TLegend *legend;
  float xpos;
  std::string tfileName;
  std::string processName;
  TString legendItem;
  //-------------------------------------------------------------------//
  
  
  // Loop over all keys in this directory
  
  TIter nextkey (current_sourcedir->GetListOfKeys());
  TKey *key;
  
  while((key = (TKey*)nextkey())) {
    
    current_sourcedir->cd();

    TObject *obj = key->ReadObj();
    
    print_debug_message("combineHistograms : location 1 : ");
    print_debug_message(obj->GetName());
    print_debug_message(gDirectory->GetPath());
    
    if( obj->IsA()->InheritsFrom("TH1D")) {
      

      // Add legend to histograms
      xpos = 0.4;
      legend = new TLegend(xpos,0.70,xpos+0.30,0.95);
      setLegendOptions(legend);
      tfileName = getTFileName(gDirectory->GetPath());
      processName = getProcessName(tfileName.c_str()); 
      legendItem = TString(processName.c_str());
      legend->AddEntry(obj, legendItem,"p");
      
      ////////////////////////////////////
      
      plotarea = new TCanvas("plotarea","General output",0,0,100,80);
      onepad = new TPad("pad","General output",0.01,0.04,0.99,0.99,10);
      plotarea->Draw();
      plotarea->cd();
      onepad->Draw();
      plotarea->Update();
      
      histoname = TString(obj->GetName()) + TString(".eps");
      onepad->cd();
      obj->Draw();
      
      TFile *nextsource = (TFile*)sourcelist->After(sourcelist->First());
      
      TH1D *h2 = new TH1D();
      
      while (nextsource) {
	
	if(nextsource->cd(path)) {
	  
	  h2 = (TH1D*)gDirectory->Get(obj->GetName());
	  
	  if( h2 ) {
	    
	    print_debug_message("combineHistograms : location 2 : ");
	    print_debug_message(obj->GetName());
	    print_debug_message(gDirectory->GetPath());
	    
	    tfileName = getTFileName(gDirectory->GetPath());
	    processName = getProcessName(tfileName.c_str()); 
	    legendItem = TString(processName.c_str());
	    legend->AddEntry(obj, legendItem,"p");
	    h2->Draw("same");
	  }
	}
	nextsource = (TFile*)sourcelist->After(nextsource);
      }
      
      legend->Draw();
      plotarea->Update();
      gErrorIgnoreLevel = 1;
      plotarea->Print(histoname);
      
      delete h2;
      delete plotarea;
      delete legend;
    }    
      
    else if( obj->IsA()->InheritsFrom("TH2D")) {   
      
      
      // Add legend to histograms
      xpos = 0.4;
      legend = new TLegend(xpos,0.70,xpos+0.30,0.95);
      setLegendOptions(legend);
      tfileName = getTFileName(gDirectory->GetPath());
      processName = getProcessName(tfileName.c_str()); 
      legendItem = TString(processName.c_str());
      legend->AddEntry(obj, legendItem,"p");
      ////////////////////////////////////
      
      plotarea = new TCanvas("plotarea","General output",0,0,100,80);
      onepad = new TPad("pad","General output",0.01,0.04,0.99,0.99,10);
      
      plotarea->cd();
      onepad->Draw();
      plotarea->Update();
      
      histoname = TString(obj->GetName()) + TString(".eps");
      onepad->cd();
      obj->Draw();
      
      TFile *nextsource = (TFile*)sourcelist->After(sourcelist->First());
      
      TH2D *h2 = new TH2D();
      
      while (nextsource) {
	if(nextsource->cd(path)) {
	  
	  h2 = (TH2D*)gDirectory->Get(obj->GetName());
	  
	  if( h2 ) {
	    
	    print_debug_message("combineHistograms : location 2 : ");
	    print_debug_message(obj->GetName());
	    print_debug_message(gDirectory->GetPath());
	    
	    tfileName = getTFileName(gDirectory->GetPath());
	    processName = getProcessName(tfileName.c_str()); 
	    legendItem = TString(processName.c_str());
	    legend->AddEntry(h2, legendItem,"p");
	    h2->Draw("same");
	  }
	}
	nextsource = (TFile*)sourcelist->After(nextsource);
      }
      legend->Draw();
      plotarea->Update();
      gErrorIgnoreLevel = 1;
      plotarea->Print(histoname);
      
      delete h2;
      delete plotarea;
      delete legend;
    } 
      
    else if ( obj->IsA()->InheritsFrom("TDirectory") ) {
      
      print_debug_message("combineHistograms : location 3 : ");
      print_debug_message(obj->GetName());
      print_debug_message(gDirectory->GetPath());
      
      mkdir(obj->GetName(),S_IRWXU);
      chdir(obj->GetName());
      
      basesource->cd(obj->GetName());
      TDirectory *newdir = gDirectory;
      
      combineHistograms(newdir, sourcelist);
      
      chdir("../");
      
    }
    
    else std::cout << "combineHistograms> What is this object?" << std::endl;
    
  }
  
  print_debug_message("combineHistograms : Done.");
  
}

////////////////////////////////////////////////////////////////////
// Combine histograms - Signal and Background

void combineSignalBack( TDirectory    *basesource , 
			TList         *sourcelist , 
			Int_t          option     , 
			std::ifstream *afile       ) 
{
  
  TString path( (char*)strstr( basesource->GetPath(), ":"));
  path.Remove(0,2);
  
  TDirectory *current_sourcedir = gDirectory;
  
  /////////////////////////////////////
  //TObjArray *signal = NULL;
  std::vector<TObjArray *> signal;
  std::vector<TObjArray *> background;
  
  bool isSignal = false;
  bool isBack = false;
  bool isdone = false;
  
  /////
  TCanvas *plotarea;
  TPad *onepad;
  TString histoname;
  TLegend *legend;
  Int_t backgroundColours[10] = {1,3,6,41,42,45,50,33,32,7};
  float xpos;
  std::string tfileName;
  std::vector<std::string> processName_signal;
  std::vector<std::string> processName_background;
  TString legendItem;
  ///////////////////////////////////////////
  //Loop over all keys in this directory
  
  TIter nextkey (current_sourcedir->GetListOfKeys());
  TKey *key;
  
  while((key = (TKey*)nextkey())) {
    
    current_sourcedir->cd();
    
    TObject *obj = key->ReadObj();
    
    if ( obj->IsA()->InheritsFrom("TDirectory") ) {
  
      print_debug_message(obj->GetName());
    
      basesource->cd(obj->GetName());
      TDirectory *newdir = gDirectory;
      
      if (isDirNamed("nunuWW",obj->GetName())) 
	{
	  isSignal = true;
	  signal.push_back(readHistograms(newdir));
	  tfileName = getTFileName(gDirectory->GetPath());
	  //processName_signal.push_back(getProcessName(tfileName.c_str())); 
	  processName_signal.push_back(std::string("WW#nu#nu (signal)")); 
	  
	}
      
      else if(isDirNamed("nunuZZ",obj->GetName())) 
	{
	  isSignal = true;
	  signal.push_back(readHistograms(newdir));
	  tfileName = getTFileName(gDirectory->GetPath());
	  //processName_signal.push_back(getProcessName(tfileName.c_str()));
	  processName_signal.push_back(std::string("ZZ#nu#nu (signal)"));
	}
      
      else if (isDirNamed("Background",obj->GetName()))
	{
	  isBack = true;
	  
	  background.push_back(readHistograms(newdir));
	  tfileName = getTFileName(gDirectory->GetPath());
	  processName_background.push_back(std::string("qqqq#nu#nu"));
	  
	  TFile *nextsource = (TFile*)sourcelist->After(sourcelist->First());
	  TString path( (char*)strstr( newdir->GetPath(), ":"));
	  path.Remove(0,2);
	  
	  TDirectory *extdir;
	  
	  while(nextsource) {
	    if(nextsource->cd(path)) {
	      extdir = gDirectory;
	      background.push_back(readHistograms(extdir));
	      tfileName = getTFileName(gDirectory->GetPath());
	      processName_background.push_back(getProcessName(tfileName.c_str())); 
	    }
	    nextsource = (TFile*)sourcelist->After(nextsource);
	  }
	}
      
      else {
	
	mkdir(obj->GetName(),S_IRWXU);
	chdir(obj->GetName());
	
	combineSignalBack(newdir, sourcelist, option, afile);
	
	chdir("../");
      }
      
    }
    
    //////////////////////////////////////////////////////////
    //Plot histograms

    if( isSignal && isBack && !isdone) {
      
      plotarea = NULL;
      onepad = NULL;
      legend = NULL;
      Color_t color = 0;

      char cpName[50];
      
      TObject *histos1 = signal[0]->First();
      
      while(histos1) {
	
	TObjArray *inputH1markers = new TObjArray();
	TObjArray *inputH1lines = new TObjArray();
	
	std::string h1name = std::string(histos1->GetName());
	
	plotarea = new TCanvas("plotarea","General output",0,0,100,80);
	onepad = new TPad("pad","General output",0.01,0.01,0.99,0.99,10);
	//onepad = new TPad("pad","General output",0.01,0.04,0.99,0.99,10);//only if date is on
	plotarea->cd();
	onepad->Draw();
	if(option == 1) onepad->SetLogy(1);
	else onepad->SetLogy(0);
	plotarea->Update();
	
	if( histos1->IsA()->InheritsFrom("TH1D") ) {
	  
	  TH1D *h1 = (TH1D*)histos1;
	  ////////////////////////////////////
	  // Add legend to histograms
	  xpos = findBestQuadrant(h1)*0.10+0.15;
	  legend = new TLegend(xpos,0.70,xpos+0.30,0.95);
	  setLegendOptions(legend);
	  legendItem = TString(processName_signal[0].c_str());
	  legend->AddEntry(histos1, legendItem,"p");
	  ////////////////////////////////////
	  
	  histoname = TString(histos1->GetName()) + TString(".eps");
	  setXaxisTitle(h1,afile);
	  color = h1->GetLineColor();
	  TH1D *h2 = (TH1D*)histos1->Clone();
	  h2->SetName("copyOne");
	  h2->SetFillColor(color);
	  h2->SetFillStyle(3004);
	  setXaxisTitle(h2,afile);
	  inputH1markers->AddLast(h1);
	  inputH1lines->AddLast(h2);
	  
	  std::vector<TObjArray*>::iterator itr;
	  
	  int k = 1;
	  
	  plotarea->Update();
	  
	  itr = signal.begin();
	  ++itr;
	  
	  for(itr=itr; itr != signal.end(); ++itr) {
	    
	    TObject *histos2 = (*itr)->First();
	    
	    legendItem = TString(processName_signal[k].c_str());
	    
	    while(histos2) {
	      
	      std::string h2name = std::string(histos2->GetName());
	      
	      if(h1name == h2name) {
		
		sprintf(cpName,"copyNameb_%d",k);
		
		TH1D *h1b = (TH1D*)histos2;
		setXaxisTitle(h1b,afile);
		color = h1b->GetLineColor();
		legend->AddEntry(histos2, legendItem,"p");
		TH1D *h2b = (TH1D*)h1b->Clone();
		h2b->SetName(cpName);
		h2b->SetFillColor(color);
		h2b->SetFillStyle(3004);
		setXaxisTitle(h2b,afile);
		inputH1markers->AddLast(h1b);
		inputH1lines->AddLast(h2b);
		
	      }
	      histos2 = (*itr)->After(histos2);
	    }
	    k++;
	  }
	  
	  
	  k=0;
	  
	  for(itr = background.begin(); itr != background.end(); ++itr) {
	    
	    TObject *histos2 = (*itr)->First();
	    legendItem = TString(processName_background[k].c_str()) + TString(" (background)");
	    
	    while(histos2) {
	      
	      std::string h2name = std::string(histos2->GetName());
	      
	      if(h1name == h2name) {
		sprintf(cpName,"copyNamec_%d",k);
		color = backgroundColours[k];
		TH1D *h1c = (TH1D*)histos2;
		setXaxisTitle(h1c,afile);
		h1c->SetMarkerColor(color);
		legend->AddEntry(h1c, legendItem,"p");
		TH1D *h2c = (TH1D*)h1c->Clone();
		h2c->SetName(cpName);
		setXaxisTitle(h2c,afile);
		h2c->SetFillColor(color);
		h2c->SetLineColor(color);
		h2c->SetFillStyle(3003);
		inputH1markers->AddLast(h1c);
		inputH1lines->AddLast(h2c);
		      
	      }
	      
	      histos2 = (*itr)->After(histos2);
	    }
	    k++;
	  }
	  
	  //////////////////////////////////////////////////////////
	  // Sort histograms and plot them acording to area criteria
	
	  TObjArray *sortedH1_markers;
	  TObjArray *sortedH1_lines;

	  sortedH1_markers = sortH1(inputH1markers, option);
	  sortedH1_lines = sortH1(inputH1lines, option);
	  
	  TObject *obj1 = sortedH1_markers->First();
	  TObject *obj2 = sortedH1_lines->First();
	  
	  TString option1 = TString("");
	  TString option2 = TString("histsame");

	  onepad->cd();
	  	  
	  ////////////////////
	  /////Loop over all the 1D histograms
	  /////Draw 
	  
	    while(obj1) {
	      std::cout << "combineHistograms> drawing ... " 
			<< std::string(obj1->GetName())
			<< std::endl;
	      obj1->Draw(option1);
	      obj2->Draw(option2);
	      obj1 = sortedH1_markers->After(obj1);
	      obj2 = sortedH1_lines->After(obj2);
	      option1 = TString("same");
	    }
	    
	    legend->Draw();
	    plotarea->Update();
	    gErrorIgnoreLevel = 1;
	    plotarea->Print(histoname);
	
	    print_message("combineHistograms> print Histogram ... done!");
	    
	    //createGIF(TString(histos1->GetName()));
	    
	    delete plotarea;
	    delete legend;
	    delete inputH1markers;
	    delete inputH1lines;
	}
	
	//////////////////////////////////////////////////////////////
	/// 2D histograms - Draw()
	
	if( histos1->IsA()->InheritsFrom("TH2D")) {
	  
	  // Add legend to histograms
	  xpos = 0.4;
	  legend = new TLegend(xpos,0.70,xpos+0.30,0.95);
	  setLegendOptions(legend);
	  legendItem = TString(processName_signal[0].c_str());
	  legend->AddEntry(histos1, legendItem,"f");
	  ////////////////////////////////////
	  
	  histoname = TString(histos1->GetName()) + TString(".eps");
	  onepad->cd();
	  onepad->SetLogy(0);
	  
	  TH2D *h2 = (TH2D*)histos1->Clone();
	  //setXaxisTitle(h2,afile);
	  color = h2->GetMarkerColor();
    // set2dHistoOptions(h2,color);
	  h2->SetName("copyOne");
	  h2->Draw();
	  
	  std::vector<TObjArray*>::iterator itr;
	  
	  int k = 1;
	  
	  plotarea->Update();
	  
	  itr = signal.begin();
	  ++itr;
	  
	  for(itr=itr; itr != signal.end(); ++itr) {
	    
	    TObject *histos2 = (*itr)->First();
	    
	    legendItem = TString(processName_signal[k].c_str());
	    
	    while(histos2) {
	      
	      std::string h2name = std::string(histos2->GetName());
	      
	      if(h1name == h2name) {
		
		sprintf(cpName,"copyNameb_%d",k);
		legend->AddEntry(histos2, legendItem,"f");
		TH2D *h2b = (TH2D*)histos2->Clone();
		//setXaxisTitle(h2b,afile);
		color = h2b->GetMarkerColor();
		//set2dHistoOptions(h2b,color);
		h2b->SetName(cpName);
		h2b->Draw("same");
		
	      }

	      histos2 = (*itr)->After(histos2);
	    }
	    k++;
	  }
	  
	  plotarea->Update();
	  k=0;

	  for(itr = background.begin(); itr != background.end(); ++itr) {
	    
	    TObject *histos2 = (*itr)->First();
	    legendItem = TString(processName_background[k].c_str()) + TString(" (background)");
	    
	    while(histos2) {
	      
	      std::string h2name = std::string(histos2->GetName());
	      
	      if(h1name == h2name) {
		sprintf(cpName,"copyNamec_%d",k);
		legend->AddEntry(histos2, legendItem,"f");
		TH2D *h2c = (TH2D*)histos2;
    //	setXaxisTitle(h2c,afile);
		color = backgroundColours[k];
		//set2dHistoOptions(h2c,color);
		h2c->SetName(cpName);
		h2c->Draw("same");
	      }
	      
	      histos2 = (*itr)->After(histos2);
	    }
	    k++;
	  }
	  	    
	  legend->Draw();
	  plotarea->Update();
    plotarea->Print(histoname);
	  //createGIF(TString(histos1->GetName()));
	  delete plotarea;
	  delete legend;
	}
	
	histos1 = signal[0]->After(histos1);
      
      }
            
      isdone = true;
      
      //////////////
      // clean data 
      
      std::vector<TObjArray *>::iterator tobji = signal.begin();
      while(tobji != signal.end()) {
	delete *tobji;
	++tobji;
      }
      
      tobji = background.begin();
      while(tobji != background.end()) {
	delete *tobji;
	++tobji;
      }
      
      processName_signal.clear();
      processName_background.clear();
      
    }
    
  }
  
  print_debug_message("combineSignalBack: Done.");
  
}


