#include <TROOT.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TString.h>
#include <TFrame.h>
#include <TLatex.h>
#include <TLegend.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include <statsStruct.h>
#include <procStats.h>
#include <plotUtil.h>

////////////////////////////
// Main program starts here

int main(int iargv, const char **argv) {
  
  if(iargv <= 2) {
    std::cout << "usage : " << argv[0] 
	      << " inputFileName " 
	      << " option " << std::endl;
    return 1;
  }
  
  const char *fileName = argv[1];
  const char *option = argv[2];
  
  if (iargv == 3) {
    
    std::cout << "statsReader> reading file: " << fileName << std::endl;
    std::cout << "statsReader> andres@hep.man.ac.uk" << std::endl;
  }
  
  //load input file
  statsLoader *sdl = new statsLoader(fileName);
  
  //read first event on it
  AllStats *statsCut = sdl->next_file();
  
  int i = 0;
  
  TCanvas *c1 = new TCanvas("nM_WWvsCuts","",200,10,600,400);
  c1->SetFillColor(10);
  TCanvas *c2 = new TCanvas("nEntriesTotal","",200,10,600,400);  
  c2->SetFillColor(10);
  TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.05,0.95,0.95,10);
  TPad *pad2 = new TPad("pad2","This is pad1",0.05,0.05,0.95,0.95,10);
  
  c1->cd();
  pad1->Draw();
  c2->cd();
  pad2->Draw();
  
  TGraph *g1;
  TGraph *g2;
  TGraph *g3;
  TGraph *f1;
  TGraph *f2;
  TGraph *f3;
  
  std::vector<double> nEvt1;
  std::vector<double> nEvt2;
  std::vector<double> nEvt3;
  std::vector<int> nEntries1;
  std::vector<int> nEntries2;
  std::vector<int> nEntries3;

  int ymaximum = 0;
  int zmaximum = 0;
  
  bool isSignal = false;
  bool isBackground = false;

  std::string cutnames[10];
  setCutNames(cutnames);
    
  // loop over all cuts
  
  while(statsCut) {
    
    std::cout << i+1 << "::" << std::endl;
    
    procStats *adir = new procStats(statsCut);
    
    //   adir->printStats();

    adir->loadData();
    
    adir->getResults();
    
    adir->produceTable();

    // to plot results
    if(adir->isSignal) {
      nEvt1.push_back(adir->nEventsWWdetAv);
      nEvt2.push_back(adir->nEventsZZdetAv);
      nEvt3.push_back(adir->nEventsBKdetAv);
      nEntries1.push_back(adir->wwEntriesTot);
      nEntries2.push_back(adir->zzEntriesTot);
      nEntries3.push_back(adir->bkEntriesTot);
      isSignal = true;
    }
    else if(adir->isBackground) {
      nEvt1.push_back(adir->nEventsBKdetAv);
      isBackground = true;
    }
    
    delete statsCut;
    
    delete adir;
    
    statsCut = sdl->next_file();
    
    i++;    
    
  }

  
  //plot entries on a graph and print it
  
  TString epsFile1 = TString("StudyOnCuts") 
    + TString(option) + TString("1.eps");
  TString epsFile2 = TString("StudyOnCuts") 
    + TString(option) + TString("2.eps");
  
  int nCuts = i;
  double x[20];
  TPaveLabel *boxes[10];

  if(isSignal) {

    double y1[20], y2[20], y3[20];
    double z1[20], z2[20], z3[20];
    for(int k = 0; k < nCuts; k++) {
      x[k] = k+1.0;
      y1[k] = nEvt1[k];
      y2[k] = nEvt2[k];
      y3[k] = nEvt3[k];
      z1[k] = nEntries1[k];
      z2[k] = nEntries2[k];
      z3[k] = nEntries3[k];
    }
    
    ymaximum = findMaximum(0,nEvt1);
    ymaximum = findMaximum(ymaximum,nEvt2);
    ymaximum = findMaximum(ymaximum,nEvt3);
    ymaximum = ymaximum + int (std::ceil(ymaximum*0.30));
    zmaximum = findMaximum(0,nEntries1);
    zmaximum = findMaximum(zmaximum,nEntries2);
    zmaximum = findMaximum(zmaximum,nEntries3);
    zmaximum = zmaximum + int (std::ceil(zmaximum*0.30));

    pad1->cd();
    pad1->DrawFrame(0.0,0.0,11.0,ymaximum);
    setFrameOptions(pad1, "Cuts", "Events_{process}");

    g1 = new TGraph(nCuts,x,y1);
    g1->SetLineColor(4);
    g1->SetMaximum(ymaximum);
    setGraphOptions(g1,"Cuts", "Events_{process}");
    
    g2 = new TGraph(nCuts,x,y2);
    g2->SetLineColor(5);
    setGraphOptions(g2,"","");
    
    g3 = new TGraph(nCuts,x,y3);
    g3->SetLineColor(6);
    setGraphOptions(g3,"","");

    TLegend *legend;
    legend = new TLegend(0.3,0.77,0.3+0.18,0.88);
    TString legendItem = TString("WW#nu#nu");
    legend->AddEntry(g1, legendItem,"l");
    legendItem = TString("ZZ#nu#nu");
    legend->AddEntry(g2, legendItem,"l");
    legendItem = TString("qqqq#nu#nu (background)");
    legend->AddEntry(g3, legendItem,"l");
    legend->SetFillColor(10); 
    legend->SetBorderSize(0); 
    legend->SetTextSize(0.04);

    pad1->cd();
    g1->Draw("CP");
    g2->Draw("CP");
    g3->Draw("CP");
    createLabels(pad1, boxes, cutnames, 10);
    legend->Draw();
    
    c1->Update();
    
    ///////////////////////////////

    f1 = new TGraph(nCuts,x,z1);
    f1->SetLineColor(4);
    f1->SetMaximum(zmaximum);
    setGraphOptions(f1,"Cuts","Entries_{process}");
    
    f2 = new TGraph(nCuts,x,z2);
    f2->SetLineColor(5);
    setGraphOptions(f2,"","");
    
    f3 = new TGraph(nCuts,x,z3);
    f3->SetLineColor(6);
    setGraphOptions(f3,"","");

    pad2->cd();
    f1->Draw("CP");
    f2->Draw("CP");
    f3->Draw("CP");

    c1->Print(epsFile1);
    c2->Print(epsFile2);

  }
  
  else if(isBackground) {
    
    double y1[20];
    
    for(int k = 0; k < nCuts; k++) {
      x[k] = k+1.0;
      y1[k] = nEvt1[k];
    }
    
    ymaximum = findMaximum(0,nEvt1);
    //ymaximum = ymaximum + int (std::ceil(ymaximum*0.30));
    ymaximum = ymaximum + ymaximum*10;

    /////////////////////////////////////
      //
      pad1->cd();
      pad1->SetLogy();
      pad1->DrawFrame(0.001,1,nCuts+1.0,ymaximum);
      
      g1 = new TGraph(nCuts,x,y1);
      g1->SetLineColor(4);

    setGraphOptions(g1,"Cuts", "Events_{process}");
    setFrameOptions(pad1, "Cuts", "Events_{process}");

    g1->Draw("CP");
    createLabels(pad1, boxes, cutnames, 10);

    TLegend *legend;
    legend = new TLegend(0.3,0.77,0.3+0.18,0.88);
    TString legendItem = TString("Total 4-2 fermions background");
    legend->AddEntry(g1, legendItem,"l");
    legend->SetFillColor(10); 
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);

    legend->Draw();
    pad1->Update();
    c1->Update();
    c1->Print(epsFile1);

    
  }
  
  
  std::cout << "statsReader> Total cuts: " << i << std::endl;
  
  delete sdl;
   
  return 0;
  
}


