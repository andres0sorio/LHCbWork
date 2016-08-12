#include "plot_util.hpp"

/**********************************************/
/* utilities for histograms, presentation etc */
/**********************************************/

void setAxeOptions(TAxis * ax) {
  
  ax->SetLabelSize(0.03);
  ax->SetLabelFont(42);
  ax->SetLabelOffset(0.006);
  
  ax->SetTitleSize(0.03);
  ax->SetTitleFont(42);
  ax->SetTitleOffset(1.5);
  
}

void setStyleOptions(TStyle * style) {
  
  style->SetCanvasColor(36);
  style->SetTextFont(22);
  style->SetOptDate(22); 
  style->GetAttDate()->SetTextFont(42);
  style->GetAttDate()->SetTextSize(0.030);
  style->GetAttDate()->SetTextColor();
  style->SetOptStat(10);
  style->SetStatFont(42);
  style->SetStatFontSize(0.030);
  style->SetStatColor(10);
  style->SetStatBorderSize(1);
  style->SetStatW(0.25);
  style->SetStatH(0.25);
  style->SetStatX(0.93);
  style->SetStatY(0.95);
  style->SetPaperSize(14.0,12.0);

  // "What're quantum mechanics?"
  // "I don't know. People who repair quantums I suppose."
  //--Rincewind, Terry Pratchett "Eric"

}

void createGIF(TString name) {

  TString firstcommand = TString ("pstopnm -ppm -xborder 0 -yborder 0 -xsize 450 -portrait ");
  TString secondcommand = TString ("ppmtogif ");
  
  TString execone = firstcommand + name + TString(".eps");
  TString exectwo = secondcommand + name + TString(".eps001.ppm > ") + name + TString(".gif");

  if (gROOT->IsBatch())  {
    
    gSystem->Exec(execone);
    gSystem->Exec(exectwo);
    gSystem->Exec("rm *.ppm");
  } 
  
  else std::cout << "createGIF> Cannot create gif file!" << std::endl;
  

}

