#include "root_utilities.h"

/******************************************/
/* File name, labels and string utilities */
/******************************************/

TString fileNameStub(const Char_t *fileName) {
  TString fileNameString = TString(fileName);
  int iStart = fileNameString.Last('/')+1;
  int iEnd   = fileNameString.Last('.');
  // use this is case there is a double extension like sim.gz
  TString temp = fileNameString(iStart,iEnd-iStart);
  iStart = 0;
  iEnd   = temp.Last('.');
  if(iEnd <= 0) iEnd = temp.Length();
  return temp(iStart,iEnd-iStart);  
}

TString fileNameProc(const Char_t *fileName) {
  TString fileNameString = TString(fileName);
  return fileNameString(10,6);
}

/**********************************************/
/* utilities for histograms, presentation etc */
/**********************************************/

void setAxeOptions(TAxis * ax) {
    
  ax->SetLabelSize(0.031);
  ax->SetLabelFont(62);
  ax->SetLabelOffset(0.007);
  
  ax->SetTitleSize(0.035);
  ax->SetTitleFont(62);
  ax->SetTitleOffset(1.4);
   
}

void set2dHistoOptions(TH2D * ahisto, Int_t acolor) {
  
  ahisto->SetFillColor(acolor);
  ahisto->SetMarkerColor(acolor);
  ahisto->SetLineColor(acolor);
  ahisto->SetMarkerSize(0.6);
  ahisto->SetMarkerStyle(6);
  
   
}

void setHistogramsOptions(TH1D *h) 
{
  TAxis::TAxis *axis1;
  TAxis::TAxis *axis2;
  
  h->SetMarkerStyle(6);
  h->SetMarkerSize(0.4);
  h->SetMarkerColor(2);
  h->SetLineColor(2);

  axis1 = h->GetXaxis();
  axis2 = h->GetYaxis();
  axis1->SetTitle(h->GetName());
  
  setAxeOptions(axis1);
  setAxeOptions(axis2);
}

void setHistogramsOptions(TH2D *hh) 
{
  
  TAxis::TAxis *axis1;
  TAxis::TAxis *axis2;

  hh->SetMarkerStyle(4);
  hh->SetMarkerSize(1.2);
  hh->SetMarkerColor(1);
  hh->SetLineColor(1);
  
  axis1 = hh->GetXaxis();
  axis2 = hh->GetYaxis();
  axis1->SetTitle(hh->GetName());
  
  setAxeOptions(axis1);
  setAxeOptions(axis2);
  
}

void setHistogramsOptionsCD(TH2D *hh, int colour) 
{
  TAxis::TAxis *axis1;
  TAxis::TAxis *axis2;
  
  hh->SetMarkerStyle(8);
  hh->SetMarkerSize(0.6);
  hh->SetMarkerColor(colour);
  hh->SetLineColor(2);
  
  axis1 = hh->GetXaxis();
  axis2 = hh->GetYaxis();
  axis1->SetTitle(hh->GetName());
  
  setAxeOptions(axis1);
  setAxeOptions(axis2);
  
}

void setHistogramsOptionsObjects(TH2D *hh, int colour, int marker) 
{
  
  TAxis::TAxis *axis1;
  TAxis::TAxis *axis2;
  
  hh->SetMarkerStyle(marker);
  hh->SetMarkerSize(0.4);
  hh->SetMarkerColor(colour);
  hh->SetLineColor(2);
  
  axis1 = hh->GetXaxis();
  axis2 = hh->GetYaxis();
  axis1->SetTitle(hh->GetName());
  
  setAxeOptions(axis1);
  setAxeOptions(axis2);
  
}

void setLegendOptions(TLegend * alegend) {
  
  alegend->SetFillColor(10);
  alegend->SetBorderSize(1);
  alegend->SetTextSize(0.03);
  
}

void setXaxisTitle(TH1D * ahisto, std::ifstream * optionFile) {
  
  std::string hname;
  std::string xtitle;
  std::string ytitle;
  std::string ahistoName;

  std::vector<std::string> names;
  std::vector<std::string> xnames;
  std::vector<std::string> ynames;
  std::vector<std::string>::iterator pos;
  
  optionFile->clear();
  optionFile->seekg(0, ios::beg);
  
  while (1) {
    
    if (optionFile->eof()) break;
    
    (*optionFile) >> hname >> xtitle >> ytitle;
    
    std::cout << hname << " " << xtitle << std::endl;
    
    names.push_back(hname);
    xnames.push_back(xtitle);
    ynames.push_back(ytitle);
    
  }
  
  
  ahistoName = std::string(ahisto->GetName());
  
  pos = std::find(names.begin(),
		  names.end(),
		  ahistoName);
  
  int loc(0);
  
  loc = pos - names.begin();
  
  if (pos != names.end() ) {
    
    TAxis *x = ahisto->GetXaxis();
    x->SetTitle( xnames[loc].c_str() );
    setAxeOptions(x);
    
    TAxis *y = ahisto->GetYaxis();
    y->SetTitle( ynames[loc].c_str() );
    setAxeOptions(y);
    
  }

}

void setXaxisTitle(TH2D * ahisto, std::ifstream * optionFile) {
  
  std::string hname;
  std::string xtitle;
  std::string ytitle;
  std::string ahistoName;

  std::vector<std::string> names;
  std::vector<std::string> xnames;
  std::vector<std::string> ynames;
  std::vector<std::string>::iterator pos;

  optionFile->clear();
  optionFile->seekg(0, ios::beg);
  
  while (1) {
    
    if (!optionFile->good()) break;
    
    (*optionFile) >> hname >> xtitle >> ytitle;
    
    std::cout << hname << " " << xtitle << ytitle << std::endl;

    names.push_back(hname);
    xnames.push_back(xtitle);
    ynames.push_back(ytitle);

  }

  ahistoName = std::string(ahisto->GetName());

  pos = std::find(names.begin(),
		  names.end(),
		  ahistoName);
  
  int loc(0);
  
  loc = pos - names.begin();
  
  if (pos != names.end() ) {

    TAxis *x = ahisto->GetXaxis();
    x->SetTitle( xnames[loc].c_str() );
    setAxeOptions(x);

    TAxis *y = ahisto->GetYaxis();
    y->SetTitle( ynames[loc].c_str() );
    setAxeOptions(y);
    
  }
  
}

void setStyleOptions(TStyle * style) {
  
  style->SetCanvasColor(36);
  style->SetTextFont(22);
  style->SetTitleBorderSize(0);
  //style->SetOptDate(22); 
  //style->GetAttDate()->SetTextFont(42);
  //style->GetAttDate()->SetTextSize(0.030);
  //style->GetAttDate()->SetTextColor();
  //style->SetOptStat(10);
  //style->SetOptStat(0);
  //style->SetStatFont(42);
  //style->SetStatFontSize(0.030);
  //style->SetStatColor(10);
  style->SetStatBorderSize(1);
  // style->SetStatW(0.25);
  // style->SetStatH(0.25);
  // style->SetStatX(0.93);
  //style->SetStatY(0.95);
  style->SetPaperSize(10.0,10.0);
  gErrorIgnoreLevel = 1;

  // "What're quantum mechanics?"
  // "I don't know. People who repair quantums I suppose."
  //--Rincewind, Terry Pratchett "Eric"

}

Int_t findBestQuadrant(TH1D *ahisto) {
  
  Int_t quadrant(0);
  Float_t qarea(0.0);
  Float_t qqarea[4];
  Float_t max_x(0.0), min_x(0.0), tot_x(0.0);
  Float_t sfactor(0.0);
  Int_t max_bin(0), sep_bin(0), k;
  
  Float_t value = 0;
  Float_t maxvalue =0;
  Float_t maximum = ahisto->GetMaximum();

  max_x=(ahisto->GetXaxis())->GetXmax();
  min_x=(ahisto->GetXaxis())->GetXmin();
  tot_x= max_x - min_x;
  
  qarea=(tot_x*maximum)/4;
  
  max_bin=ahisto->GetNbinsX();
  sep_bin=max_bin/4;

  sfactor=tot_x/max_bin;

  qqarea[0]=qarea-sfactor*(ahisto->Integral(1,sep_bin));
  qqarea[1]=qarea-sfactor*(ahisto->Integral(sep_bin+1,2*sep_bin));
  qqarea[2]=qarea-sfactor*(ahisto->Integral(2*sep_bin+1,3*sep_bin));
  qqarea[3]=qarea-sfactor*(ahisto->Integral(3*sep_bin+1,max_bin));
  
  for(k=0; k<4 ; k++) {
    value=qqarea[k];
    if(value > maxvalue) {
      maxvalue = value;
      quadrant=k;
    }
  }
  
  return quadrant;

}

 
void createGIF(TString name) {

  //If Imagemagick is not installed:
  // TString firstcommand = TString ("pstopnm -ppm -xborder 0 -yborder 0 -nocrop -portrait ");
  // TString secondcommand = TString ("ppmtogif ");  
  // TString execone = firstcommand + name + TString(".eps");
  // TString exectwo = secondcommand + name + TString(".eps001.ppm > ") + name + TString(".gif");
  
  TString firstcommand = TString ("convert -depth 16 -resize 768x512 ");
  TString execone = firstcommand + name + TString(".eps ") + name + TString(".tiff");
  
  if (gROOT->IsBatch())  {    
    
    gSystem->Exec(execone);
    //gSystem->Exec(exectwo);
    //gSystem->Exec("rm *.ppm");
 
  } 
  
  else print_message("createGIF> Cannot create gif file!");
  
}


TObjArray * readHistograms(TDirectory *dir) {
  
  TDirectory *current_sourcedir = gDirectory;
  TIter nextkey (current_sourcedir->GetListOfKeys());
  TKey *key;
  TObjArray *temp = new TObjArray();
  
  while((key = (TKey*)nextkey())) {
    
    current_sourcedir->cd();
    TObject *obj = key->ReadObj();
    
    if( obj->IsA()->InheritsFrom("TH1D")) temp->AddLast(obj);
    else if( obj->IsA()->InheritsFrom("TH2D")) temp->AddLast(obj);
    else print_message("Whoaaat's this object? ");
    
  }
  
  return temp;
  
}

TObjArray* sortH1(TObjArray *inHist, Int_t option)
{

  Int_t nh = inHist->GetEntriesFast();
  Float_t maximumY(0.0);
  Float_t maxA[nh];
  TObject *objects = inHist->First();
  Int_t i = 0;
  
  while(objects) {
    TH1D *h1 = (TH1D*)objects;
    maxA[i] = h1->Integral();
    ++i;
    objects = inHist->After(objects);   
  }
  
  //sort histograms. Store in a new TObjArray
  Int_t index[nh];
  TMath::Sort(nh,maxA,index);
  Int_t pos =0;
  TObjArray *temp = new TObjArray(nh);
  for( int k = 0; k < nh; k++) {
    //    printf("i=%d, index=%d, max=%f\n",k,index[k],maxA[index[k]]);
    pos = index[k];
    TH1D *h1 = (TH1D*)inHist->At(pos);
    maximumY = h1->GetMaximum();
    
    //We need to set maximum and minimum now so log scale plots look nicer
    //h1->SetMaximum(maximumY+maximumY*0.2);
    if(option == 1) {
      h1->SetMaximum(maximumY+maximumY*10.0);
      h1->SetMinimum(1.0);
    }
    else {
      h1->SetMaximum(maximumY+maximumY*0.3);
      h1->SetMinimum(0.0);
    }
    temp->Add(h1);
  }
  
  return temp;
  
}      


void returnDiffCrossSection(TH1D *baseHisto, TH1D *targetHisto, double totLum)
{
  
  ////////////////////////////////////////////
  int i(0);
  int n_bins(0);
  double scalefactor(0.0);
  double delta_m(0.0);
  double f_m(0.0);
  double value(0.0);
  double errorbin(0.0);
  
  scalefactor = 1/totLum;
  n_bins = baseHisto->GetNbinsX();
  delta_m = baseHisto->GetBinWidth(0);
  
  for( i = 1; i <= n_bins; i++) {
    
    f_m = baseHisto->GetBinContent(i);
    value = f_m * (scalefactor/delta_m);
    errorbin = baseHisto->GetBinError(i);
    errorbin = errorbin * (scalefactor/delta_m);
    targetHisto->SetBinContent(i,value);
    targetHisto->SetBinError(i,errorbin);
    
  }


}
