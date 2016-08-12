#include "PlotUtility.h"

const int MAXPARS=17;

void PlotResults(const char *infile,const char *results )
{
  
  TFile *f1 = new  TFile (infile, "READ");
  
  if( f1 ) {
    std::cout << "getData> File: " << std::string(infile) 
	      << " opened." << std::endl;
  }
  else {
    std::cout << "getData> Error: File: " << std::string(infile) 
	      << " not accessible" << std::endl;
    exit(1);
  }
  
  f1->cd();
  TTree *T = (TTree*)gROOT->FindObject("data");
  
  Double_t var[6];
  
  T->SetBranchAddress("ptime"   ,&var[0]);
  T->SetBranchAddress("theta"   ,&var[1]);
  T->SetBranchAddress("psi"     ,&var[2]);
  T->SetBranchAddress("phi"     ,&var[3]);
  T->SetBranchAddress("ptwres"  ,&var[4]);
  T->SetBranchAddress("thwres"  ,&var[5]);
  
  int n_events = T->GetEntriesFast();
  
  TCanvas *cv1 = new TCanvas("cv1","Plots",200,10,700,700);
  cv1->Divide(2,2);
  TCanvas *cv2 = new TCanvas("cv2","Plots",250,10,700,700);
  cv2->Divide(2,2);
  
  TH1D *dist[10];
  dist[0] = new TH1D("h1","distribution1"      ,60,   0.0, 10.0);
  dist[1] = new TH1D("h2","distribution1"      ,60,  -0.1, 3.16);
  dist[2] = new TH1D("h3","distribution1"      ,60,  -0.1, 3.16);
  dist[3] = new TH1D("h4","distribution1"      ,60,  -0.1, 6.20);
  dist[4] = new TH1D("h5","distribution1"      ,50,  -1.0, 1.00);
  dist[0]->Sumw2();
  dist[1]->Sumw2();
  dist[2]->Sumw2();
  dist[3]->Sumw2();
  dist[4]->Sumw2();
  //TH2D *dist2d[5];
  //dist2d[0] = new TH2D("h6","distribution1"      ,50,-1.0, 1.0,50,-1.0, 1.0);
  
  for(int i=0; i <  n_events; ++i){
    T->GetEntry(i);
    dist[0]->Fill(var[4]);
    dist[1]->Fill(var[1]);
    dist[2]->Fill(var[2]);
    dist[3]->Fill(var[3]);
    dist[4]->Fill(TMath::Sin(var[1])*TMath::Sin(var[3]));
    //dist2d[0]->Fill(TMath::Cos(var[1]),TMath::Sin(var[2]));
  }
  
  double sfactor1 = dist[0]->Integral("width");
  double sfactor2 = dist[1]->Integral("width");
  double sfactor3 = dist[2]->Integral("width");
  double sfactor4 = dist[3]->Integral("width");
  double sfactor5 = dist[4]->Integral("width");
  
  //std::cout << sfactor5 << std::endl;
  
  double par[MAXPARS];
  double errpar[MAXPARS];
  std::vector<double> v1;
  std::vector<double> v2;

  readData(results, v1, v2);

  //needed to insert two extra parameters :(
  //not clear and clear
  //could change this in the future
  
  if( (int)v1.size() < (MAXPARS-2) ) {std::cout << "stepTwo> Error readingData" << '\n'; exit(1);}
  
  for(int i =0; i < 9; ++i) {
    par[i] = v1[i];
    errpar[i] = (int) v2[i];
  }
  par[9]  = 0.0;
  par[10] = 1.0;
  for(int i = 9; i < (MAXPARS-2); ++i) {
    par[i+2] = v1[i];
    errpar[i+2] = (int) v2[i];
  }

  for(int i = 0; i < MAXPARS; ++i) {
    std::cout << i << '\t' << par[i] << std::endl;
  }

  ///////////////////////////////
  
  double x[1];
  double area1 = integrate( properTimePDF, x , par,-2.0,10.0);
  double area2 = integrate( thetaPDF, x , par, 0 , 3.1415);
  double area3 = integrate( psiPDF, x , par, 0 , 3.1415);
  double area4 = integrate( phiPDF, x , par, 0 , 6.30);
  double area5 = integrate( CosthetaTrTotal , x , par, -1.0, 1.0);
  
  double pi     =TMath::Pi();
  
  TF1 *g[15];
  g[0]  = new TF1("pdf1", properTimePDF  ,      0.0,   10.0      ,17); //Proper time projection
  g[1]  = new TF1("pdf2", thetaPDF        ,     0.0,   pi        ,11); //Theta of Jpsi
  g[2]  = new TF1("pdf3", psiPDF          ,     0.0,   pi        ,11); //Theta of phi 
  g[3]  = new TF1("pdf4", phiPDF          ,     0.0,   6.20      ,11); //phi (phi1+phi2)
  g[4]  = new TF1("proj2Even", thetaEven   ,     0.0,   pi        ,11); //Theta of Jpsi
  g[5]  = new TF1("proj2Odd", thetaOdd     ,     0.0,   pi        ,11); //Theta of Jpsi
  g[6]  = new TF1("cosprj2Even", CosthetaTrEven     ,  -1.0,   1.0      ,11); //Theta of Jpsi
  g[7]  = new TF1("cosproj2Odd", CosthetaTrOdd     ,  -1.0,   1.0      ,11); //Theta of Jpsi
  g[8]  = new TF1("cosproj2Tot", CosthetaTrTotal   ,  -1.0,   1.0      ,11); //Theta of Jpsi
  g[9]  = new TF1("ptimeOdd"   , properTimePDFOdd ,  0.0,   10.0      ,11); //Theta of Jpsi
  g[10]  = new TF1("ptimeEven" , properTimePDFEven ,   0.0,   10.0      ,11); //Theta of Jpsi
  g[11]  = new TF1("psiEven"   , psiPDFEven         ,     0.0,   pi        ,11); //Theta of Jpsi
  g[12]  = new TF1("psiOdd"    , psiPDFOdd          ,     0.0,   pi        ,11); //Theta of Jpsi
  g[13]  = new TF1("phiEven"   , phiPDFEven         ,     0.0,   6.20      ,11); //phi (phi1+phi2)
  g[14]  = new TF1("phiOdd"    , phiPDFOdd          ,     0.0,   6.20      ,11); //phi (phi1+phi2)



  for( int k = 0; k<11; ++k) {
    g[0]->FixParameter(k, par[k]);      
    g[1]->FixParameter(k, par[k]);
    g[2]->FixParameter(k, par[k]);
    g[3]->FixParameter(k, par[k]);
    g[4]->FixParameter(k, par[k]);
    g[5]->FixParameter(k, par[k]);
    g[6]->FixParameter(k, par[k]);
    g[7]->FixParameter(k, par[k]);
    g[8]->FixParameter(k, par[k]);
    g[9]->FixParameter(k, par[k]);
    g[10]->FixParameter(k, par[k]);
    g[11]->FixParameter(k, par[k]);
    g[12]->FixParameter(k, par[k]);
    g[13]->FixParameter(k, par[k]);
    g[14]->FixParameter(k, par[k]);
  }
  
  g[0]->FixParameter(11, par[11]);
  g[0]->FixParameter(12, par[12]);
  g[0]->FixParameter(13, par[13]);
  g[0]->FixParameter(14, par[14]);
  g[0]->FixParameter(16, par[15]);
  g[0]->FixParameter(16, par[16]);
  g[0]->FixParameter(10, sfactor1 / area1 );
  g[0]->FixParameter(10, sfactor1 / area1 );

  g[9]->FixParameter(10, sfactor1 / area1 );
  g[10]->FixParameter(10, sfactor1 / area1 );
  
  g[1]->FixParameter(9,  0.0);
  g[1]->FixParameter(10, sfactor2 / area2 );
  g[4]->FixParameter(9,  0.0);
  g[4]->FixParameter(10, sfactor2 / area2 );
  g[5]->FixParameter(9,  0.0);
  g[5]->FixParameter(10, sfactor2 / area2 );

  g[6]->FixParameter(9,  0.0);
  g[6]->FixParameter(10, sfactor5 / area5 );

  g[7]->FixParameter(9,  0.0);
  g[7]->FixParameter(10, sfactor5 / area5 );

  g[8]->FixParameter(9,  0.0);
  g[8]->FixParameter(10, sfactor5 / area5 );


  g[2]->FixParameter(9,  0.0);
  g[2]->FixParameter(10, sfactor3 / area3 );
  g[11]->FixParameter(9,  0.0);
  g[11]->FixParameter(10, sfactor3 / area3 );
  g[12]->FixParameter(9,  0.0);
  g[12]->FixParameter(10, sfactor3 / area3 );
  
  
  g[3]->FixParameter(9,  0.0);
  g[3]->FixParameter(10, sfactor4 / area4 );
  g[13]->FixParameter(9,  0.0);
  g[13]->FixParameter(10, sfactor4 / area4 );
  g[14]->FixParameter(9,  0.0);
  g[14]->FixParameter(10, sfactor4 / area4 );

  //  
  std::cout << "PlotUtility> drawing projections ..." << std::endl;
  cv1->cd(1);
  g[0]->SetLineColor(2);
  g[0]->SetLineWidth(1);
  dist[0]->Draw();
  g[0]->Draw("same");
  g[10]->SetLineColor(4); 
  g[10]->SetLineWidth(1);
  g[10]->Draw("same");
  g[9]->SetLineColor(4);
  g[9]->SetLineWidth(1);
  g[9]->SetLineStyle(2);
  g[9]->Draw("same");

  cv1->cd(2);
  g[1]->SetLineColor(2);
  g[1]->SetLineWidth(1);
  dist[1]->Draw();
  g[1]->Draw("same");
  g[4]->SetLineColor(4); 
  g[4]->SetLineWidth(1);
  g[4]->Draw("same");    //Even
  g[5]->SetLineColor(4);
  g[5]->SetLineWidth(1);
  g[5]->SetLineStyle(2);
  g[5]->Draw("same");    //Odd

  cv1->cd(3);
  g[2]->SetLineColor(2);
  g[2]->SetLineWidth(1);
  dist[2]->Draw();
  g[2]->Draw("same");
  g[11]->SetLineColor(4); 
  g[11]->SetLineWidth(1);
  g[11]->Draw("same");    //Even
  g[12]->SetLineColor(4);
  g[12]->SetLineWidth(1);
  g[12]->SetLineStyle(2);
  g[12]->Draw("same");    //Odd

  
  cv1->cd(4);
  g[3]->SetLineColor(2); 
  g[3]->SetLineWidth(1);
  dist[3]->Draw();
  g[3]->Draw("same");
  g[13]->SetLineColor(4); 
  g[13]->SetLineWidth(1);
  g[13]->Draw("same");    //Even
  g[14]->SetLineColor(4);
  g[14]->SetLineWidth(1);
  g[14]->SetLineStyle(2);
  g[14]->Draw("same");    //Odd


  cv1->cd();
  
  cv2->cd(1);
  dist[4]->SetMinimum(0.);
  dist[4]->Draw();
  g[8]->SetLineColor(2); //Total
  g[8]->SetLineWidth(1);
  g[8]->Draw("same");
  g[6]->SetLineColor(4); //Even
  g[6]->SetLineWidth(1);
  g[6]->Draw("same");
  g[7]->SetLineColor(4); //Odd
  g[7]->SetLineWidth(1);
  g[7]->SetLineStyle(2);
  g[7]->Draw("same");
  
  cv2->cd(2);
  //dist2d[0]->SetOption("Box");
  //dist2d[0]->Draw();
  g[7]->SetLineColor(4); //Odd
  g[7]->SetLineWidth(1);
  g[7]->SetLineStyle(2);
  g[7]->Draw();

  //cv2->cd(3);
  //dist[4]->SetMinimum(0.);
  
  //dist[4]->Draw();
  //g[6]->Draw("same");
  //g[7]->Draw("same");
  // g[8]->SetLineColor(2); //Total
  //g[8]->SetLineWidth(1);
  //g[8]->Draw("same");
  
  cv2->cd();
  
}


void PlotUtility(double xmin, double xmax)
{
  //
  gROOT->Reset();

  TCanvas *c1 = new TCanvas("c1","Weak phase analysis",100,10,900,600);
  
  c1->SetFillColor(21);
  c1->SetGrid();
  c1->Divide(5,2);
  
  const int n = 20;
  double x[n];
  double y[n];
  
  double dx= (xmax-xmin)/ (double)n;

  const int ng = 10;
  double tm= 10.0 / (double)ng;
  
  TGraph *go[ng];
  
  double xx[4];
  double par[9];
  
  xx[0]=0.0;
  xx[1]=0.3;
  xx[2]=0.3;
  xx[3]=0.6;

  par[0]=0.64;         //Gamma
  par[1]=0.10;         //DGrate
  par[2]=0.14;         //Rt
  par[3]=0.64;         //Rp
  par[4]=0.000;        //tphase1
  par[5]=TMath::Pi();  //tphase2
  par[6]=-0.04;        //w phase
  par[7]=20.00;        //dms
  par[8]=0.5;          //omega

  int j = 0;
  
  while( j < ng) {
    
    double time = 0.5 + j*tm;
    
    xx[0] = time; //fix the time
    
    for(int i = 0;i < n; ++i) {
      
      double p = xmin + i*dx;
      par[6] = p; //wphase is the variable
      y[i]   = jpsiphipdf(xx,par);
      x[i]   = p;
      
    }
   
    go[j] = new TGraph(n,x,y);
    go[j]->SetLineColor(2);
    go[j]->SetLineWidth(1);
    go[j]->SetMarkerColor(4);
    go[j]->SetMarkerStyle(7);
    
    ++j;
    
  }
  

  for( int i = 0; i < ng; ++i) {
    c1->cd(i+1);
    gPad->SetFillColor(10);
    setGraphOptions(*go[i]);
    go[i]->Draw("ACP");
  }

  c1->cd();
  

}


double projectTo1D(double *x, double *par)
{

  //double Gamma    = par[0];
  //double DGrate   = par[1];
  //double Rt       = par[2];
  //double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  //double wphase   = par[6];
  //double dms      = par[7];
  //double omega    = par[8];

  double low[3]; //limits of integration
  double upp[3];
  
  //use par[8] to control the variable to which the fucntion is projected
  low[0] = par[10];
  upp[0] = par[11];
  low[1] = par[12];
  upp[1] = par[13];
  low[2] = par[14];
  upp[2] = par[15];
  
  double result = projectTo1D (jpsiphipdfVar, x, par, low, upp, 50000, 4);
  
  return result;
  
}


void PlotPDFs()
{

  double par[20];
  
  //double Gamma    = par[0]; //BarGamma
  //double DGrate   = par[1]; //DeltaGamma / BarGamma
  //double Rt       = par[2];
  //double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  //double wphase   = par[6];
  //double dms      = par[7];
  
  par[0] = 0.70;
  par[1] = 0.15;
  par[2] = 0.20;
  par[3] = 0.60;
  par[4] = 0.0;
  par[5] = TMath::Pi();
  par[6] = -0.04;
  par[7] = 20.0;
  par[8] = 0.30;
  
  TCanvas *cv1 = new TCanvas("cv1","Plots",200,10,700,700);
  cv1->Divide(2,2);
  
  //Do some plots!
  std::cout << "PlotUtility> plotting the projection of the llh on the data" << std::endl;
  
  double pi     =TMath::Pi();
  
  TF1 *g[8];
  g[0]  = new TF1("pdf1", projectTo1D ,  0.0,   5.0       ,16); //Proper time projection
  g[1]  = new TF1("pdf2", projectTo1D ,  0.0,   pi        ,16); //Theta of Jpsi
  g[2]  = new TF1("pdf3", projectTo1D ,  0.0,   pi        ,16); //Theta of phi 
  g[3]  = new TF1("pdf4", projectTo1D ,  -pi,   pi        ,16); //phi (phi1+phi2)
  g[4]  = new TF1("pdf1b", properTimePDF ,  0.0,   10.0       ,11); //Proper time projection
  g[5]  = new TF1("pdf2b", thetaPDF      ,  0.0,   pi        ,11); //Theta of Jpsi
  g[6]  = new TF1("pdf3b", psiPDF        ,  0.0,   pi        ,11); //Theta of phi 
  g[7]  = new TF1("pdf4b", phiPDF        ,  -pi,   pi        ,11); //phi (phi1+phi2)
    
  for( int k = 0; k<9; ++k) {
    g[0]->FixParameter(k, par[k]);      
    g[1]->FixParameter(k, par[k]);
    g[2]->FixParameter(k, par[k]);
    g[3]->FixParameter(k, par[k]);
    g[4]->FixParameter(k, par[k]);      
    g[5]->FixParameter(k, par[k]);
    g[6]->FixParameter(k, par[k]);
    g[7]->FixParameter(k, par[k]);
    
  }
  
  /// time
  g[0]->FixParameter(9 , 1.0);
  g[0]->FixParameter(10, 0.0);//theta limits
  g[0]->FixParameter(11, pi);
  g[0]->FixParameter(12, 0.0);//psi limits
  g[0]->FixParameter(13, pi);
  g[0]->FixParameter(14, -pi);//phi limits
  g[0]->FixParameter(15,  pi);
  g[4]->FixParameter(9,  0.0);
  g[4]->FixParameter(10, 1.0);
  
  
  //// theta 1
  g[1]->FixParameter(9,   2.0);
  g[1]->FixParameter(10,  0.0);//time limits
  g[1]->FixParameter(11, 20.0);
  g[1]->FixParameter(12,  0.0);//psi limits
  g[1]->FixParameter(13,   pi);
  g[1]->FixParameter(14,  -pi);// phi limits
  g[1]->FixParameter(15,   pi);
  g[5]->FixParameter(9,  0.0);
  g[5]->FixParameter(10, 1.0);
  
  //// theta 2
  g[2]->FixParameter(9,   3.0);
  g[2]->FixParameter(10,  0.0);//time limits
  g[2]->FixParameter(11, 20.0);
  g[2]->FixParameter(12,  0.0);//psi
  g[2]->FixParameter(13,   pi);
  g[2]->FixParameter(14,  -pi);//phi
  g[2]->FixParameter(15,   pi);
  g[6]->FixParameter(9,  0.0);
  g[6]->FixParameter(10, 1.0);
  
  //// phi
  g[3]->FixParameter(9,   4.0);
  g[3]->FixParameter(10,  0.0);//time limits
  g[3]->FixParameter(11, 20.0);
  g[3]->FixParameter(12,  0.0);//theta
  g[3]->FixParameter(13,   pi);
  g[3]->FixParameter(14,  0.0);//psi
  g[3]->FixParameter(15,   pi);
  g[7]->FixParameter(9,   0.0);
  g[7]->FixParameter(10,  1.0);
  
  //  
  std::cout << "PlotUtility> drawing projections ..." << std::endl;
  cv1->cd(1);
  g[0]->SetLineColor(2);
  g[0]->SetLineWidth(1);
  g[4]->SetLineColor(1);
  g[4]->SetLineWidth(1);
  g[0]->Draw();
  g[4]->Draw("same");
  
  cv1->cd(2);
  g[1]->SetLineColor(2);
  g[1]->SetLineWidth(1);
  g[5]->SetLineColor(1);
  g[5]->SetLineWidth(1);
  g[5]->Draw();
  g[1]->Draw("same");
  
  cv1->cd(3);
  g[2]->SetLineColor(2);
  g[2]->SetLineWidth(1);
  g[6]->SetLineColor(1);
  g[6]->SetLineWidth(1);
  g[6]->Draw();
  g[2]->Draw("same");

  cv1->cd(4);
  g[3]->SetLineColor(2);
  g[3]->SetLineWidth(1);
  g[7]->SetLineColor(1);
  g[7]->SetLineWidth(1);
  g[7]->Draw();
  g[3]->Draw("same");
  //g[3]->Draw();
  cv1->cd();
  
  //cv1->Print("projections.eps");
  
  
  


}

double ptConvolMod1(double *xx, double *par)
{

  double convpars[1];
  convpars[0] = 0.3;
  
  double tiwres = convolve( properTimePDF, withGaussian, par, xx, convpars);
  
  return tiwres;
  
}

double ptConvolMod2(double *xx, double *par)
{

  double convpars[6];
  convpars[0] = 0.500;
  convpars[1] = 0.80;     //f1
  convpars[2] = 0.0022;   //mu1
  convpars[3] = 0.0584;   //s1
  convpars[4] = -0.1914;  //mu2
  convpars[5] = 0.471;    //s2
  

  double tiwres = convolve( properTimePDF, with2Gaussians, par, xx, convpars);
  
  return tiwres;
  
}

double logptConvolMod1(double *xx, double *par)
{

  double convpars[1];
  convpars[0] = 0.3;
  
  double tiwres = convolve( properTimePDF, withGaussian, par, xx, convpars);

  if( tiwres == 0.0 ) return 1.0;

  return TMath::Log10(tiwres);
  
}

double logptConvolMod2(double *xx, double *par)
{

  double convpars[6];
  convpars[0] = 0.500;
  convpars[1] = 0.80;     //f1
  convpars[2] = 0.0022;   //mu1
  convpars[3] = 0.0584;   //s1
  convpars[4] = -0.1914;  //mu2
  convpars[5] = 0.471;    //s2
  

  double tiwres = convolve( properTimePDF, with2Gaussians, par, xx, convpars);
  
  if ( tiwres == 0.0 ) return 1.0;

  return TMath::Log10(tiwres);
  
}



void checkResModel()
{
  
  double par[20];
  
  //double Gamma    = par[0]; //BarGamma
  //double DGrate   = par[1]; //DeltaGamma / BarGamma
  //double Rt       = par[2];
  //double Rp       = par[3];
  //double tphase1  = par[4];
  //double tphase2  = par[5];
  //double wphase   = par[6];
  //double dms      = par[7];

  par[0] = 0.66;
  par[1] = 0.22;
  par[2] = 0.14;
  par[3] = 0.60;
  par[4] = 0.0;
  par[5] = TMath::Pi();
  par[6] = 0.04;
  par[7] = 20.0;
  par[8] = 0.03;
  par[9] = 0.3;
  par[10]= 1.0;

  TCanvas *cv1 = new TCanvas("cv1","Plots",200,10,700,700);
  cv1->Divide(2,2);
  
  //Do some plots!
  std::cout << "PlotUtility> looking at resolution model" << std::endl;
  
  //double pi     =TMath::Pi();
  //double piover2=TMath::Pi()/2.0;

  TF1 *g[8];
  g[0]  = new TF1("pdf1b", ptConvolMod1,  -2.0,   10                 ,11); 
  g[1]  = new TF1("pdf2b", ptConvolMod2,  -2.0,   5.0                ,11); 
  g[2]  = new TF1("logpdf1b", logptConvolMod1,  -1.0,   10.0         ,11); 
  g[3]  = new TF1("logpdf2b", logptConvolMod2,  -1.0,   10.0         ,11); 
  
  for( int k = 0; k<11; ++k) {
    g[0]->FixParameter(k, par[k]);
    g[1]->FixParameter(k, par[k]);
    g[2]->FixParameter(k, par[k]);
    g[3]->FixParameter(k, par[k]);
  }
  
  //  
  std::cout << "PlotUtility> drawing projections ..." << std::endl;
  cv1->cd(1);
  g[0]->SetLineColor(2);
  g[0]->SetLineWidth(1);
  g[0]->Draw();

  cv1->cd(2);
  g[1]->SetLineColor(2);
  g[1]->SetLineWidth(1);
  g[1]->Draw();

  cv1->cd(3);
  g[2]->SetLineColor(2);
  g[2]->SetLineWidth(1);
  g[2]->Draw();

  cv1->cd(4);
  g[3]->SetLineColor(2);
  g[3]->SetLineWidth(1);
  g[3]->Draw();
  
  cv1->cd();
  
  

}


// void testingToy( int maxev )
// {

//   std::vector<double> xi[6];

//   gROOT->Reset();
  
//   TCanvas *cv1 = new TCanvas("cv1","Plots",200,10,600,600);
//   cv1->Divide(3,2);
  
//   TH1D *cosine  = new TH1D("cosine","Angular distribution",100,-1.1,1.1);
//   cosine->Sumw2();
  
//   double par[10];
//   double xmin[4];
//   double xmax[4];
  
//   std::ofstream *os = new std::ofstream("data.out",ofstream::out);
  
//   par[0]=0.72;         //Gamma
//   par[1]=0.17;         //DGrate
//   par[2]=0.14;         //Rt
//   par[3]=0.64;         //Rp
//   par[4]=0.000;        //tphase1
//   par[5]=TMath::Pi();  //tphase2
//   par[6]=-0.04;        //w phase
//   par[7]=20.00;        //dms
    
//   RandomDists *exp = new RandomDists("testingToyData.root", maxev);
//   exp->RunExperiment(par,8);
  
//   for(int i = 0; i<4; ++i) {
//     xi[i]   = exp->data[i];
//     xmin[i] = exp->xmin[i];
//     xmax[i] = exp->xmax[i];
//   }
  
//   for(int i = 0; i < maxev; ++i) {

//     double xx[4];

//     xx[0]=xi[0][i];
//     xx[1]=xi[1][i];
//     xx[2]=xi[2][i];
//     xx[3]=xi[3][i];
    
//     (*os) << jpsiphipdfA(xx,par) << '\t' 
// 	  << jpsiphipdfB(xx,par) << '\t'
// 	  << jpsiphipdf(xx,par)  << '\t'
// 	  << std::endl;

//     cosine->Fill(TMath::Cos(xi[1][i]));
    
//   }
  

//   cv1->cd();
//   cv1->cd(1);
//   exp->dist[0]->Draw();
//   cv1->cd(2);
//   exp->dist[1]->Draw();
//   cv1->cd(3);
//   exp->dist[2]->Draw();
//   cv1->cd(4);
//   exp->dist[3]->Draw();
//   cv1->cd(5);
//   cosine->Draw();
  
//   TF1 *g1    = new TF1("pdf",decayratepdfw,-0.1,0.1,7);
//   g1->FixParameter(0,par[0]);
//   g1->FixParameter(1,par[1]);
//   g1->FixParameter(2,par[2]);
//   g1->FixParameter(3,par[3]);
//   g1->FixParameter(4,par[4]);
//   g1->FixParameter(5,par[5]);
//   g1->FixParameter(6,par[7]);
  
//   TCanvas *cv2 = new TCanvas("cv2","Plots",200,10,600,600);
//   cv2->cd();
//   g1->Draw();
  
// }
