// $Id: Projections.C,v 1.4 2007/04/12 11:11:20 aosorio Exp $
// Include files 



// local
#include "Projections.h"

//-----------------------------------------------------------------------------
// Implementation file for class : Projections
//
// 2007-03-18 : Andres Osorio
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
Projections::Projections( const char *distname, int bins, double min, double max  ) {

  mbins = bins;
  xmin  = min;
  xmax  = max;
  name  = std::string( distname );

  gDirectory->mkdir(distname)->cd();
  
  dist = new TH1D(distname ,"distribution"      , mbins, xmin , xmax);
  dist->Sumw2();
  
  gDirectory->cd("../");
  
  char cname[30];
  sprintf ( cname, "c1_%s", distname);
  acanv  = new TCanvas(cname, "Combined plot", 200, 10, 600, 420);
  acanv->Divide(1,1);
  
  colours[0]=1;
  colours[1]=2;
  colours[2]=4;
  linesty[0]=1;
  linesty[1]=4;
  linesty[2]=7;

  cstyle = new TStyle();
  setStyleOptions(cstyle);
  cstyle->SetStatBorderSize(0);
  cstyle->cd();
  
    
}

//=============================================================================
// Destructor
//=============================================================================

Projections::~Projections() {
  
  if (dist) delete dist;
  if (acanv) delete acanv;
  
  std::vector<TF1*>::iterator itr1;
  for( itr1=curves.begin();itr1!=curves.end();++itr1) delete (*itr1);
  
  if (cstyle) delete cstyle;
  
} 

//=============================================================================

void Projections::Initialise( double data , bool end)
{
  
  if ( !end ) {
    //std::cout << "filling histogram" << std::endl;
    dist->Fill( data );
  
  } else {
    //std::cout << "no more data " << dist << std::endl;
    sfactor1 = dist->Integral("width");
    //std::cout << "no more data" << std::endl;
    dist->Scale(1.0/sfactor1);
    //std::cout << "no more data" << std::endl;
    dist->SetMinimum(0.0);
    //std::cout << "no more data" << std::endl;
  }
  
}

void Projections::AddFunction( TF1 * pdf , int npar)
{
  
  curves.push_back( pdf );
  std::cout << "AddFucntion> TF1 function added \n";
  maxparameters = npar;
  
}

void Projections::AddFunction( const char *funcname, int npar, ptr2Function pdf )
{
  
  ff  = new TF1( funcname, pdf , xmin , xmax, npar );
  maxparameters = npar;
  curves.push_back( ff );
  std::cout << "AddFucntion> new TF1 function added \n";
  
}

void Projections::SetParameters( const char * filename ,  double qfac )
{
  readData(filename, v1, v2);
  
  std::cout << "SetParameters> total number of parameters read:"
            << v1.size() << " " << v2.size() << '\n';

  std::cout << "SetParameters> qfactor: " << qfac << '\n';
    
  int maxpars = v1.size();
  
  std::vector<TF1*>::iterator itr;
  for ( itr = curves.begin(); itr != curves.end(); ++itr ) {
    for (int k=0; k < maxpars; ++k) {
      (*itr)->FixParameter(k, v1[k]);
    }
  }
  
  for ( itr = curves.begin(); itr != curves.end(); ++itr ) {
    (*itr)->FixParameter(15, qfac);
    (*itr)->FixParameter(18, xmin);
    (*itr)->FixParameter(19, xmax);
  }
  
}

void Projections::SetTimeLimits( double t1, double t2 )
{
  std::vector<TF1*>::iterator itr;
  for ( itr = curves.begin(); itr != curves.end(); ++itr ) {
    (*itr)->FixParameter(18, t1);
    (*itr)->FixParameter(19, t2);
  }

  
}


void Projections::Draw( )
{
  
  std::cout << "PlotUtility> drawing projections ... " <<  curves.size() << std::endl;


  std::vector<TF1*>::iterator itr;

  acanv->cd(1);
  dist->SetTitle("");
  dist->Draw();
  
  int cc = 0;
  
   for ( itr = curves.begin(); itr != curves.end(); ++itr ) {

     std::cout << "Draw> drawing function: " << cc+1 << std::endl;
          
     ff = (TF1*)(*itr);
    
     ff->SetLineColor(colours[cc]);
     ff->SetLineStyle(linesty[cc]);
     ff->SetLineWidth(2);
     ff->SetNpx(150);
     ff->Draw("same");
     acanv->Update();
    
     ++cc;
    
   }
  
  std::cout << "PlotUtility> drawing projections ... Done! " << std::endl;
  
}

void Projections::DrawLog( )
{
  
  std::cout << "PlotUtility> drawing projections ... " <<  curves.size() << std::endl;


  std::vector<TF1*>::iterator itr;

  acanv->cd(1);
  dist->Scale(sfactor1);
  gPad->SetLogy(1);
  dist->SetTitle("");
  dist->Draw( );
  
  int cc = 0;
  
  for ( itr = curves.begin(); itr != curves.end(); ++itr ) {

    std::cout << "Draw> drawing function: " << cc+1 << std::endl;
    
    ff = (TF1*)(*itr);
    
    ff->FixParameter(20, sfactor1 );
    ff->SetLineColor(colours[cc]);
    ff->SetLineStyle(linesty[cc]);
    ff->SetLineWidth(2);
    ff->SetNpx(100);
    ff->Draw("same");
    acanv->Update();
    
    ++cc;
    
  }
  
  std::cout << "PlotUtility> drawing projections ... Done! " << std::endl;
  
}
