// $Id: RandomDists.h,v 1.10 2006/12/09 18:17:56 aosorio Exp $
#ifndef RANDOMDISTS_H 
#define RANDOMDISTS_H 1

// Include files

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <cmath>
#include <vector>


#include "Riostream.h"
#include <TROOT.h>
#include <TFile.h>
#include <TList.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TDirectory.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TString.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TMath.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <TPaveLabel.h>
#include <TKey.h>
#include <TLegend.h>
#include <TLine.h>
#include <TError.h>
#include <TFrame.h>
#include <TRandom.h>
#include <TRandom3.h>

#include "ResModels.h"
#include "PDFs.h"
#include "PDFsWacc.h"

//
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/RootFinder.h"
#include "Math/RootFinderAlgorithms.h"
#include "Math/WrappedFunction.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

/** @class RandomDists RandomDists.h RandomDists.h
 *  
 *
 *  Toy Monte Carlo for Bs -> Jpsi phi
 *  
 *  
 *
 *  @author Andres Osorio Oliveros
 *  @date   2006-07-04
 */

typedef double (*Ptr2pdf)  (double * , double *);
typedef double (*Ptr2Dpdf) (double   , void   *);

class RandomDists {

public: 
  /// Standard constructor
  RandomDists( ); 
  
  RandomDists(const char* , int); 
  
  virtual ~RandomDists( ); ///< Destructor
  
  TH1D *dist[10];
  std::vector<double> data[8];
  
  double xmin[4];
  double xmax[4];

  double par[20];
  
  int max_events;

  void initialise(double *, int );
  
  void SelectTagging();
  
  void RunExperiment();
  
  void WriteToTree();
  
  struct pdf_params * params;
  
protected:
  
private:
  
  TFile *outfile;
  TTree *t1;
  
  double pi;
  double var[9];
  double tempoval[4];
  bool   foundval[4];
  
  Ptr2pdf astPDF;
  double  fmax[2];
  double  tagfrac[2][2];
  double  pdfmax;
    
  bool generatingBs;

  TRandom3 *rdn[4]; // one random number for each variable
  TRandom3 *reg;    // rescaling number
  TRandom3 *tag;    // consider tagging
  
  void   findMaxima      ( );
  double findMaximum     ( int , Ptr2pdf , Ptr2Dpdf );
  void   updateMaxima    ( );
  void   setParameters   ( double * );
  
  ///////////////////////////////
  //Using gsl
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  
  
};

#endif // RANDOMDISTS_H
