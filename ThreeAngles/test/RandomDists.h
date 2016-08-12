// $Id: RandomDists.h,v 1.3 2006/11/24 22:14:57 aosorio Exp $
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


/** @class RandomDists RandomDists.h RandomDists.h
 *  
 *
 *  This class is ready for production and it is intended to be
 *  used as part of a Toy Monte Carlo study
 *
 *  @author Andres Osorio Oliveros
 *  @date   2006-07-04
 */

typedef double (*Ptr2pdf) (double *, double *);

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
  
  void RunExperiment(double *, int );
  
  void CopyDataTo(int , double *);
  
  struct pdf_params * params;

protected:
  
  double ptMax;
  double thetaMax;
  double psiMax;
  double phiMax;
  
private:

  TFile *outfile;
  
  TTree *t1;
  
  //....
  Double_t var[9];
  
  double pi;
  
  void findMaxima( double * );
  
  void setParameters( double * );
  
  
};
#endif // RANDOMDISTS_H
