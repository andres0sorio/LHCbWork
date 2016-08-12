// $Id: TestToyMC.h,v 1.1 2006/11/24 22:15:17 aosorio Exp $
#ifndef TESTMC_TESTTOYMC_H 
#define TESTMC_TESTTOYMC_H 1

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

/** @class TestToyMC TestToyMC.h testmc/TestToyMC.h
 *  
 *
 *  @author Andres Felipe OSORIO OLIVEROS
 *  @date   2006-11-23
 */
class TestToyMC {
public: 
  /// Standard constructor
  TestToyMC( ) { }; 
  
  TestToyMC( const char * );
  
  virtual ~TestToyMC( ); ///< Destructor
  

  TH1D *hh[10];
  
  void GetData(const char *);
  
  void DrawDistributions();
  

protected:
  
private:
  
  TFile *f1;
  
  TTree *T;
  
  TCanvas *acanv;
  
  std::vector<double> xi[6];

  int maxhistos;

  Double_t var[7];
  
  
};
#endif // TESTMC_TESTTOYMC_H
