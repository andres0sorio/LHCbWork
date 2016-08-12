#ifndef ONEANGLEFIT_H
#define ONEANGLEFIT_H

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
#include <TMinuit.h>

#include "Utilities.h"

#include "PDFs.h"



void fcn         (int &, double *, double &, double *, int );

int  getData     (const char *);

void fitData     (const char * , const char *, int);

///////
#endif
