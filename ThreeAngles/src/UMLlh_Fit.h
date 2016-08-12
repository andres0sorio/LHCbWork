#ifndef UMLLH_FIT_H
#define UMLLH_FIT_H

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
#include "RandomDists.h"
#include "PDFs.h"
#include "Integrate.h"



///////////////////


void stepOne(const char *, int );

void stepTwo(const char *, int );

void getData(const char *);

void fcn(int &npar, double *grad, double &fval, double *par, int iflag);

#endif
