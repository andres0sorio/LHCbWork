#ifndef UMLLH_FIT_H
#define UMLLH_FIT_H

#include "ThreeAnglesCommon.h"

#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TKey.h>

#include <TString.h>
#include <TMath.h>
#include <TMinuit.h>

#include "Utilities.h"
#include "PDFs.h"
#include "PDFsWRes.h"

///////////////////

void stepOne(const char *, int , int );

int  getData(const char *);

//........ 

void fcnwresol(int &npar, double *grad, double &fval, double *par, int iflag);


#endif
