#include <TROOT.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TSystem.h"

#include <iostream>
#include <sstream>
#include <cstdio>
#include <iomanip>
#include <cmath>

void setAxeOptions(TAxis *);

void setStyleOptions(TStyle *);

void createGIF(TString );
