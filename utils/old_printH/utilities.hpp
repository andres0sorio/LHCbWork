#include <TROOT.h>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TKey.h"
#include "TString.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TObjArray.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include <sys/stat.h>
#include <sys/unistd.h>

#include "plot_util.hpp"

std::string getBaseName(const char *);

void printHistograms ( TDirectory * );


