#include "TROOT.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFrame.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TPaveLabel.h"
#include "TLine.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

int findMaximum( int , const std::vector<double> &);
int findMaximum( int , const std::vector<int> &);
void setGraphOptions( TGraph *, char *, char *);
void setAxeOptions(TAxis *);
void setFrameOptions( TPad *, char *, char *);
void createLabels ( TPad *, TPaveLabel **, std::string *, int);
void joinLabels ( TPad * , TPaveLabel *, TLine **, int ); 
void setLabelOptions ( TPaveLabel *);
void setLineOptions ( TLine *);
void setCutNames(std::string *);
