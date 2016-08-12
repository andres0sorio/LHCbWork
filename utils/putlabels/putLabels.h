// $Id: putLabels.h,v 1.1 2006/03/26 15:03:26 aooliver Exp $
#ifndef PUTLABELS_PUTLABELS_H 
#define PUTLABELS_PUTLABELS_H 1

// Include files
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TFile.h"
#include "TObject.h"
#include "TKey.h"

#include "TNamed.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPaveLabel.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

/** @class putLabels putLabels.h putLabels/putLabels.h
 * With this little class I hope I can read a Gaudi histograms
 * and add labels 
 *
 *  @author Andres Osorio Oliveros
 *  @date   2006-03-26
 */

struct Options {
  friend std::fstream& operator>>(std::fstream &istr, Options &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, Options &rhs);
  
  int id;
  std::string xLabel;
  std::string yLabel;
  
  Options () {};
  ~Options () {};
  
};


class putLabels {
 
 public: 
  /// Standard constructor
  putLabels(const char*, const char* ); 
  
  virtual ~putLabels( ); ///< Destructor

  void loopOverHistograms(TDirectory *,TDirectory *); 

  void processFile();

  inline void setOptions(TH1D *);
  inline void setOptions(TH2D *);

  protected:
  
 private:
  
  std::fstream *opts;
  
  TFile *input;
  
  TFile *output;
  
  std::vector<Options*> vopts;
  
};
#endif // PUTLABELS_PUTLABELS_H
