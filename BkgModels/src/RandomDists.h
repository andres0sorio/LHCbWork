// $Id: RandomDists.h,v 1.4 2007/02/18 23:53:03 aosorio Exp $
#ifndef RANDOMDISTS_H 
#define RANDOMDISTS_H 1

// Include files
#include "ThreeAnglesCommon.h"

#include <TFile.h>
#include <TDirectory.h>
#include <TList.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom3.h>

#include "PDFsBkg.h"

/** @class RandomDists RandomDists.h RandomDists.h
 *  
 *
 *  Toy Monte Carlo for Bs(anti-Bs) -> Jpsi phi
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

  void SelectSignalType();
    
  void RunExperiment();
  
  void WriteToTree();
  
  Ptr2pdf astPDF;
  double  fmax[2];
  double  pdfmax;
  
  double  qfactor;
  double  tag_eff;
  double  tag_frac;
  double  omega;
  bool    generatingBs;
  double  bkg_ov_sig;  // -> B/S ratio
  
  TRandom3 *rdn[4];    // one random number for each variable
  TRandom3 *reg;       // rescaling number
  TRandom3 *tag[4];    // consider tagging
  
  TFile *outfile;
  TTree *t1;
  
  double pi;
  double var[10];
  double tempoval[10];
  
  void   updateMaxima    ( double   );
  
protected:
  
private:
  
  
  
};

#endif // RANDOMDISTS_H
