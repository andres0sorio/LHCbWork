// $Id: GetResults.h,v 1.6 2006/12/17 01:24:30 aosorio Exp $
#ifndef GETRESULTS_H 
#define GETRESULTS_H 1

// Include files
#include "Utilities.h"
#include "FitResults.h"
#include "TMatrixD.h"

/** @class GetResults GetResults.h
 *  
 *
 *  @author Andres Felipe OSORIO OLIVEROS
 *  @date   2006-11-21
 */


class GetResults {
public: 
  
  /// Standard constructor
  GetResults( ); 
  
  virtual ~GetResults( ); ///< Destructor
  
  void addDataFile( const char * ); //< Results from the fit
  
  void addDataList( const char * );
  
  void addDataSpec( const char * ); //< ToyMC input data
  
  void addHistogram( const char * , int , double , double);
  
  void addFitFunction( const char * , double , double);
  
  void addSelection( int , const char *, double, double);
    
  void getFitResults();

  void evalGCorrel( int );
  
  void drawHistograms();
    
  void plotFittedValue( int , int );
  
  void plotFittedError( int , int );
  
  void plotPull       ( int , int );

  void GetFailedExps();
  
  //void saveHistograms();
  
  FitResults * fitres;
  
  TCanvas *acanv;

  TCanvas *plotarea[10];
    
protected:
  

  
private:
  
  TH1D * h1d;
  TF1 * ff;
  
  std::vector<TObject *> histograms;
  std::vector<TF1*> fitfunctions;
  std::vector<std::string> fitresults;
  std::vector<double> tmcinputs;
  std::vector<int> selpar;
  std::vector<std::string> selparname;
  std::vector<double> selparmin;
  std::vector<double> selparmax;
  
  int npars;
    
};
#endif // GETRESULTS_H
