// $Id: GetResults.h,v 1.4 2006/11/27 14:26:02 aosorio Exp $
#ifndef GETRESULTS_H 
#define GETRESULTS_H 1

// Include files
#include "Utilities.h"
#include "FitResults.h"

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
  
  void getFitResults();
  
  void plotFittedValue( int , int );
  
  void plotFittedError( int , int );
  
  void plotPull       ( int , int );

  void GetFailedExps();
  
  //void saveHistograms();
  
  FitResults * fitres;
  
  TCanvas *acanv;
  
protected:
  

  
private:
  
  TH1D * h1d;
  std::vector<TObject *> histograms;
  std::vector<std::string> fitresults;
  std::vector<double> tmcinputs;
  
  int hpos;
    
};
#endif // GETRESULTS_H
