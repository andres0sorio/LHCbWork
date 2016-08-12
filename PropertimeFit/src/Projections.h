// $Id: Projections.h,v 1.4 2007/04/12 11:11:21 aosorio Exp $
#ifndef PROJECTIONS_H 
#define PROJECTIONS_H 1

// Include files
#include "ThreeAnglesCommon.h"
#include "Utilities.h"

/** @class Projections Projections.h
 *  
 *
 *  @author Andres Osorio
 *  @date   2007-03-18
 */

typedef double (*ptr2Function) (double *, double *);

class Projections {
public: 
  /// Standard constructor
  Projections( ) { }; 
  Projections( const char * , int , double, double); 

  virtual ~Projections( ); ///< Destructor
  
  void Initialise (  double , bool );
  
  void AddFunction( TF1 * , int );
    
  void AddFunction( const char *, int , ptr2Function );
  
  void SetParameters( const char * , double );
  
  void SetTimeLimits( double , double );
  
  void Draw( );
  
  void DrawLog( );
    
  TH1D *dist;
  TF1 *ff;
  TCanvas *acanv;
  std::vector<TF1 *> curves;
  TStyle *cstyle;
      
protected:
  
private:
  
  std::string name;
  double sfactor1;
  double area1;
  double xmin;
  double xmax;
  int mbins;
  int maxparameters;
  
  std::vector<double> v1;
  std::vector<double> v2;

  int colours[10];
  int linesty[10];
  
};
#endif // PROJECTIONS_H
