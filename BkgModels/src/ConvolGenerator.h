// $Id: ConvolGenerator.h,v 1.3 2007/02/11 23:50:18 aosorio Exp $
#ifndef CONVOLGENERATOR_H 
#define CONVOLGENERATOR_H 1

// Include files
#include "RandomDists.h"
#include "ResModels.h"

/** @class ConvolGenerator ConvolGenerator.h
 *  
 *  implements Convolution to the already existing class RandomDists
 *  @author Andres Osorio
 *  @date   2007-01-17
 */

typedef double (*Ptr2ResModel) ( double (*) (double *, double *), 
                                 double *, double *, double *);
typedef double (*Ptr2pdfwconv) ( Ptr2ResModel, double *, double *, double *);

class ConvolGenerator : public RandomDists {
public: 
  /// Standard constructor
  ConvolGenerator( ); 
  
  ConvolGenerator(const char *, int ); 
  
  virtual ~ConvolGenerator( ); ///< Destructor
  
  void initialise(double *, int );
  
  void WriteToTree();
  
  void RunExperiment();
  

protected:
  


private:

  double convpars[8];
  double sigmat;

  Ptr2pdfwconv astPDFwConv;

   
};

#endif // CONVOLGENERATOR_H
