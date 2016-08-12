#include "Integrate.h"

//std::ofstream *os = new std::ofstream("data.out",ofstream::out);

double integrate( Ptr2Function func , double *x , double *pp , double xmin , double xmax ) 
{
  
  std::vector<double> fvalues;
  
  double xx = xmin; 
  int np    = 100; //
  double step = (xmax-xmin) / (double)np;
  
  for( int i=0; i < np; ++i) {
    
    xx = xmin + i * step;
    x[0] = xx;
    double val = func( x , pp );
    fvalues.push_back(val);
    
    //     (*os) << xx  << '\t' 
    // 	  << val
    // 	  << std::endl;
    
  }
  
  double val = simpson( fvalues, xmin, xmax);
  
  return val;
  
}

double projectTo1D( Ptr2Function func , double *x , double *pp , double *xmin , double *xmax ) 
{
  
  //Projection: performs double integral on a three variable function to 
  //Notice: we leave x[0] free

  std::vector<double> fxvalues;
  std::vector<double> fyvalues;
  
  double xx = xmin[0];
  double yy = xmin[1];
  int np    = 100; //
  
  double stepx = (xmax[0]-xmin[0]) / (double)np;
  double stepy = (xmax[1]-xmin[1]) / (double)np;
  
  double val(0.);
  
  for( int i=0; i < np; ++i) {
    
    xx   = xmin[0] + i * stepx;
    x[1] = xx;
    
    for ( int j=0; j < np; ++j ) {
      yy = xmin[1] + j * stepy;
      x[2] = yy;
      
      double val = func( x , pp );
      fyvalues.push_back(val);
      //    (*os) << xx  << '\t' 
      //      << val
      //      << std::endl;
      
    }
    
    val = simpson( fyvalues, xmin[1], xmax[1]);
    fxvalues.push_back(val);
    fyvalues.clear();
  }
  
  val = simpson( fxvalues, xmin[0], xmax[0]);
  
  //std::cout << "With simpson: " << val << std::endl;
  
  return val;
  
}

double integrate3d( Ptr2Function func , double *x , double *pp , double *xmin , double *xmax ) 
{
  
  //Extending integrate to three dimensions
  
  std::vector<double> fxvalues;
  std::vector<double> fyvalues;
  std::vector<double> fzvalues;
  
  double xx = xmin[0];
  double yy = xmin[1];
  double zz = xmin[2];
  
  int np    = 10; //
  
  double stepx = (xmax[0]-xmin[0]) / (double)np;
  double stepy = (xmax[1]-xmin[1]) / (double)np;
  double stepz = (xmax[1]-xmin[1]) / (double)np;
  
  double val(0.);
  
  for( int i=0; i < np; ++i) {
    xx   = xmin[0] + i * stepx;
    x[0] = xx;
    
    for ( int j=0; j < np; ++j ) {
      yy = xmin[1] + j * stepy;
      x[1] = yy;
      
      for ( int k=0; k < np; ++k ) {
        zz = xmin[2] + k * stepz;
        x[2] = zz;
        double val = func( x , pp );
        fzvalues.push_back(val);
      }
      
      val = simpson( fzvalues, xmin[2], xmax[2]);
      fyvalues.push_back(val);
      fzvalues.clear();
    }
    val = simpson( fyvalues, xmin[1], xmax[1]);
    fxvalues.push_back(val);
    fyvalues.clear();
  }
  
  val = simpson( fxvalues, xmin[0], xmax[0]);
  
  return val;
  
}


double simpson(std::vector<double>& y, double xlo, double xhi)
{
  
  int n = y.size();
  double val;
  if (n < 3)
    val = -1;
  else {
    double dx = (xhi-xlo)/n;
    // the integral between the first two points is calculated by
    // (5*y[0]+8*y[1]-y[2])*dx/12
    val = (5*y[0]+8*y[1]-y[2])*dx/12;
    // the integral between the last two points is calculated by
    // (-y[n-3]+8*y[n-2]+5*y[n-1])*dx/12
    val = val + (-y[n-3]+8*y[n-2]+5*y[n-1])*dx/12;
    // the rest is calculated by the well-known Simpson's rule
    // (y[i]+4*y[i+1]+y[i+2])*dx/3
    for (int i=0; i<n-2; i++) {
      val = val + (y[i]+4*y[i+1]+y[i+2])*dx/3;
    }    
  }
  
  // divide val by 2, since every interval is integrated twice.
  val = val/2.0;
  
  return val;
  
}


