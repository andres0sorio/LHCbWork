#define FUNC(x) ((*func)(x))
#define EPS 1.0e-6
#define JMAX 2

double func (double y) {
  
  return TMath::Sin(y);

}


double trapzd(double (*func)(double), double a, double b, int n)
{
 
  double x(0.),tnm(0.),sum(0.),del(0.);
  
  static double s(0.);

  int it(0.),j(0.);

  if (n == 1) {
    return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));

  } else {

    for (it=1,j=1;j<n-1;j++) it <<= 1;
    
    printf("tnm: %d \n",it);

    tnm=it;
    del=(b-a)/tnm; // This is the spacing of the points to be added.
    x=a+0.5*del;
    
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
    
    s=0.5*(s+(b-a)*sum/tnm); // This replaces s by its refined value.
    return s;
  
  }
  
}


double qsimp(double (*func)(double), double a, double b)
{
  
  double s, st;
  double ost = 0.0;
  double os  = 0.0;

  for(int j = 1; j <= JMAX; ++j) {
    
    st = trapzd(func,a,b,j);
    
    s = (4.0 * st -ost) / 3.0;

    if ( j > 5 )
      if ( fabs(s-os) < EPS*fabs(os) || 
	   (s == 0.0 && os == 0.0)) return s;
    
    os  =  s;
    ost = st;

  }
  
  printf("Too many steps in routine qsimp\n");

  return 0.0;
  
  
}
