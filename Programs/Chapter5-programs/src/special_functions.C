#include <math.h>
#include "special_functions.h"
#include "constants.h"
#include "utils.h" 
double HermiteH(int n, double x) // Returns h_n(x)
{

  double p1 = PIM4;
  double p2 = 0.0;
  double p3;

  for (int j=1; j<=n; j++) {
    p3 = p2;
    p2 = p1;
    p1 = x*sqrt(2.0/j)*p2 - sqrt((double)(j-1)/j)*p3;
  }

  return p1;

}

void HermiteH(double* h, int n, double x) { // Returns array h = { h_0(x), h_1(x), h_2(x), ...., h_{n-1}(x) }

  double a = 0.0;
  double b;

  h[0] = PIM4;
  for (int j=1; j<n; j++) {
    b = a;
    a = h[j-1];
    h[j] = x*sqrt(2.0/j)*a - sqrt((double)(j-1)/j)*b;
  }

}




#define NMAX 6
#define H 0.4
#define A1 (2.0/3.0)
#define A2 0.4
#define A3 (2.0/7.0)
static double dsqrarg;

double DawsonF(double x)
// Numerical recipes in c;
// Returns Dawson’s integral function for any real x.
// This function is limited to single precision
{
  int i, n0;
  double d1,d2,e1,e2,sum,x2,xp,xx,ans;
  static double c[NMAX+1];
  static int init = 0; //  Flag is 0 if we need to initialize, else 1.
  if (init == 0) {
    init=1;
    for (i=1;i<=NMAX;i++) c[i]=exp(-SQR((2.0*i-1.0)*H));
  }
  if (fabs(x) < 0.2) { //  Use series expansion.
    x2=x*x;
    ans=x*(1.0-A1*x2*(1.0-A2*x2*(1.0-A3*x2)));
  } else {  //  Use sampling theorem representation.
    xx=fabs(x);
    n0=2*(int)(0.5*xx/H+0.5);
    xp=xx-n0*H;
    e1=exp(2.0*xp*H);
    e2=e1*e1;
    d1=n0+1;
    d2=d1-2.0;
    sum=0.0;
    for (i=1;i<=NMAX;i++,d1+=2.0,d2-=2.0,e1*=e2)
      sum += c[i]*(e1/d1+1.0/(d2*e1));
    ans=SIGN(exp(-xp*xp),x)*sum/SQRTPI; // Constant is 1/sqrt(pi)
  }
  return ans;
}
#undef NMAX
#undef H
#undef A1
#undef A2
#undef A3

