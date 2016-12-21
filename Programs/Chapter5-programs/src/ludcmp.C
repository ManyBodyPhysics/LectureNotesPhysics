#include <iostream>
#include <complex>
#include <stdlib.h>
using namespace std;
//#include "math.h"
#define TINY 1.0e-20

void ludcmp(complex<double> *a, int n, int *indx, double *d)
{
// Modified version of Numerical Recipes in C ludcmp code.
  int i,imax,j,k;
  double big;
  complex<double> dum;
  complex<double> sum;
  double temp;
  double *vv;

  vv=(double *)malloc(n*sizeof(double));
  *d=1.0;
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++)
      if ((temp=abs(a[n*i+j])) > big) big=temp;
    if (big == 0.0) {
      cout << "ludcmp(complex<double>*, int, int*, double*): singular matrix\n";
      exit(EXIT_FAILURE);
    }
    vv[i]=1.0/big;
  }
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      sum=a[n*i+j];
      for (k=0;k<i;k++) sum -= a[n*i+k]*a[n*k+j];
      a[n*i+j]=sum;
    }
    big=0.0;
    for (i=j;i<n;i++) {
      sum=a[n*i+j];
      for (k=0;k<j;k++)
        sum -= a[n*i+k]*a[n*k+j];
      a[n*i+j]=sum;
      if ( real(dum=vv[i]*abs(sum)) >= big) {
        big=real(dum);
        imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<n;k++) {
        dum=a[n*imax+k];
        a[n*imax+k]=a[n*j+k];
        a[n*j+k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[n*j+j] == 0.0) a[n*j+j]=TINY;
    if (j != n) {
      dum=complex<double>(1.0,0.0)/a[n*j+j];
      for (i=j+1;i<n;i++) a[n*i+j] *= dum;
    }
  }
  free(vv);
}
