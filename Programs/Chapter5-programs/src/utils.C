#include "utils.h"

int is_power(int x)
{
    return !((x-1) & x);
}

double InnerProduct(double* a, double* b, int len)
// Given two real vectors a and b of length len, returns the inner product (a.b)
{
  double sum = 0.0;
  for(int i=0; i<len; i++){
    sum += a[i] * b[i];
  }
  return sum;
}

