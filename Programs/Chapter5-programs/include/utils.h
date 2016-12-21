#ifndef INCLUDED_UTILS
#define INCLUDED_UTILS

#define SQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

int is_power(int);
double InnerProduct(double*, double*, int);

#endif
