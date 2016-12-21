#include <complex>
#include "arg.h"

#ifndef INCLUDED_TWO_BODY
#define INCLUDED_TWO_BODY

class TwoBody
{
  private:
    double* wavefunc;
    int* opposite_parity_index;
    int vol;

  public:
    TwoBody(TwoBodyArg);
    ~TwoBody();
    void Deform(string); // Set source based on parameters specified in a file
    complex<double> Run(double*,double*);
};

#endif
