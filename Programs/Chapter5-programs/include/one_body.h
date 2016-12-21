#include <complex>
#include "arg.h"

#ifndef INCLUDED_ONE_BODY
#define INCLUDED_ONE_BODY

class OneBody
{
  private:
    SourceType source_type;
    double lambda1;
    double lambda2;
    double lambda3;
    double* psi; // Source vector
    int fermion_size;

  public:
    OneBody(OneBodyArg); // Specify basis, lambda parameters, etc
    ~OneBody();
    void Set(string); // Set source based on parameters specified in a file
    void Set(); // Set source based on what is stored in psi and specified source_type.
    void Set(int, int, int); // Set source to basis vector labeled by ( coord1, coord2, coord3)
    void Set(double*); // Set source to external vector
    complex<double> Project(double*); // perform projection of argument vertor onto source vector
    double* Get(); // Get the address to the source vector
};

#endif
