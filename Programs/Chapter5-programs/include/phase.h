using namespace std;
#include "arg.h"
#include <fftw3.h> 

#ifndef INCLUDED_PHASE
#define INCLUDED_PHASE

class Phase
{
  private:
    int vol;
    double* phase;

  public:
    Phase(double*); // Use for arbitrary phases (for instance, to position SHO wave function in the center of the lattice)
    Phase(); // Use for APBC phases
    ~Phase();
    void PhaseGen(double*);

    void Add(double*);
    void Subtract(double*);

    void Add(fftw_complex*);
    void Subtract(fftw_complex*);

};

#endif
