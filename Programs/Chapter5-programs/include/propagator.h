#include <vector>
#include "hamiltonian.h"
#include "arg.h"

#ifndef INCLUDED_FERMION_PROP
#define INCLUDED_FERMION_PROP

class Propagator
{
  private:
    double* b;
    int vol;

  public:
    Propagator(PropagatorArg);
    ~Propagator();
    void Run( vector<Hamiltonian*> );
    void SetSource(double*);
    double* Get();
    void Normalize();
    void Normalize(double);
};

#endif
