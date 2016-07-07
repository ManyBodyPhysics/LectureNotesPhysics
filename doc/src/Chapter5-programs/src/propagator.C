#include <iostream>
#include <vector>
#include <math.h>
using namespace std;
#include "dispersion.h"
#include "hamiltonian.h"
#include "propagator.h"
#include "arg.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"
#include "utils.h"

Propagator::Propagator(PropagatorArg propagator_arg)
{

  char* fname = "Propagator::Propagator(PropagatorArg)";

  vol = GJP.Vol();
  b = (double *) malloc(2*vol*sizeof(double));  

}

Propagator::~Propagator()
{

  free(b);

}

void Propagator::Run(vector<Hamiltonian*> hamiltonians)
{
  char* fname = "void Propagator::Run(vector<Hamiltonian*>)";
  VRB.Func(fname);

  //---- Apply product ofinteraction matrices, b -> e^{-H}.b
  for (int i=0; i<hamiltonians.size(); i++) {
    hamiltonians[i]->Evolve(b);
  }

}

void Propagator::SetSource(double* psi)
{
  char* fname = "void Propagator::SetSource(double*)";
  VRB.Func(fname);

  for (int i=0; i<2*vol; i++) {
    b[i] = psi[i];
  }

}

double* Propagator::Get()
{

  char* fname = "double* Propagator::Get()";
  VRB.Func(fname);

  return b;

}

void Propagator::Normalize()
{

  char* fname = "void Propagator::Normalize()";
  VRB.Func(fname);

  double norm = 0.0;

  for (int i=0; i<2*vol; i++) {
    norm += b[i]*b[i];;
  }

  norm = sqrt(norm);

  for (int i=0; i<2*vol; i++) {
    b[i] /= norm;
  }


}

void Propagator::Normalize(double norm)
{

  char* fname = "void Propagator::Normalize(double)";
  VRB.Func(fname);

  for (int i=0; i<2*vol; i++) {
    b[i] /= norm;
  }

}
