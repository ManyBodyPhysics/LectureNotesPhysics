#ifndef MBPT_H
#define MBPT_H
#include <cmath>
#include "twobody-potential.hpp"
#include "modelspace.hpp"
#include "mbSolver.hpp"

class mbpt2V05 : public mbSolver {
public:
  mbpt2V05(twobodyPotential &potential, Model_Space &space);
  double getEnergy();
private:
  double * setup_twobody_potential_hhpp();
  double * allocate_and_initialize_double_vector(double value);
  double * setup_vector_of_energy_denominators_hhpp ();
  double energy_from_vectors(double *v, double *denominators);

  twobodyPotential *potential;
  Model_Space modelspace;
  int number_of_elements;
  int number_of_occupied_states;
  int number_of_unoccupied_states;
};

#endif
