#ifndef MBPT2V00_H
#define MBPT2V00_H
#include "twobody-potential.hpp"
#include "modelspace.hpp"
#include "mbSolver.hpp"

class mbpt2V00 : public mbSolver {
public:
  mbpt2V00(twobodyPotential &potential, Model_Space &space);
  double getEnergy();
private:
  twobodyPotential *potential;
  Model_Space modelspace;
  int number_of_occupied_states;
  int number_of_unoccupied_states;
};

#endif
