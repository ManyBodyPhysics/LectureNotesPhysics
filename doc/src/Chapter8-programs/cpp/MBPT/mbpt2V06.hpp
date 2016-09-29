#ifndef MBPT2V06_H
#define MBPT2V06_H
#include "parameters.hpp"
#include "modelspace.hpp"
#include "mbSolver.hpp"

class mbpt2V06 : public mbSolver {
public:
  mbpt2V06(Input_Parameters &parameters, Model_Space &space);
  double getEnergy();
private:
  Model_Space modelspace;
  int number_of_occupied_states;
  int number_of_unoccupied_states;
  double L;
};

#endif
