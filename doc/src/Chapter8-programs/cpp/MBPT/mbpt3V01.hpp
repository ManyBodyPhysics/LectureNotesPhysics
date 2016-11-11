#ifndef MBPT3V01_H
#define MBPT3V01_H
#include "parameters.hpp"
#include "modelspace.hpp"
#include "mbSolver.hpp"

class mbpt3V01 : public mbSolver {
public:
  mbpt3V01(Input_Parameters &parameters, Model_Space &space);
  double getEnergy();
private:
  Model_Space modelspace;
  int number_of_occupied_states;
  int number_of_unoccupied_states;
  int number_of_states;
  double L;
};

#endif
