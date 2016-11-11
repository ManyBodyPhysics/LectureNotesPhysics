#ifndef MBPT3V02_H
#define MBPT3V02_H
#include "parameters.hpp"
#include "modelspace.hpp"
#include "channels.hpp"
#include "mbSolver.hpp"

class mbpt3V02 : public mbSolver {
public:
  mbpt3V02(Input_Parameters &parameters, Model_Space &space, Channels &channels);
  double getEnergy();
private:
  Model_Space modelspace;
  Channels channels;
  int number_of_occupied_states;
  int number_of_unoccupied_states;
  double L;
};

#endif
