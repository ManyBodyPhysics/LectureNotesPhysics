#ifndef MBPT3V04_H
#define MBPT3V04_H
#include "parameters.hpp"
#include "modelspace.hpp"
#include "channels.hpp"
#include "mbSolver.hpp"

class mbpt3V04 : public mbSolver {
public:
  mbpt3V04(Input_Parameters &parameters, Model_Space &space, Channels &channels);
  double getEnergy();
private:
  Model_Space modelspace;
  Channels channels;
  int number_of_occupied_states;
  int number_of_unoccupied_states;
  double L;
};

#endif
