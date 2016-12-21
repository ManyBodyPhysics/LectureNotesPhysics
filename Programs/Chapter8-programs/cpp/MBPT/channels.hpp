#ifndef CHANNELS_H
#define CHANNELS_H

#include "parameters.hpp"
#include "state.hpp"
#include "modelspace.hpp"

typedef struct Channels{
  int size;

  int* nhh;
  int* npp;
  int** hhvec;
  int** ppvec;
} Channels;

void Setup_Channels_MBPT(const Input_Parameters &Parameters, const Model_Space &Space, Channels &Chan);
int Chan_2bInd(const Model_Space &Space, const State &State);
#endif
