#ifndef MODELSPACE_H
#define MODELSPACE_H

#include <iostream>
#include <cmath>
#include "parameters.hpp"
#include "state.hpp"

typedef struct Model_Space{
  int indp; //number of proton orbits
  int indn; //number of neutron orbits
  int indpar; //number of particle orbits
  int indhol; //number of hole orbits
  int indtot; //number of total orbits

  State *qnums;
  struct State qmins;
  struct State qsizes;

  int Nmax;
  int nmax;
  int *map_2b;
  int size_2b;
} Model_Space;

void Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space);

#endif
