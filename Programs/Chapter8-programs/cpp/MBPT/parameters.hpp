#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdlib.h>

typedef struct Input_Parameters{
  int Nshells; //number of neutrons shells
  int Pshells; //number of protons shells
  int N; //number of neutrons
  int P; //number of protons
  double density;
  int Shells; // total number of shells
  int MBPT_Approx; // 2 for MBPT2, 3 for MBPT3_pp
  int MBPT_Function; // 0 for serial, 1 for parallel, 2 for block serial, 3 for block parallel, 4 for block M-M

  int numberOfProtons;
  int numberOfNeutrons;
  int numberOfParticles;

} Input_Parameters;

Input_Parameters parse_commandline_arguments(int argc, char *argv[]);
#endif
