#include "parameters.hpp"

Input_Parameters parse_commandline_arguments(int argc, char *argv[]) {
  Input_Parameters Parameters;

  Parameters.density = atof(argv[1]);
  Parameters.Shells = atoi(argv[2]);
  Parameters.Pshells = atoi(argv[3]);
  Parameters.Nshells = atoi(argv[4]);
  Parameters.MBPT_Approx = atoi(argv[5]);
  Parameters.MBPT_Function = atoi(argv[6]);

  return Parameters;
}

