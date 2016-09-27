#include "mbpt2_factory.hpp"
#include "modelspace.hpp"
#include "channels.hpp"
#include "minnesota_potential.hpp"
#include "mbpt2V00.hpp"
#include "mbpt2V01.hpp"
#include "mbpt2V02.hpp"
#include "mbpt2V03.hpp"
#include "mbpt2V04.hpp"
#include "mbpt2V05.hpp"

bool needChannels(Input_Parameters &parameters) {
  return  parameters.MBPT_Function > 1 && parameters.MBPT_Function < 5;
}

mbSolver* mbpt2Factory(Input_Parameters &parameters) {

  mbSolver *solver;
  Model_Space modelspace;
  Build_Model_Space(parameters, modelspace);
  calculate_sp_energies(parameters, modelspace);

  Channels channels;
  if ( needChannels(parameters) ) {
    Setup_Channels_MBPT(parameters, modelspace, channels);
  }

  minnesotaPotential* potential = new minnesotaPotential(parameters.numberOfParticles, parameters.density);

  switch (parameters.MBPT_Function) {
  case 0:
  {
    solver = (mbSolver*) new mbpt2V00(*potential, modelspace);
    break;
  }
  case 1:
  {
    solver = (mbSolver*) new mbpt2V01(*potential, modelspace);
    break;
  }
  case 2:
  {
    solver = (mbSolver*) new mbpt2V02(parameters, modelspace, channels);
    break;
  }
  case 3:
  {
    solver = (mbSolver*) new mbpt2V03(parameters, modelspace, channels);
    break;
  }
  case 4:
  {
    solver = (mbSolver*) new mbpt2V04(parameters, modelspace, channels);
    break;
  }
  case 5:
  {
    solver = (mbSolver*) new mbpt2V05(*potential, modelspace);
    break;
  }
  }

  return solver;
}
