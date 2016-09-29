#include "mbpt_factory.hpp"
#include "modelspace.hpp"
#include "channels.hpp"
#include "minnesota_potential.hpp"
#include "mbpt2V00.hpp"
#include "mbpt2V01.hpp"
#include "mbpt2V02.hpp"
#include "mbpt2V03.hpp"
#include "mbpt2V04.hpp"
#include "mbpt2V05.hpp"
#include "mbpt2V06.hpp"
#include "mbpt3V00.hpp"
#include "mbpt3V01.hpp"
#include "mbpt3V02.hpp"
#include "mbpt3V03.hpp"
#include "mbpt3V04.hpp"
#include "timer.hpp"

bool needChannels(Input_Parameters &parameters) {
  if (parameters.MBPT_Approx == 3) {
    return parameters.MBPT_Function > 1;
  } else {
    return parameters.MBPT_Function > 1 && parameters.MBPT_Function < 5;
  }
}

mbSolver* mbptFactory(Input_Parameters &parameters) {

  timer time = timer();
  time.start();
  Model_Space modelspace;
  Build_Model_Space(parameters, modelspace);
  calculate_sp_energies(parameters, modelspace);
  double elapsed = time.get_elapsed();
  std::cout << "Time to setup modelspace: " << elapsed << std::endl;

  std::cout << "# of states: " << modelspace.indtot << ", # of protons = " << parameters.numberOfProtons
                                  << ", # of neutrons = " << parameters.numberOfNeutrons << std::endl; 
  time.start();
  Channels channels;
  if ( needChannels(parameters) ) {
    Setup_Channels_MBPT(parameters, modelspace, channels);
  }
  elapsed = time.get_elapsed();
  std::cout << "Time to setup channels: " << elapsed << std::endl;

  time.start();
  minnesotaPotential* potential = new minnesotaPotential(parameters.numberOfParticles, parameters.density);
  elapsed = time.get_elapsed();
  std::cout << "Time to setup potential: " << elapsed << std::endl;

  if ( parameters.MBPT_Approx == 2) {
    return mbpt2Factory(parameters, modelspace, channels, potential);
  } else if (parameters.MBPT_Approx == 3) {
    return mbpt3Factory(parameters, modelspace, channels);
  } else {
    return NULL;
  }
    
}

mbSolver* mbpt2Factory(Input_Parameters &parameters, Model_Space &modelspace, Channels &channels, twobodyPotential *potential) {

  mbSolver *solver = NULL;

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
    case 6:
    {
      solver = (mbSolver*) new mbpt2V06(parameters, modelspace);
      break;
    }
  }

  return solver;
}

mbSolver* mbpt3Factory(Input_Parameters &parameters, Model_Space &modelspace, Channels &channels) {
  mbSolver *solver = NULL;

  switch (parameters.MBPT_Function) {
    /* Too slow for test case
    case 0:
    {
      solver = (mbSolver*) new mbpt3V00(parameters, modelspace);
      break;
    }
    case 1:
    {
      solver = (mbSolver*) new mbpt3V01(parameters, modelspace);
      break;
    } */
    case 2:
    {
      solver = (mbSolver*) new mbpt3V02(parameters, modelspace, channels);
      break;
    }
    case 3:
    {
      solver = (mbSolver*) new mbpt3V03(parameters, modelspace, channels);
      break;
    }
    case 4:
    {
      solver = (mbSolver*) new mbpt3V04(parameters, modelspace, channels);
      break;
    }
  }
  return solver;
}
