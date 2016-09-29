#ifndef MBPT2FACTORY_H
#define MBPT2FACTORY_H
#include "parameters.hpp"
#include "modelspace.hpp"
#include "channels.hpp"
#include "mbSolver.hpp"
#include "twobody-potential.hpp"

mbSolver* mbptFactory(Input_Parameters &parameters);
mbSolver* mbpt2Factory(Input_Parameters &parameters, Model_Space &modelspace,
                              Channels &channels, twobodyPotential *potential);
mbSolver* mbpt3Factory(Input_Parameters &parameters, Model_Space &modelspace, Channels &channels);

#endif
