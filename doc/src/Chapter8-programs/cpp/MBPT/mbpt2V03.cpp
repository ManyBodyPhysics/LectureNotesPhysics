#include "mbpt2V03.hpp"
#include "minnesota_potential.hpp"

mbpt2V03::mbpt2V03(Input_Parameters &parameters, Model_Space &space, Channels &channels) {
  mbpt2V03::modelspace = space;
  mbpt2V03::channels = channels;

  number_of_occupied_states = space.indhol;
  number_of_unoccupied_states = space.indpar;
  L = pow((parameters.numberOfParticles)/parameters.density, 1./3.);

}

double mbpt2V03::getEnergy() {
  double energy = 0.0;
  double energy0;
  int nhh, npp, i, j, a, b;
  for(int chan = 0; chan < channels.size; ++chan){
    nhh = channels.nhh[chan];
    npp = channels.npp[chan];
    if(nhh*npp == 0){ continue; }

    #pragma omp parallel private(i, j, a, b, energy0)
    {
      #pragma omp for schedule(static) reduction(+:energy)
      for(int hh = 0; hh < nhh; ++hh){
        i = channels.hhvec[chan][2*hh];
        j = channels.hhvec[chan][2*hh + 1];
        for(int pp = 0; pp < npp; ++pp){
          a = channels.ppvec[chan][2*pp];
          b = channels.ppvec[chan][2*pp + 1];
          energy0 = V_Minnesota(modelspace, i, j, a, b, L);
          energy0 *= energy0;
          energy0 /= (modelspace.qnums[i].energy + modelspace.qnums[j].energy -
                        modelspace.qnums[a].energy - modelspace.qnums[b].energy);
          energy += energy0;
        }
      }
    }
  }
  energy *= 0.25;
  return energy;

}
