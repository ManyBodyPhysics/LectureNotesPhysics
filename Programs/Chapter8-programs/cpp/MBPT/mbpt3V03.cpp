#include "mbpt3V03.hpp"
#include "minnesota_potential.hpp"

mbpt3V03::mbpt3V03(Input_Parameters &parameters, Model_Space &space, Channels &channels) {
  mbpt3V03::name = "mbpt3V03";
  mbpt3V03::modelspace = space;
  mbpt3V03::channels = channels;

  number_of_occupied_states = space.indhol;
  number_of_unoccupied_states = space.indpar;
  L = pow((parameters.numberOfParticles)/parameters.density, 1./3.);

}

double mbpt3V03::getEnergy() {
  double energy = 0.0;
  double energy0, energy1;
  int nhh, npp, i, j, a, b, c, d;
  for(int chan = 0; chan < channels.size; ++chan){
    nhh = channels.nhh[chan];
    npp = channels.npp[chan];
    if(nhh*npp == 0){ continue; }

    #pragma omp parallel private(i, j, a, b, c, d, energy0, energy1)
    {
      #pragma omp for schedule(static) reduction(+:energy)
      for(int hh = 0; hh < nhh; ++hh){
        i = channels.hhvec[chan][2*hh];
        j = channels.hhvec[chan][2*hh + 1];
        for(int ab = 0; ab < npp; ++ab){
          a = channels.ppvec[chan][2*ab];
          b = channels.ppvec[chan][2*ab + 1];
          energy0 = V_Minnesota(modelspace, i, j, a, b, L);
          energy0 /= (modelspace.qnums[i].energy + modelspace.qnums[j].energy -
                        modelspace.qnums[a].energy - modelspace.qnums[b].energy);
          for(int cd = 0; cd < npp; ++cd){
            c = channels.ppvec[chan][2*cd];
            d = channels.ppvec[chan][2*cd + 1];
            energy1 = V_Minnesota(modelspace, a, b, c, d, L);
            energy1 *= V_Minnesota(modelspace, c, d, i, j, L);
            energy1 /= (modelspace.qnums[i].energy + modelspace.qnums[j].energy -
                          modelspace.qnums[c].energy - modelspace.qnums[d].energy);
            energy += energy0*energy1;
          }
        }
      }
    }
  }
  energy *= 0.125;
  return energy;
}
