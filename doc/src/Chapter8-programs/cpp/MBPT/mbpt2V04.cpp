#include "blas.hpp"
#include "mbpt2V04.hpp"
#include "minnesota_potential.hpp"

mbpt2V04::mbpt2V04(Input_Parameters &parameters, Model_Space &space, Channels &channels) {
  mbpt2V04::name = "mbpt2V04";
  mbpt2V04::modelspace = space;
  mbpt2V04::channels = channels;

  number_of_occupied_states = space.indhol;
  number_of_unoccupied_states = space.indpar;
  L = pow((parameters.numberOfParticles)/parameters.density, 1./3.);

}

double mbpt2V04::getEnergy() {
  double energy = 0.0;
  double *V1, *V2, *S1;
  char N = 'N';
  double fac0 = 0.0;
  double fac1 = 1.0;
  int nhh, npp, i, j, a, b, idx;
  double energy0;
  for(int chan = 0; chan < channels.size; ++chan){
    nhh = channels.nhh[chan];
    npp = channels.npp[chan];
    if(nhh*npp == 0){ continue; }

    V1 = new double[nhh * npp];
    V2 = new double[npp * nhh];
    S1 = new double[nhh * nhh];
    #pragma omp parallel shared(V1, V2) private(i, j, a, b, idx, energy0)
    {
      #pragma omp for schedule(static)
      for(int hh = 0; hh < nhh; ++hh){
        i = channels.hhvec[chan][2*hh];
        j = channels.hhvec[chan][2*hh + 1];
        for(int pp = 0; pp < npp; ++pp){
          a = channels.ppvec[chan][2*pp];
          b = channels.ppvec[chan][2*pp + 1];
          idx = hh * npp + pp;
          energy0 = V_Minnesota(modelspace, i, j, a, b, L);
          V1[idx] = energy0 / (modelspace.qnums[i].energy + modelspace.qnums[j].energy -
                modelspace.qnums[a].energy - modelspace.qnums[b].energy);
          idx = pp * nhh + hh;
          V2[idx] = energy0;
        }
      }
    }
    
    RM_dgemm(V1, V2, S1, &nhh, &nhh, &npp, &fac1, &fac0, &N, &N);
    delete V1; delete V2;

    for(int hh = 0; hh < nhh; ++hh){
      energy += S1[nhh*hh + hh];
    }
    delete S1;

  }
  energy *= 0.25;
  return energy;
}
