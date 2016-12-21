#include "mbpt3V05.hpp"
#include "minnesota_potential.hpp"
#include "blas.hpp"

mbpt3V05::mbpt3V05(Input_Parameters &parameters, Model_Space &space, Channels &channels) {
  mbpt3V05::name = "mbpt3V05";
  mbpt3V05::modelspace = space;
  mbpt3V05::channels = channels;

  number_of_occupied_states = space.indhol;
  number_of_unoccupied_states = space.indpar;
  L = pow((parameters.numberOfParticles)/parameters.density, 1./3.);

}

double mbpt3V05::getEnergy() {
  double energy = 0.0;
  double *V1, *S1, *V2, *S2;
  char N = 'N';
  char T = 'T';
  double fac0 = 0.0;
  double fac1 = 1.0;
  int nhh, npp, i, j, a, b, c, d, idx;
  double energy0;
  for(int chan = 0; chan < channels.size; ++chan){
    nhh = channels.nhh[chan];
    npp = channels.npp[chan];
    if(nhh*npp == 0){ continue; }

    V1 = new double[nhh * npp];
    S1 = new double[nhh * npp];
    V2 = new double[npp * npp];
    S2 = new double[nhh * nhh];
    for(int ab = 0; ab < npp; ++ab){
      a = channels.ppvec[chan][2*ab];
      b = channels.ppvec[chan][2*ab + 1];
      for(int hh = 0; hh < nhh; ++hh){
        i = channels.hhvec[chan][2*hh];
        j = channels.hhvec[chan][2*hh + 1];
        idx = hh * npp + ab;
        energy0 = V_Minnesota(modelspace, i, j, a, b, L);
        energy0 /= (modelspace.qnums[i].energy + modelspace.qnums[j].energy -
                    modelspace.qnums[a].energy - modelspace.qnums[b].energy);
        V1[idx] = energy0;
      }
      for(int cd = 0; cd < npp; ++cd){
        c = channels.ppvec[chan][2*cd];
        d = channels.ppvec[chan][2*cd + 1];
        idx = ab * npp + cd;
        V2[idx] = V_Minnesota(modelspace, a, b, c, d, L);
      }
    }

    RM_dgemm(V1, V2, S1, &nhh, &npp, &npp, &fac1, &fac0, &N, &N);
    RMT_dgemm(S1, V1, S2, &nhh, &nhh, &npp, &fac1, &fac0, &N, &T);
    delete V1; delete V2; delete S1;
    for(int hh = 0; hh < nhh; ++hh){
      energy += S2[nhh*hh + hh];
    }
    delete S2;
  }
  energy *= 0.125;

  return energy;
}
