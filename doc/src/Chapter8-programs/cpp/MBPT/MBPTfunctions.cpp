#include "MBPTfunctions.hpp"
#include "blas.hpp"
#include "minnesota_potential.hpp"

void MBPT3_0(const Input_Parameters &Parameters, const Model_Space &Space)
{
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);
  double energy0, energy1;
  double energy = 0.0;
  for(int i = 0; i < Space.indhol; ++i){
    for(int j = 0; j < Space.indhol; ++j){
      if(i == j){ continue; }
      for(int a = Space.indhol; a < Space.indtot; ++a){
	for(int b = Space.indhol; b < Space.indtot; ++b){
	  if(a == b){ continue; }
	  energy0 = V_Minnesota(Space, i, j, a, b, L);
	  energy0 /= (Space.qnums[i].energy + Space.qnums[j].energy - Space.qnums[a].energy - Space.qnums[b].energy);
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    for(int d = Space.indhol; d < Space.indtot; ++d){
	      if(c == d){ continue; }
	      energy1 = V_Minnesota(Space, a, b, c, d, L);
	      energy1 *= V_Minnesota(Space, c, d, i, j, L);
	      energy1 /= (Space.qnums[i].energy + Space.qnums[j].energy - Space.qnums[c].energy - Space.qnums[d].energy);
	      energy += energy0 * energy1;
	    }
	  }
	}
      }
    }
  }
  energy *= 0.125;
  std::cout << "DeltaE for MBPT(3)_0 = " << energy << std::endl;
}

void MBPT3_1(const Input_Parameters &Parameters, const Model_Space &Space)
{
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);
  double energy = 0.0;
  double energy0, energy1;
  #pragma omp parallel private(energy0, energy1)
  {
    #pragma omp for schedule(static) reduction(+:energy)
    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	if(i == j){ continue; }
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    if(a == b){ continue; }
	    energy0 = V_Minnesota(Space, i, j, a, b, L);
	    energy0 /= (Space.qnums[i].energy + Space.qnums[j].energy - Space.qnums[a].energy - Space.qnums[b].energy);
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(c == d){ continue; }
		energy1 = V_Minnesota(Space, a, b, c, d, L);
		energy1 *= V_Minnesota(Space, c, d, i, j, L);
		energy1 /= (Space.qnums[i].energy + Space.qnums[j].energy - Space.qnums[c].energy - Space.qnums[d].energy);
		energy += energy0 * energy1;
	      }
	    }
	  }
	}
      }
    }
  }
  energy *= 0.125;
  std::cout << "DeltaE for MBPT(3)_1 = " << energy << std::endl;
}

void MBPT3_2(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);
  double energy = 0.0;
  double energy0, energy1;
  int nhh, npp, i, j, a, b, c, d;
  for(int chan = 0; chan < Chan.size; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh*npp == 0){ continue; }

    for(int hh = 0; hh < nhh; ++hh){
      i = Chan.hhvec[chan][2*hh];
      j = Chan.hhvec[chan][2*hh + 1];
      for(int pp0 = 0; pp0 < npp; ++pp0){
	a = Chan.ppvec[chan][2*pp0];
	b = Chan.ppvec[chan][2*pp0 + 1];
	energy0 = V_Minnesota(Space, i, j, a, b, L);
	energy0 /= (Space.qnums[i].energy + Space.qnums[j].energy - Space.qnums[a].energy - Space.qnums[b].energy);
	for(int pp1 = 0; pp1 < npp; ++pp1){
	  c = Chan.ppvec[chan][2*pp1];
	  d = Chan.ppvec[chan][2*pp1 + 1];
	  energy1 = V_Minnesota(Space, a, b, c, d, L);
	  energy1 *= V_Minnesota(Space, c, d, i, j, L);
	  energy1 /= (Space.qnums[i].energy + Space.qnums[j].energy - Space.qnums[c].energy - Space.qnums[d].energy);
	  energy += energy0 * energy1;
	}
      }
    }
  }
  energy *= 0.125;
  std::cout << "DeltaE for MBPT(3)_2 = " << energy << std::endl;
}

void MBPT3_3(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);
  double energy = 0.0;
  double energy0, energy1;
  int nhh, npp, i, j, a, b, c, d;
  for(int chan = 0; chan < Chan.size; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh*npp == 0){ continue; }

    #pragma omp parallel private(i, j, a, b, c, d, energy0, energy1)
    {
      #pragma omp for schedule(static) reduction(+:energy)
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	for(int pp0 = 0; pp0 < npp; ++pp0){
	  a = Chan.ppvec[chan][2*pp0];
	  b = Chan.ppvec[chan][2*pp0 + 1];
	  energy0 = V_Minnesota(Space, i, j, a, b, L);
	  energy0 /= (Space.qnums[i].energy + Space.qnums[j].energy - Space.qnums[a].energy - Space.qnums[b].energy);
	  for(int pp1 = 0; pp1 < npp; ++pp1){
	    c = Chan.ppvec[chan][2*pp1];
	    d = Chan.ppvec[chan][2*pp1 + 1];
	    energy1 = V_Minnesota(Space, a, b, c, d, L);
	    energy1 *= V_Minnesota(Space, c, d, i, j, L);
	    energy1 /= (Space.qnums[i].energy + Space.qnums[j].energy - Space.qnums[c].energy - Space.qnums[d].energy);
	    energy += energy0 * energy1;
	  }
	}
      }
    }
  }
  energy *= 0.125;
  std::cout << "DeltaE for MBPT(3)_3 = " << energy << std::endl;
}

void MBPT3_4(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);
  double energy = 0.0;
  double *V1, *S1, *V2, *S2;
  char N = 'N';
  char T = 'T';
  double fac0 = 0.0;
  double fac1 = 1.0;
  int nhh, npp, i, j, a, b, c, d, idx;
  double energy0;
  for(int chan = 0; chan < Chan.size; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh*npp == 0){ continue; }

    V1 = new double[nhh * npp];
    S1 = new double[nhh * npp];
    V2 = new double[npp * npp];
    S2 = new double[nhh * nhh];
    #pragma omp parallel shared(V1, V2) private(i, j, a, b, c, d, idx, energy0)
    {
      #pragma omp for schedule(static)
      for(int pp0 = 0; pp0 < npp; ++pp0){
	a = Chan.ppvec[chan][2*pp0];
	b = Chan.ppvec[chan][2*pp0 + 1];
	for(int hh = 0; hh < nhh; ++hh){
	  i = Chan.hhvec[chan][2*hh];
	  j = Chan.hhvec[chan][2*hh + 1];
	  idx = hh * npp + pp0;
	  energy0 = V_Minnesota(Space, i, j, a, b, L);
	  energy0 /= (Space.qnums[i].energy + Space.qnums[j].energy - Space.qnums[a].energy - Space.qnums[b].energy);
	  V1[idx] = energy0;
	}
	for(int pp1 = 0; pp1 < npp; ++pp1){
	  c = Chan.ppvec[chan][2*pp1];
	  d = Chan.ppvec[chan][2*pp1 + 1];
	  idx = pp0 * npp + pp1;
	  V2[idx] = V_Minnesota(Space, a, b, c, d, L);
	}
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
  std::cout << "DeltaE for MBPT(3)_4 = " << energy << std::endl;
}
