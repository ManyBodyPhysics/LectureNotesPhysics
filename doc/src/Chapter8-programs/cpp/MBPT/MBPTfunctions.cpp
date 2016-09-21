#include "MBPTfunctions.hpp"

void MBPT2_0(const Input_Parameters &Parameters, const Model_Space &Space)
{
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);
  double energy0;
  double energy = 0.0;
  for(int i = 0; i < Space.indhol; ++i){
    for(int j = 0; j < Space.indhol; ++j){
      if(i == j){ continue; }
      for(int a = Space.indhol; a < Space.indtot; ++a){
	for(int b = Space.indhol; b < Space.indtot; ++b){
	  if(a == b){ continue; }
	  energy0 = V_Minnesota(Space, i, j, a, b, L);
	  energy0 *= energy0;
	  energy0 /= (Space.qnums[i].energy + Space.qnums[j].energy - Space.qnums[a].energy - Space.qnums[b].energy);
	  energy += energy0;
	}
      }
    }
  }
  energy *= 0.25;
  std::cout << "DeltaE for MBPT(2)_0 = " << energy << std::endl;
}

void MBPT2_1(const Input_Parameters &Parameters, const Model_Space &Space)
{
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);
  double energy = 0.0;
  double energy0;
  #pragma omp parallel private(energy0)
  {
    #pragma omp for schedule(static) reduction(+:energy)
    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	if(i == j){ continue; }
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    if(a == b){ continue; }
	    energy0 = V_Minnesota(Space, i, j, a, b, L);
	    energy0 *= energy0;
	    energy0 /= (Space.qnums[i].energy + Space.qnums[j].energy - Space.qnums[a].energy - Space.qnums[b].energy);
	    energy += energy0;
	  }
	}
      }
    }
  }
  energy *= 0.25;
  std::cout << "DeltaE for MBPT(2)_1 = " << energy << std::endl;
}

void MBPT2_2(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);
  double energy = 0.0;
  double energy0;
  int nhh, npp, i, j, a, b;
  for(int chan = 0; chan < Chan.size; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh*npp == 0){ continue; }

    for(int hh = 0; hh < nhh; ++hh){
      i = Chan.hhvec[chan][2*hh];
      j = Chan.hhvec[chan][2*hh + 1];
      for(int pp = 0; pp < npp; ++pp){
	a = Chan.ppvec[chan][2*pp];
	b = Chan.ppvec[chan][2*pp + 1];
	energy0 = V_Minnesota(Space, i, j, a, b, L);
	energy0 *= energy0;
	energy0 /= (Space.qnums[i].energy + Space.qnums[j].energy - Space.qnums[a].energy - Space.qnums[b].energy);
	energy += energy0;
      }
    }
  }
  energy *= 0.25;
  std::cout << "DeltaE for MBPT(2)_2 = " << energy << std::endl;
}

void MBPT2_3(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);
  double energy = 0.0;
  double energy0;
  int nhh, npp, i, j, a, b;
  for(int chan = 0; chan < Chan.size; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh*npp == 0){ continue; }
    
    #pragma omp parallel private(i, j, a, b, energy0)
    {
      #pragma omp for schedule(static) reduction(+:energy)
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	for(int pp = 0; pp < npp; ++pp){
	  a = Chan.ppvec[chan][2*pp];
	  b = Chan.ppvec[chan][2*pp + 1];
	  energy0 = V_Minnesota(Space, i, j, a, b, L);
	  energy0 *= energy0;
	  energy0 /= (Space.qnums[i].energy + Space.qnums[j].energy - Space.qnums[a].energy - Space.qnums[b].energy);
	  energy += energy0;
	}
      }
    }
  }
  energy *= 0.25;
  std::cout << "DeltaE for MBPT(2)_3 = " << energy << std::endl;
}

void MBPT2_4(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);
  double energy = 0.0;
  double *V1, *V2, *S1;
  char N = 'N';
  double fac0 = 0.0;
  double fac1 = 1.0;
  int nhh, npp, i, j, a, b, idx;
  double energy0;
  for(int chan = 0; chan < Chan.size; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh*npp == 0){ continue; }

    V1 = new double[nhh * npp];
    V2 = new double[npp * nhh];
    S1 = new double[nhh * nhh];
    #pragma omp parallel shared(V1, V2) private(i, j, a, b, idx, energy0)
    {
      #pragma omp for schedule(static)
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	for(int pp = 0; pp < npp; ++pp){
	  a = Chan.ppvec[chan][2*pp];
	  b = Chan.ppvec[chan][2*pp + 1];
	  idx = hh * npp + pp;
	  energy0 = V_Minnesota(Space, i, j, a, b, L);
	  V1[idx] = energy0 / (Space.qnums[i].energy + Space.qnums[j].energy - Space.qnums[a].energy - Space.qnums[b].energy);
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
  std::cout << "DeltaE for MBPT(2)_4 = " << energy << std::endl;
}

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

void plus(State &S, const State &S1, const State &S2){
  S.t = S1.t + S2.t;
  S.m = S1.m + S2.m;
  S.nx = S1.nx + S2.nx;
  S.ny = S1.ny + S2.ny;
  S.nz = S1.nz + S2.nz;
}

void Setup_Channels_MBPT(const Input_Parameters &Parameters, const Model_Space &Space, Channels &Chan)
{
  State state0;
  int idx;
  Chan.size = Space.size_2b; // size1 = number of direct channels

  Chan.hhvec = new int*[Chan.size];
  Chan.ppvec = new int*[Chan.size];
  Chan.nhh = new int[Chan.size];
  Chan.npp = new int[Chan.size];

  for(int chan = 0; chan < Chan.size; ++chan){
    Chan.nhh[chan] = 0;
    Chan.npp[chan] = 0;
  }

  // count number of hh and pp pairs in each channel
  for(int i = 0; i < Space.indhol; ++i){
    for(int j = 0; j < Space.indhol; ++j){
      if(i == j){ continue; }
      plus(state0, Space.qnums[i], Space.qnums[j]);
      idx = Chan_2bInd(Space, state0);
      ++Chan.nhh[idx];
    }
  }
  for(int a = Space.indhol; a < Space.indtot; ++a){
    for(int b = Space.indhol; b < Space.indtot; ++b){
      if(a == b){ continue; }
      plus(state0, Space.qnums[a], Space.qnums[b]);
      idx = Chan_2bInd(Space, state0);
      ++Chan.npp[idx];
    }
  }

  // allocate memory for hhvec and ppvec which list hh and pp pairs
  for(int chan = 0; chan < Chan.size; ++chan){
    Chan.hhvec[chan] = new int[2 * Chan.nhh[chan]];
    Chan.ppvec[chan] = new int[2 * Chan.npp[chan]];
    Chan.nhh[chan] = 0;
    Chan.npp[chan] = 0;
  }

  // fill hhvec and ppvec in each channel
  for(int i = 0; i < Space.indhol; ++i){
    for(int j = 0; j < Space.indhol; ++j){
      if(i == j){ continue; }
      plus(state0, Space.qnums[i], Space.qnums[j]);
      idx = Chan_2bInd(Space, state0);
      Chan.hhvec[idx][2 * Chan.nhh[idx]] = i;
      Chan.hhvec[idx][2 * Chan.nhh[idx] + 1] = j;
      ++Chan.nhh[idx];
    }
  }
  for(int a = Space.indhol; a < Space.indtot; ++a){
    for(int b = Space.indhol; b < Space.indtot; ++b){
      if(a == b){ continue; }
      plus(state0, Space.qnums[a], Space.qnums[b]);
      idx = Chan_2bInd(Space, state0);
      Chan.ppvec[idx][2 * Chan.npp[idx]] = a;
      Chan.ppvec[idx][2 * Chan.npp[idx] + 1] = b;
      ++Chan.npp[idx];
    }
  }
}

void Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space)
{
  double E;
  double neutron_prefac = hbarc_MeVfm * hbarc_MeVfm / (2.0 * m_neutronc2);
  double proton_prefac = hbarc_MeVfm * hbarc_MeVfm / (2.0 * m_protonc2);
  
  int count = 0; // total state count
  int pcount = 0; // proton state count
  int ncount = 0; // neutron state count
  int holcount = 0; // hole state count
  int parcount = 0; // particle state count
  int phcount = 0; // proton hole state count
  int nhcount = 0; // neutron hole state count
  
  int shellnums [] = {1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 
	21, 22, 23, 25, 26, 27, 28, 30, 31, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46};

  // Find appropriate Nmax for the number of shells using shellnums
  if(Parameters.Shells > 40){ std::cerr << "Nmax too big!" << std::endl; exit(1); }
  Space.Nmax = shellnums[Parameters.Shells - 1];

  // Find maximum number of states (indtot) depending on Nmax and whether or not there are protons/neutrons
  if(Parameters.Pshells != 0 && Parameters.Nshells != 0){ 
    Space.indtot = (2*Space.Nmax + 1) * (2*Space.Nmax + 1) * (2*Space.Nmax + 1) * 2 * 2;
    Space.qmins.t = -1, Space.qsizes.t = 2;
  }
  else if(Parameters.Pshells != 0 && Parameters.Nshells == 0){
    Space.indtot = (2*Space.Nmax + 1) * (2*Space.Nmax + 1) * (2*Space.Nmax + 1) * 2;
    Space.qmins.t = -1, Space.qsizes.t = 1;
  }
  else if(Parameters.Pshells == 0 && Parameters.Nshells != 0){
    Space.indtot = (2*Space.Nmax + 1) * (2*Space.Nmax + 1) * (2*Space.Nmax + 1) * 2;
    Space.qmins.t = 1, Space.qsizes.t = 1;
  }
  else{ std::cerr << "No Protons or Neutrons Entered!!!" << std::endl; exit(1); }

  // allocate memory for quantum numbers for each state
  Space.qnums = new State[Space.indtot];
  Space.nmax = std::floor(std::sqrt(Space.Nmax)); // nmax^2 <= Nmax
  for(int shell = 0; shell <= Space.Nmax; ++shell){
    for(int nx = -Space.nmax; nx <= Space.nmax; ++nx){
      for(int ny = -Space.nmax; ny <= Space.nmax; ++ny){
        for(int nz = -Space.nmax; nz <= Space.nmax; ++nz){
          if(shell != nx*nx + ny*ny + nz*nz){ continue; }
          for(int sz = -1; sz <= 1; sz = sz+2){ // spin x2
            for( int tz = -1; tz <= 1; tz = tz+2){ // isospin x2
	      if(tz == -1){
                ++pcount;
                if(Parameters.Pshells == 0){ continue; }
                E = 4.0*(nx*nx + ny*ny + nz*nz);
                Space.qnums[count].energy = E;
		if(shell < Parameters.Pshells){ Space.qnums[count].type = "hole"; ++holcount; ++phcount; }
		else{ Space.qnums[count].type = "particle"; ++parcount; }
	      }
	      else if(tz == 1){
		++ncount;
		if(Parameters.Nshells == 0){ continue; }
		E = 4.0*(nx*nx + ny*ny + nz*nz);
		Space.qnums[count].energy = E;
		if(shell < Parameters.Nshells){ Space.qnums[count].type = "hole"; ++holcount; ++nhcount; }
		else{ Space.qnums[count].type = "particle"; ++parcount; }
	      }
	      Space.qnums[count].nx = nx;
	      Space.qnums[count].ny = ny;
	      Space.qnums[count].nz = nz;
	      Space.qnums[count].m = sz;
	      Space.qnums[count].t = tz;
	      count++;
	    }
	  }   
	} 
      }
    }
  }
  Space.indtot = count;
  Space.qsizes.m = 3; // -2, 0, 2
  Space.qsizes.nx = 4*Space.nmax + 1;
  Space.qsizes.ny = 4*Space.nmax + 1;
  Space.qsizes.nz = 4*Space.nmax + 1;
  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = nhcount;

  // With the number of hole states, the length scale and state energies can be calculated
  double L = pow(holcount/Parameters.density, 1.0/3.0);
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].t == -1){ Space.qnums[i].energy *= proton_prefac*M_PI*M_PI/(L*L); }
    else if(Space.qnums[i].t == 1){ Space.qnums[i].energy *= neutron_prefac*M_PI*M_PI/(L*L); }
  }

  // Change energies to Hartree-Fock energies, E_p = E_p + 2*V_pipi
  for(int p = 0; p < Space.indtot; ++p){
    for(int i = 0; i < Space.indhol; ++i){
      if(p == i){ continue; }
      Space.qnums[p].energy += 2*V_Minnesota(Space, p, i, p, i, L);
    }
  }

  // Order two-body states with respect to Nx, Ny, Nz for channel index function
  int count1 = 0;
  int N2max = 4*Space.Nmax;
  Space.map_2b = new int[Space.qsizes.nx * Space.qsizes.ny * Space.qsizes.nz];
  for(int nx = -2 * Space.nmax; nx <= 2 * Space.nmax; ++nx){
    for(int ny = -2 * Space.nmax; ny <= 2 * Space.nmax; ++ny){
      for(int nz = -2 * Space.nmax; nz <= 2 * Space.nmax; ++nz){
	if(nx*nx + ny*ny + nz*nz <= N2max){
	  Space.map_2b[(nx + 2*Space.nmax) * Space.qsizes.ny*Space.qsizes.nz + (ny + 2*Space.nmax) * Space.qsizes.nz + (nz + 2*Space.nmax)] = count1;
	  ++count1;
	}
      }
    }
  }
  Space.size_2b = count1 * Space.qsizes.t * Space.qsizes.m;
}

// Function that maps two-body quantum numbers onto a unique index
int Chan_2bInd(const Model_Space &Space, const State &State)
{
  return Space.map_2b[(State.nx + 2*Space.nmax) * Space.qsizes.ny*Space.qsizes.nz + (State.ny + 2*Space.nmax) * Space.qsizes.nz +
		      (State.nz + 2*Space.nmax)] * Space.qsizes.m*Space.qsizes.t + int((State.m + 2)/2) * Space.qsizes.t + 
                      int((State.t - 2*Space.qmins.t)/2);
}

int kron_del(const int &i, const int &j)
{
  if(i != j){ return 0; }
  return 1;
}

int SpinExchange(const int &i, const int &j, const int &k, const int &l)
{
  if(i == l && j == k){ return 1; }
  else{ return 0; }
}

double V_Minnesota(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L)
{
  double V_R1, V_T1, V_S1, V_R2, V_T2, V_S2;
  double V_0R, V_0T, V_0S;
  double kappa_R, kappa_T, kappa_S;
  double kX1, kY1, kZ1, kX2, kY2, kZ2;
  double qSquared1, spinEx1, isoSpinEx1, qSquared2, spinEx2, isoSpinEx2;
  double IsIt1, PsIt1, PsPt1, IsPt1, IsIt2, PsIt2, PsPt2, IsPt2;
  V_0R = 200; //MeV
  V_0T = 178; //MeV
  V_0S = 91.85; //MeV
  kappa_R = 1.487; //fm^-2
  kappa_T = 0.639; //fm^-2
  kappa_S = 0.465; //fm^-2

  if(Space.qnums[qi].nx + Space.qnums[qj].nx != Space.qnums[qk].nx + Space.qnums[ql].nx){ return 0.0; }
  if(Space.qnums[qi].ny + Space.qnums[qj].ny != Space.qnums[qk].ny + Space.qnums[ql].ny){ return 0.0; }
  if(Space.qnums[qi].nz + Space.qnums[qj].nz != Space.qnums[qk].nz + Space.qnums[ql].nz){ return 0.0; }
  if(Space.qnums[qi].m + Space.qnums[qj].m != Space.qnums[qk].m + Space.qnums[ql].m){ return 0.0; }
  if(Space.qnums[qi].t + Space.qnums[qj].t != Space.qnums[qk].t + Space.qnums[ql].t){ return 0.0; }

  kX1 = (M_PI/L) * (Space.qnums[qi].nx - Space.qnums[qj].nx - Space.qnums[qk].nx + Space.qnums[ql].nx);
  kY1 = (M_PI/L) * (Space.qnums[qi].ny - Space.qnums[qj].ny - Space.qnums[qk].ny + Space.qnums[ql].ny);
  kZ1 = (M_PI/L) * (Space.qnums[qi].nz - Space.qnums[qj].nz - Space.qnums[qk].nz + Space.qnums[ql].nz);

  kX2 = (M_PI/L) * (Space.qnums[qi].nx - Space.qnums[qj].nx - Space.qnums[ql].nx + Space.qnums[qk].nx);
  kY2 = (M_PI/L) * (Space.qnums[qi].ny - Space.qnums[qj].ny - Space.qnums[ql].ny + Space.qnums[qk].ny);
  kZ2 = (M_PI/L) * (Space.qnums[qi].nz - Space.qnums[qj].nz - Space.qnums[ql].nz + Space.qnums[qk].nz);

  qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
  qSquared2 = kX2 * kX2 + kY2 * kY2 + kZ2 * kZ2;
  
  V_R1 = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5) * exp(-qSquared1/(4*kappa_R));
  V_T1 = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5) * exp(-qSquared1/(4*kappa_T));
  V_S1 = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5) * exp(-qSquared1/(4*kappa_S));

  V_R2 = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5) * exp(-qSquared2/(4*kappa_R));
  V_T2 = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5) * exp(-qSquared2/(4*kappa_T));
  V_S2 = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5) * exp(-qSquared2/(4*kappa_S));
  
  spinEx1 = SpinExchange(Space.qnums[qi].m, Space.qnums[qj].m, Space.qnums[qk].m, Space.qnums[ql].m);
  isoSpinEx1 = SpinExchange(Space.qnums[qi].t, Space.qnums[qj].t, Space.qnums[qk].t, Space.qnums[ql].t);

  spinEx2 = SpinExchange(Space.qnums[qi].m, Space.qnums[qj].m, Space.qnums[ql].m, Space.qnums[qk].m);
  isoSpinEx2 = SpinExchange(Space.qnums[qi].t, Space.qnums[qj].t, Space.qnums[ql].t, Space.qnums[qk].t);
  
  IsIt1 = kron_del(Space.qnums[qi].m, Space.qnums[qk].m) * kron_del(Space.qnums[qj].m, Space.qnums[ql].m) * 
    kron_del(Space.qnums[qi].t, Space.qnums[qk].t) * kron_del(Space.qnums[qj].t, Space.qnums[ql].t);
  PsIt1 = spinEx1 * kron_del(Space.qnums[qi].t, Space.qnums[qk].t) * kron_del(Space.qnums[qj].t, Space.qnums[ql].t);
  PsPt1 = spinEx1 * isoSpinEx1;
  IsPt1 = kron_del(Space.qnums[qi].m, Space.qnums[qk].m)*kron_del(Space.qnums[qj].m, Space.qnums[ql].m) * isoSpinEx1;

  IsIt2 = kron_del(Space.qnums[qi].m, Space.qnums[ql].m) * kron_del(Space.qnums[qj].m, Space.qnums[qk].m) * 
    kron_del(Space.qnums[qi].t, Space.qnums[ql].t) * kron_del(Space.qnums[qj].t, Space.qnums[qk].t);
  PsIt2 = spinEx2 * kron_del(Space.qnums[qi].t, Space.qnums[ql].t) * kron_del(Space.qnums[qj].t, Space.qnums[qk].t);
  PsPt2 = spinEx2 * isoSpinEx2;
  IsPt2 = kron_del(Space.qnums[qi].m, Space.qnums[ql].m) * kron_del(Space.qnums[qj].m, Space.qnums[qk].m) * isoSpinEx2;

  return 0.5 * (V_R1 + 0.5*V_T1 + 0.5*V_S1) * IsIt1 + 
    0.25 * (V_T1 - V_S1) * PsIt1 - 
    0.5 * (V_R1 + 0.5*V_T1 + 0.5*V_S1) * PsPt1 - 
    0.25 * (V_T1 - V_S1) * IsPt1 -
    0.5 * (V_R2 + 0.5*V_T2 + 0.5*V_S2) * IsIt2 - 
    0.25 * (V_T2 - V_S2) * PsIt2 +
    0.5 * (V_R2 + 0.5*V_T2 + 0.5*V_S2) * PsPt2 + 
    0.25 * (V_T2 - V_S2) * IsPt2;
}
