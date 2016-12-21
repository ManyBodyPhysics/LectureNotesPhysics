#include "modelspace.hpp"

void Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space) {
  double E;
  
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
  if (Parameters.Shells > 40) {
    std::cerr << "Nmax too big!" << std::endl;
    exit(1);
  }
  Space.Nmax = shellnums[Parameters.Shells - 1];

  Parameters.numberOfProtons = 0;
  Parameters.numberOfNeutrons = 0;
  // Find maximum number of states (indtot) depending on Nmax and whether or not there are protons/neutrons
  if (Parameters.Pshells != 0 && Parameters.Nshells != 0) { 
    Space.indtot = (2*Space.Nmax + 1) * (2*Space.Nmax + 1) * (2*Space.Nmax + 1) * 2 * 2;
    Space.qmins.t = -1, Space.qsizes.t = 2;
  } else if(Parameters.Pshells != 0 && Parameters.Nshells == 0) {
    Space.indtot = (2*Space.Nmax + 1) * (2*Space.Nmax + 1) * (2*Space.Nmax + 1) * 2;
    Space.qmins.t = -1, Space.qsizes.t = 1;
  } else if(Parameters.Pshells == 0 && Parameters.Nshells != 0) {
    Space.indtot = (2*Space.Nmax + 1) * (2*Space.Nmax + 1) * (2*Space.Nmax + 1) * 2;
    Space.qmins.t = 1, Space.qsizes.t = 1;
  } else {
    std::cerr << "No Protons or Neutrons Entered!!!" << std::endl;
    exit(1);
  }

  // allocate memory for quantum numbers for each state
  Space.qnums = new State[Space.indtot];
  Space.nmax = std::floor(std::sqrt(Space.Nmax)); // nmax^2 <= Nmax
  for (int shell = 0; shell <= Space.Nmax; ++shell) {
    for (int nx = -Space.nmax; nx <= Space.nmax; ++nx) {
      for (int ny = -Space.nmax; ny <= Space.nmax; ++ny) {
        for (int nz = -Space.nmax; nz <= Space.nmax; ++nz) {
          if (shell != nx*nx + ny*ny + nz*nz) { continue; }
          for (int sz = -1; sz <= 1; sz = sz+2) { // spin x2
            for( int tz = -1; tz <= 1; tz = tz+2) { // isospin x2
              if (tz == -1) {
                if (Parameters.Pshells == 0) { continue; }
                ++pcount;
                E = 4.0*(nx*nx + ny*ny + nz*nz);
                Space.qnums[count].energy = E;
                if (shell < Parameters.Pshells) {
                  Space.qnums[count].type = "hole";
                  ++holcount;
                  ++phcount;
                  Parameters.numberOfProtons++;
                } else {
                  Space.qnums[count].type = "particle";
                  ++parcount;
                }
              } else if (tz == 1) {
                if (Parameters.Nshells == 0) { continue; }
                ++ncount;
                E = 4.0*(nx*nx + ny*ny + nz*nz);
                Space.qnums[count].energy = E;
                if (shell < Parameters.Nshells) {
                  Space.qnums[count].type = "hole";
                  ++holcount;
                  ++nhcount;
                  Parameters.numberOfNeutrons++;
                } else {
                  Space.qnums[count].type = "particle";
                  ++parcount;
                }
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
  Parameters.numberOfParticles = Parameters.numberOfProtons + Parameters.numberOfNeutrons;
}

