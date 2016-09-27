#include "minnesota_potential.hpp"

minnesotaPotential::minnesotaPotential(int numberOfParticles, double density) {

  L = pow((numberOfParticles)/density, 1./3.);

  piOverL = M_PI/L;
  VRfactor = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5);
  VTfactor = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5);
  VSfactor = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5);
}


int kroneckerDelta(const int &i, const int &j)
{
  if(i != j){ return 0; }
  return 1;
}

int spinExchangeTerm(const int &i, const int &j, const int &k, const int &l)
{
  if(i == l && j == k){ return 1; }
  else{ return 0; }
}

inline double minnesotaPotential::get_element(State *qnums, int qi, int qj, int qk, int ql) {
  double V_R1, V_T1, V_S1, V_R2, V_T2, V_S2;
  double kX1, kY1, kZ1, kX2, kY2, kZ2;
  double qSquared1, spinEx1, isoSpinEx1, qSquared2, spinEx2, isoSpinEx2;
  double IsIt1, PsIt1, PsPt1, IsPt1, IsIt2, PsIt2, PsPt2, IsPt2;

  if(qnums[qi].nx + qnums[qj].nx != qnums[qk].nx + qnums[ql].nx){ return 0.0; }
  if(qnums[qi].ny + qnums[qj].ny != qnums[qk].ny + qnums[ql].ny){ return 0.0; }
  if(qnums[qi].nz + qnums[qj].nz != qnums[qk].nz + qnums[ql].nz){ return 0.0; }
  if(qnums[qi].m + qnums[qj].m != qnums[qk].m + qnums[ql].m){ return 0.0; }
  if(qnums[qi].t + qnums[qj].t != qnums[qk].t + qnums[ql].t){ return 0.0; }

  kX1 = piOverL * (qnums[qi].nx - qnums[qj].nx - qnums[qk].nx + qnums[ql].nx);
  kY1 = piOverL * (qnums[qi].ny - qnums[qj].ny - qnums[qk].ny + qnums[ql].ny);
  kZ1 = piOverL * (qnums[qi].nz - qnums[qj].nz - qnums[qk].nz + qnums[ql].nz);

  kX2 = piOverL * (qnums[qi].nx - qnums[qj].nx - qnums[ql].nx + qnums[qk].nx);
  kY2 = piOverL * (qnums[qi].ny - qnums[qj].ny - qnums[ql].ny + qnums[qk].ny);
  kZ2 = piOverL * (qnums[qi].nz - qnums[qj].nz - qnums[ql].nz + qnums[qk].nz);

  qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
  qSquared2 = kX2 * kX2 + kY2 * kY2 + kZ2 * kZ2;
  
  V_R1 = VRfactor * exp(-qSquared1/(4*kappa_R));
  V_T1 = VTfactor * exp(-qSquared1/(4*kappa_T));
  V_S1 = VSfactor * exp(-qSquared1/(4*kappa_S));

  V_R2 = VRfactor * exp(-qSquared2/(4*kappa_R));
  V_T2 = VTfactor * exp(-qSquared2/(4*kappa_T));
  V_S2 = VSfactor * exp(-qSquared2/(4*kappa_S));
  
  spinEx1 = spinExchangeTerm(qnums[qi].m, qnums[qj].m, qnums[qk].m, qnums[ql].m);
  isoSpinEx1 = spinExchangeTerm(qnums[qi].t, qnums[qj].t, qnums[qk].t, qnums[ql].t);

  spinEx2 = spinExchangeTerm(qnums[qi].m, qnums[qj].m, qnums[ql].m, qnums[qk].m);
  isoSpinEx2 = spinExchangeTerm(qnums[qi].t, qnums[qj].t, qnums[ql].t, qnums[qk].t);
  
  IsIt1 = kroneckerDelta(qnums[qi].m, qnums[qk].m) * kroneckerDelta(qnums[qj].m, qnums[ql].m) * 
    kroneckerDelta(qnums[qi].t, qnums[qk].t) * kroneckerDelta(qnums[qj].t, qnums[ql].t);
  PsIt1 = spinEx1 * kroneckerDelta(qnums[qi].t, qnums[qk].t) * kroneckerDelta(qnums[qj].t, qnums[ql].t);
  PsPt1 = spinEx1 * isoSpinEx1;
  IsPt1 = kroneckerDelta(qnums[qi].m, qnums[qk].m)*kroneckerDelta(qnums[qj].m, qnums[ql].m) * isoSpinEx1;

  IsIt2 = kroneckerDelta(qnums[qi].m, qnums[ql].m) * kroneckerDelta(qnums[qj].m, qnums[qk].m) * 
    kroneckerDelta(qnums[qi].t, qnums[ql].t) * kroneckerDelta(qnums[qj].t, qnums[qk].t);
  PsIt2 = spinEx2 * kroneckerDelta(qnums[qi].t, qnums[ql].t) * kroneckerDelta(qnums[qj].t, qnums[qk].t);
  PsPt2 = spinEx2 * isoSpinEx2;
  IsPt2 = kroneckerDelta(qnums[qi].m, qnums[ql].m) * kroneckerDelta(qnums[qj].m, qnums[qk].m) * isoSpinEx2;

  return 0.5 * (V_R1 + 0.5*V_T1 + 0.5*V_S1) * IsIt1 +
    0.25 * (V_T1 - V_S1) * PsIt1 -
    0.5 * (V_R1 + 0.5*V_T1 + 0.5*V_S1) * PsPt1 -
    0.25 * (V_T1 - V_S1) * IsPt1 -
    0.5 * (V_R2 + 0.5*V_T2 + 0.5*V_S2) * IsIt2 -
    0.25 * (V_T2 - V_S2) * PsIt2 +
    0.5 * (V_R2 + 0.5*V_T2 + 0.5*V_S2) * PsPt2 +
    0.25 * (V_T2 - V_S2) * IsPt2;
}

void calculate_sp_energies(Input_Parameters &Parameters, Model_Space &Space ) {
  // With the number of hole states, the length scale and state energies can be calculated
  double L = pow(Space.indhol/Parameters.density, 1.0/3.0);
  for (int i = 0; i < Space.indtot; ++i) {
    if (Space.qnums[i].t == -1) { 
      Space.qnums[i].energy *= proton_prefac*M_PI*M_PI/(L*L);
    }
    else if (Space.qnums[i].t == 1) {
      Space.qnums[i].energy *= neutron_prefac*M_PI*M_PI/(L*L);
    }
  }

  // Change energies to Hartree-Fock energies, E_p = E_p + 2*V_pipi
  for (int p = 0; p < Space.indtot; ++p) {
    for (int i = 0; i < Space.indhol; ++i) {
      if (p == i) { continue; }
      Space.qnums[p].energy += 2*V_Minnesota(Space, p, i, p, i, L);
    }
  }

  // Order two-body states with respect to Nx, Ny, Nz for channel index function
  int count1 = 0;
  int N2max = 4*Space.Nmax;
  Space.map_2b = new int[Space.qsizes.nx * Space.qsizes.ny * Space.qsizes.nz];
  for (int nx = -2 * Space.nmax; nx <= 2 * Space.nmax; ++nx) {
    for (int ny = -2 * Space.nmax; ny <= 2 * Space.nmax; ++ny) {
      for (int nz = -2 * Space.nmax; nz <= 2 * Space.nmax; ++nz) {
        if (nx*nx + ny*ny + nz*nz <= N2max) {
          int idx = (nx + 2*Space.nmax) * Space.qsizes.ny*Space.qsizes.nz +
                                    (ny + 2*Space.nmax) * Space.qsizes.nz +
                                    (nz + 2*Space.nmax);
          Space.map_2b[idx] = count1;
          ++count1;
        }
      }
    }
  }
  Space.size_2b = count1 * Space.qsizes.t * Space.qsizes.m;

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
  
  spinEx1 = spinExchangeTerm(Space.qnums[qi].m, Space.qnums[qj].m, Space.qnums[qk].m, Space.qnums[ql].m);
  isoSpinEx1 = spinExchangeTerm(Space.qnums[qi].t, Space.qnums[qj].t, Space.qnums[qk].t, Space.qnums[ql].t);

  spinEx2 = spinExchangeTerm(Space.qnums[qi].m, Space.qnums[qj].m, Space.qnums[ql].m, Space.qnums[qk].m);
  isoSpinEx2 = spinExchangeTerm(Space.qnums[qi].t, Space.qnums[qj].t, Space.qnums[ql].t, Space.qnums[qk].t);
  
  IsIt1 = kroneckerDelta(Space.qnums[qi].m, Space.qnums[qk].m) * kroneckerDelta(Space.qnums[qj].m, Space.qnums[ql].m) * 
    kroneckerDelta(Space.qnums[qi].t, Space.qnums[qk].t) * kroneckerDelta(Space.qnums[qj].t, Space.qnums[ql].t);
  PsIt1 = spinEx1 * kroneckerDelta(Space.qnums[qi].t, Space.qnums[qk].t) * kroneckerDelta(Space.qnums[qj].t, Space.qnums[ql].t);
  PsPt1 = spinEx1 * isoSpinEx1;
  IsPt1 = kroneckerDelta(Space.qnums[qi].m, Space.qnums[qk].m)*kroneckerDelta(Space.qnums[qj].m, Space.qnums[ql].m) * isoSpinEx1;

  IsIt2 = kroneckerDelta(Space.qnums[qi].m, Space.qnums[ql].m) * kroneckerDelta(Space.qnums[qj].m, Space.qnums[qk].m) * 
    kroneckerDelta(Space.qnums[qi].t, Space.qnums[ql].t) * kroneckerDelta(Space.qnums[qj].t, Space.qnums[qk].t);
  PsIt2 = spinEx2 * kroneckerDelta(Space.qnums[qi].t, Space.qnums[ql].t) * kroneckerDelta(Space.qnums[qj].t, Space.qnums[qk].t);
  PsPt2 = spinEx2 * isoSpinEx2;
  IsPt2 = kroneckerDelta(Space.qnums[qi].m, Space.qnums[ql].m) * kroneckerDelta(Space.qnums[qj].m, Space.qnums[qk].m) * isoSpinEx2;

  return 0.5 * (V_R1 + 0.5*V_T1 + 0.5*V_S1) * IsIt1 + 
    0.25 * (V_T1 - V_S1) * PsIt1 - 
    0.5 * (V_R1 + 0.5*V_T1 + 0.5*V_S1) * PsPt1 - 
    0.25 * (V_T1 - V_S1) * IsPt1 -
    0.5 * (V_R2 + 0.5*V_T2 + 0.5*V_S2) * IsIt2 - 
    0.25 * (V_T2 - V_S2) * PsIt2 +
    0.5 * (V_R2 + 0.5*V_T2 + 0.5*V_S2) * PsPt2 + 
    0.25 * (V_T2 - V_S2) * IsPt2;
}
