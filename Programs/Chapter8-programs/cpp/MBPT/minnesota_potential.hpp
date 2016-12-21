#ifndef MINNESOTA_H
#define MINNESOTA_H
#include <cmath>
#include "constants.hpp"
#include "twobody-potential.hpp"
#include "modelspace.hpp"

class minnesotaPotential : public twobodyPotential {
public:
  minnesotaPotential(int numberOfParticles, double density);
  double get_element(State *qnums, int qi, int qj, int qk, int ql);
private:
  double L;
  static const double V_0R = 200; //MeV
  static const double V_0T = 178; //MeV
  static const double V_0S = 91.85; //MeV
  static const double kappa_R = 1.487; //fm^-2
  static const double kappa_T = 0.639; //fm^-2
  static const double kappa_S = 0.465; //fm^-2
  double VRfactor;
  double VTfactor;
  double VSfactor;
  double piOverL;
};

int kroneckerDelta(const int &i, const int &j);
int spinExchangeTerm(const int &i, const int &j, const int &k, const int &l);
void calculate_sp_energies(Input_Parameters &Parameters, Model_Space &Space );
double V_Minnesota(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L);
#endif
