#ifndef TWOBODYPOTENTIAL_H
#define TWOBODYPOTENTIAL_H

#include "state.hpp"

class twobodyPotential {
public:
  virtual double get_element(State *qnums, int qi, int qj, int qk, int ql) = 0;
};
#endif
