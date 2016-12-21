#ifndef STATE_H
#define STATE_H
#include <string>

typedef struct State {
  int t; // x2
  int m; // x2
  int nx;
  int ny;
  int nz;
  double energy;
  std::string type;
} State;

void plus(State &S, const State &S1, const State &S2);
#endif
