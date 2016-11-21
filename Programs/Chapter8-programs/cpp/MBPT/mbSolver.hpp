#ifndef MBSOLVER_H
#define MBSOLVER_H

#include <string>

class mbSolver {
public:
  virtual double getEnergy() = 0;
  std::string getName() {
    return name;
  }
  std::string name;
};

#endif
