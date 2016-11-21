#include "arg.h"
#include "uniform_deviates.h"

#ifndef INCLUDED_RANDOM
#define INCLUDED_RANDOM
class Random
// Random number generator base class
{
  private:
    UniformDeviate* r;
    RandomType random_type;
  public:
    Random(RandomArg);
    ~Random();
    double Uniform();
    double Gauss(double);
    int Z(int);
    void ReadState(const char*);
    void WriteState(const char*);
};

#endif
