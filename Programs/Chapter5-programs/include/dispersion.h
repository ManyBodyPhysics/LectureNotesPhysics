using namespace std;
#include "enum.h"

#ifndef INCLUDED_DISPERSION
#define INCLUDED_DISPERSION
class Dispersion
{
  private:
    double* dispersion;

  public:
    Dispersion(DispersionType, double, CutoffType);
    ~Dispersion();
    double Get(int, int, int); // Returns energy associated with momenta specifies by three integers
    double Get(int); // Returns energy associated with momenta specifies by collective index
};

#endif
