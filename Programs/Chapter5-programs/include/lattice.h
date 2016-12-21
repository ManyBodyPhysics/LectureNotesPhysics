#include "arg.h"
#include "random.h"

#ifndef INCLUDED_LATTICE
#define INCLUDED_LATTICE
class Lattice
{
  private:
    Random* rng_p;
    enum FieldType field_type;
    int field_size;
    double* field;

  public:
    Lattice(LatticeArg, Random*);
    ~Lattice();
    void Refresh();
    double* Get();
    int FieldSize();
    
};
#endif
