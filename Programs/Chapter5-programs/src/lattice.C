#include <iostream>
using namespace std;
#include "lattice.h"
#include "arg.h"
#include "random.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"

Lattice::Lattice(LatticeArg lattice_arg, Random* rng)
{

  rng_p = rng;

  field_type = lattice_arg.field_type;

  field_size = 2*GJP.Vol();
  field = (double *) malloc(field_size*sizeof(double)); //Allocate memory for field

  Refresh();

}

Lattice::~Lattice()
{
  free(field);
}


double* Lattice::Get()
{
  return field;
}

int Lattice::FieldSize()
{
  return field_size;
}

void Lattice::Refresh()
{

  char* fname = "void Lattice::Refresh()";
  VRB.Func(fname);


  if (field_type == FIELD_TYPE_GAUSS) {
    for(int i=0; i<field_size; i+=2) {
      field[i] = (*rng_p).Gauss(1.0); // Eventually specify width in Lattice Arg
      field[i+1] = 0.0;
    }
  }

  if (field_type == FIELD_TYPE_COMPLEXGAUSS) {
    for(int i=0; i<field_size; i+=2) {
      field[i] = (*rng_p).Gauss(1.0); // Eventually specify width in Lattice Arg
      field[i+1] = (*rng_p).Gauss(1.0); // Eventually specify width in Lattice Arg
    }
  }

  if (field_type == FIELD_TYPE_ZTWO) {
    for(int i=0; i<field_size; i+=2) {
      field[i] = 2.0*(*rng_p).Z(2)-1.0; // Eventually specify "p" in Lattice Arg
      field[i+1] = 0.0;
    }
  }

  if (field_type == FIELD_TYPE_ZTHREE) {
    int rand_num;
    for(int i=0; i<field_size; i+=2) {
      rand_num = (*rng_p).Z(3);
      if (rand_num==0) {
        field[i] = 1.0;
        field[i+1] = 0.0;
      } else if (rand_num==1) {
        field[i] = -0.5;
        field[i+1] = 0.8660254037844388;
      } else {
        field[i] = -0.5;
        field[i+1] = -0.8660254037844388;
      }
    }
  }

}
