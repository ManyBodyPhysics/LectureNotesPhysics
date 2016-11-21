#include <complex>
#include "config.h"
using namespace std;

#ifndef INCLUDED_SYSFUNC
#define INCLUDED_SYSFUNC

namespace Comms {

  //extern bool initialized;
  extern bool initialized;

  void Initialize();
  void Finalize();
  void RaiseError();

  int Rank();

  int Size();

  void Sync();

  complex<double> GlobalSum(complex<double>);

}

#endif
