#include <fftw3.h> 
#include "phase.h"

#ifndef INCLUDED_FOURIER
#define INCLUDED_FOURIER

namespace Fourier {

  extern int instances;


  extern fftw_complex* in;
  extern fftw_complex* out;

  extern Phase phase;

  void Initialize();
  void Finalize();

  void Forward();
  void Backward();

}

#endif
