// ran1.C
// An implementation of ran1 decribed in Numerical Recipes in C

// Edited excerpt from the numerical recipes source code:
// “Minimal” random number generator of Park and Miller with Bays-Durham shuffle and added safeguards.
// Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values).
// Call with idum a negative integer to initialize; thereafter, do not alter idum between successive deviates in a sequence.
// RNMX should approximate the largest floating value that is less than 1.

#include "arg.h"
#include "uniform_deviates.h"
#include "error.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

Ran1::Ran1()
{
  iv = (long *) malloc(NTAB*sizeof(long));  
  Init(1);
}

Ran1::~Ran1()
{
  free(iv);
}

void Ran1::Init(long seed)
{
  idum = seed;
  int j;
  long k;
  float temp;

  if (idum <= 0) idum=1; // Be sure to prevent idum = 0 and negative idum.
  for (j=NTAB+7;j>=0;j--) {
    k=idum/IQ;
    idum=IA*(idum-k*IQ)-IR*k;
    if (idum < 0) idum += IM;
    if (j < NTAB) iv[j] = idum;
  }
  iy=iv[0];

}

double Ran1::Run()
{
  int j;
  long k;
  float temp;

  k=idum/IQ;
  idum=IA*(idum-k*IQ)-IR*k;
  if (idum < 0) idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;

}

int Ran1::StateSize()
{
  return(3+NTAB);
}

void Ran1::GetState(long* state)
{
  state[0] = StateSize();
  state[1] = idum;
  state[2] = iy;
  for (int i=3; i<StateSize(); i++) {
    state[i] = iv[i-3];
  }
}

void Ran1::SetState(long* state)
{
  char* fname = "Ran1::SetState(long*)";
  if (state[0]!=StateSize()) { ERR.General(fname,"state size is inconsistent with generator type."); }

  idum = state[1];
  iy = state[2];
  for (int i=3; i<StateSize(); i++) {
    iv[i-3]=state[i];
  }


}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
