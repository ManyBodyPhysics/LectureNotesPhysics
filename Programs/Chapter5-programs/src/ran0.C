// ran0.C
// An implementation of ran0 decribed in Numerical Recipes in C
// The period of ran0 is 2^31 − 2

// Excerpt from the numerical recipes source code:
// “Minimal” random number generator of Park and Miller.
// Returns a uniform random deviate between 0.0 and 1.0.
// Set or reset idum to any integer value (except the
// unlikely value MASK) to initialize the sequence;
// idum must not be altered between calls for successive
// deviates in a sequence.

#include "arg.h"
#include "uniform_deviates.h"
#include "error.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

Ran0::Ran0()
{
  Init(1);
}

Ran0::~Ran0()
{
}

void Ran0::Init(long seed)
{
  char* fname = "Ran0::Init(long)";
  if (seed == MASK) { ERR.General(fname,"seed cannot equal to MASK."); }
  idum = seed;
}

double Ran0::Run()
{
  long k;
  float ans;
  idum ^= MASK; // XORing with MASK allows use of zero and other
  k=idum/IQ;  //   simple bit patterns for idum.
  idum=IA*(idum-k*IQ)-IR*k; // Compute idum=(IA*idum) % IM without over-
  if (idum < 0) idum += IM; //   flows by Schrage’s method.
  ans=AM*idum; //Convert idum to a floating result.
  idum ^= MASK;
  return ans; // Unmask before return.
}

int Ran0::StateSize()
{
  return(2);
}

void Ran0::GetState(long* state)
{
  state[0] = StateSize();
  state[1] = idum;
}

void Ran0::SetState(long* state)
{
  char* fname = "Ran0::SetState(long*)";
  if (state[0]!=StateSize()) { ERR.General(fname,"state size is inconsistent with generator type."); }
  idum = state[1];
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK
