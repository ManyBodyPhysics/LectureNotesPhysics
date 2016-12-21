// ran2.C
// An implementation of ran2 decribed in Numerical Recipes in C

// Edited excerpt from the numerical recipes source code:
// Long period (> 2 × 1018) random number generator of L’Ecuyer with Bays-Durham shuffle and
// added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
// the endpoint values). Call with idum a positive integer to initialize; thereafter, do not alter
// idum between successive deviates in a sequence. RNMX should approximate the largest floating
// value that is less than 1.

#include "arg.h"
#include "uniform_deviates.h"
#include "error.h"

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

Ran2::Ran2()
{
  iv = (long *) malloc(NTAB*sizeof(long));  
  Init(1);
}

Ran2::~Ran2()
{
  free(iv);
}

void Ran2::Init(long seed)
{
  idum = seed;
  idum2=123456789;
  iy=0;

  int j;
  long k;

  if (idum <= 0) idum=1; // Be sure to prevent idum = 0 and negative idum.
  idum2=idum;
  for (j=NTAB+7;j>=0;j--) { // Load the shuffle table (after 8 warm-ups).
    k=idum/IQ1;
    idum=IA1*(idum-k*IQ1)-k*IR1;
    if (idum < 0) idum += IM1;
    if (j < NTAB) iv[j] = idum;
  }   
  iy=iv[0];

}

double Ran2::Run()
{
  int j;
  long k;
  float temp;

  k=idum/IQ1; // Start here when not initializing.
  idum=IA1*(idum-k*IQ1)-k*IR1; // Compute idum=(IA1*idum) % IM1 without
  if (idum < 0) idum += IM1;   // overflows by Schrage’s method.
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; // Compute idum2=(IA2*idum) % IM2 likewise.
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;      // Will be in the range 0..NTAB-1.
  iy=iv[j]-idum2; // Here idum is shuffled, idum and idum2 are
  iv[j] = idum;  // combined to generate output.
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; // Because users don’t expect endpoint values.
  else return temp;
}


int Ran2::StateSize()
{
  return(4+NTAB);
}

void Ran2::GetState(long* state)
{
  state[0] = StateSize();
  state[1] = idum;
  state[2] = idum2;
  state[3] = iy;
  for (int i=4; i<StateSize(); i++) {
    state[i] = iv[i-4];
  }
}

void Ran2::SetState(long* state)
{

  char* fname = "Ran2::SetState(long*)";
  if (state[0]!=StateSize()) { ERR.General(fname,"state size is inconsistent with generator type."); }

  idum = state[1];
  idum2 = state[2];
  iy = state[3];
  for (int i=4; i<StateSize(); i++) {
    iv[i-4]=state[i];
  }
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
