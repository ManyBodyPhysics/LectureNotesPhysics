// ran3.C
// An implementation of ran3 decribed in Numerical Recipes in C

// Edited excerpt from the numerical recipes source code:
// Returns a uniform random deviate between 0.0 and 1.0.
// Set idum to any positive value to initialize or reinitialize the sequence.

#include "arg.h"
#include "uniform_deviates.h"
#include "error.h"

#define MBIG 1000000000  // According to Knuth, any large MBIG, and any smaller
#define MSEED 161803398  // (but still large) MSEED can be substituted for the above values.
#define MZ 0
#define FAC (1.0/MBIG)


Ran3::Ran3()
{
  ma = (long *) malloc(56*sizeof(long)); // The value 56 (range ma[1..55]) is special and
  Init(1); // default
}

Ran3::~Ran3()
{
  free(ma);
}

void Ran3::Init(long seed)
{
 
  long mj,mk;
  int i,ii,k;

  mj=MSEED-(seed < 0 ? -seed : seed); // Initialize ma[55] using the seed
  mj %= MBIG;                            //   idum and the large number MSEED.
  ma[55]=mj;
  mk=1;
  for (i=1;i<55;i++) { // Now initialize the rest of the table,
    ii=(21*i) % 55;     //   in a slightly random order,
    ma[ii]=mk;          //   with numbers that are not especially random.
    mk=mj-mk;
    if (mk < MZ) mk += MBIG;
    mj=ma[ii];
  }
  for (k=0;k<4;k++)
    for (i=1;i<56;i++) {
      ma[i] -= ma[1+(i+30) % 55];
      if (ma[i] < MZ) ma[i] += MBIG;
  }
  inext=0; // Prepare indices for our first generated number.
  inextp=31; //The constant 31 is special; see Knuth.

}

double Ran3::Run()
{
  long mj;

  // Here is where we start, except on initialization.
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;

}


int Ran3::StateSize()
{
  return(59);
}

void Ran3::GetState(long* state)
{
  state[0] = StateSize();
  state[1] = inext;
  state[2] = inextp;
  for (int i=3; i<StateSize(); i++) {
    state[i] = ma[i-3];
  }
}

void Ran3::SetState(long* state)
{
  char* fname = "Ran3::SetState(long*)";
  if (state[0]!=StateSize()) { ERR.General(fname,"state size is inconsistent with generator type."); }

  inext = state[1];
  inextp = state[2];
  for (int i=3; i<StateSize(); i++) {
    ma[i-3]=state[i];
  }
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
