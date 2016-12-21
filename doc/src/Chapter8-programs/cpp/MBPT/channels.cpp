#include "channels.hpp"

void Setup_Channels_MBPT(const Input_Parameters &Parameters, const Model_Space &Space, Channels &Chan)
{
  State state0;
  int idx;
  Chan.size = Space.size_2b; // size1 = number of direct channels

  Chan.hhvec = new int*[Chan.size];
  Chan.ppvec = new int*[Chan.size];
  Chan.nhh = new int[Chan.size];
  Chan.npp = new int[Chan.size];

  for(int chan = 0; chan < Chan.size; ++chan){
    Chan.nhh[chan] = 0;
    Chan.npp[chan] = 0;
  }

  // count number of hh and pp pairs in each channel
  for(int i = 0; i < Space.indhol; ++i){
    for(int j = 0; j < Space.indhol; ++j){
      if(i == j){ continue; }
      plus(state0, Space.qnums[i], Space.qnums[j]);
      idx = Chan_2bInd(Space, state0);
      ++Chan.nhh[idx];
    }
  }
  for(int a = Space.indhol; a < Space.indtot; ++a){
    for(int b = Space.indhol; b < Space.indtot; ++b){
      if(a == b){ continue; }
      plus(state0, Space.qnums[a], Space.qnums[b]);
      idx = Chan_2bInd(Space, state0);
      ++Chan.npp[idx];
    }
  }

  // allocate memory for hhvec and ppvec which list hh and pp pairs
  for(int chan = 0; chan < Chan.size; ++chan){
    Chan.hhvec[chan] = new int[2 * Chan.nhh[chan]];
    Chan.ppvec[chan] = new int[2 * Chan.npp[chan]];
    Chan.nhh[chan] = 0;
    Chan.npp[chan] = 0;
  }

  // fill hhvec and ppvec in each channel
  for(int i = 0; i < Space.indhol; ++i){
    for(int j = 0; j < Space.indhol; ++j){
      if(i == j){ continue; }
      plus(state0, Space.qnums[i], Space.qnums[j]);
      idx = Chan_2bInd(Space, state0);
      Chan.hhvec[idx][2 * Chan.nhh[idx]] = i;
      Chan.hhvec[idx][2 * Chan.nhh[idx] + 1] = j;
      ++Chan.nhh[idx];
    }
  }
  for(int a = Space.indhol; a < Space.indtot; ++a){
    for(int b = Space.indhol; b < Space.indtot; ++b){
      if(a == b){ continue; }
      plus(state0, Space.qnums[a], Space.qnums[b]);
      idx = Chan_2bInd(Space, state0);
      Chan.ppvec[idx][2 * Chan.npp[idx]] = a;
      Chan.ppvec[idx][2 * Chan.npp[idx] + 1] = b;
      ++Chan.npp[idx];
    }
  }
}

// Function that maps two-body quantum numbers onto a unique index
int Chan_2bInd(const Model_Space &Space, const State &State)
{
  return Space.map_2b[(State.nx + 2*Space.nmax) * Space.qsizes.ny*Space.qsizes.nz + (State.ny + 2*Space.nmax) * Space.qsizes.nz +
		      (State.nz + 2*Space.nmax)] * Space.qsizes.m*Space.qsizes.t + int((State.m + 2)/2) * Space.qsizes.t + 
                      int((State.t - 2*Space.qmins.t)/2);
}
