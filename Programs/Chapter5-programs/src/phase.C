#include <iostream>
#include <math.h>
using namespace std;
#include "phase.h"
#include "global_job_parameter.h"
#include "error.h"
#include "verbose.h"
#include "constants.h"



Phase::Phase(double* p)
{

  char* fname = "Phase::Phase(double*)";
  VRB.Func(fname);

  ERR.NotImplemented(fname, "Sorry.");

  PhaseGen(p);

}

Phase::Phase()
{

  char* fname = "Phase::Phase()";
  VRB.Func(fname);


  //---- Set the boundary conditions
  double bc[3];
  GJP.GetBcShift(bc);

  PhaseGen(bc);

}

Phase::~Phase()
{

  free(phase);

}


void Phase::PhaseGen(double* bc)
{

  vol = GJP.Vol();
  phase = (double*) malloc(2*vol*sizeof(double));

  //---- Set the dimensions
  int dim[3];
  dim[0] =  GJP.Xsites();
  dim[1] =  GJP.Ysites();
  dim[2] =  GJP.Zsites();

  //---- Compute phases for BCs 
  double w;
  int x;
  int y;
  int Y;
  int z;
  for (int i=0; i<vol; i++) {
    z = i%dim[2];
    Y = i/dim[2];
    y = Y%dim[1];
    x = Y/dim[1];

    w  = bc[0]*x/(double)dim[0];
    w += bc[1]*y/(double)dim[1];
    w += bc[2]*z/(double)dim[2];
    w *= TWOPI;

    phase[2*i]   = cos(w);
    phase[2*i+1] = sin(w);
  }


}

void Phase::Add(double* in)
{
// Replaces in[n] -> in[n]*exp(i w . n)

  double re;
  double im;

  for (int i=0; i<2*vol; i+=2) {

    re = in[i]*phase[i] - in[i+1]*phase[i+1];
    im = in[i]*phase[i+1] + in[i+1]*phase[i];

    in[i] = re; 
    in[i+1] = im; 

  }

}


void Phase::Subtract(double* in)
{
// Replaces in[n] -> in[n]*exp(-i w . n)

  double re;
  double im;

  for (int i=0; i<2*vol; i+=2) {

    re = in[i]*phase[i] + in[i+1]*phase[i+1];
    im = -in[i]*phase[i+1] + in[i+1]*phase[i];

    in[i] = re; 
    in[i+1] = im; 

  }

}

void Phase::Add(fftw_complex* in)
{
// Replaces in[n] -> in[n]*exp(i w . n)

  double re;
  double im;

  int j;
  for (int i=0; i<vol; i++) {

    j = 2*i;
    re = in[i][0]*phase[j] - in[i][1]*phase[j+1];
    im = in[i][0]*phase[j+1] + in[i][1]*phase[j];

    in[i][0] = re; 
    in[i][1] = im; 

  }

}


void Phase::Subtract(fftw_complex* in)
{
// Replaces in[n] -> in[n]*exp(-i w . n)

  double re;
  double im;

  int j;
  for (int i=0; i<vol; i++) {

    j = 2*i;
    re = in[i][0]*phase[j] + in[i][1]*phase[j+1];
    im = -in[i][0]*phase[j+1] + in[i][1]*phase[j];

    in[i][0] = re; 
    in[i][1] = im; 

  }

}
