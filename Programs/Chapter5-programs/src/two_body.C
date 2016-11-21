#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;
#include "arg.h"
#include "two_body.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"
#include "constants.h"
#include "dispersion.h"
#include "special_functions.h"

TwoBody::TwoBody(TwoBodyArg two_body_arg)
{
  char* fname = "TwoBody::TwoBody(TwoBodyArg)";

  vol = GJP.Vol();
  wavefunc = (double *) malloc(vol*sizeof(double));
  opposite_parity_index= (int *) malloc(vol*sizeof(int));

  //---- Compute opposite parity index
  int x_sites = GJP.Xsites();
  int y_sites = GJP.Ysites();
  int z_sites = GJP.Zsites();
  int i;

  int bc[3];

  switch (GJP.Xbc()) {
    case BOUNDARY_TYPE_PRD:
      bc[0] = 0;
      break;  
    case BOUNDARY_TYPE_APRD:
      bc[0] = 1; 
      break;  
    default:
      ERR.NotImplemented(fname, "Boundary condition type not supported");
  }

  switch (GJP.Ybc()) {
    case BOUNDARY_TYPE_PRD:
      bc[1] = 0;
      break;  
    case BOUNDARY_TYPE_APRD:
      bc[1] = 1; 
      break;  
    default:
      ERR.NotImplemented(fname, "Boundary condition type not supported");
  }

  switch (GJP.Zbc()) {
    case BOUNDARY_TYPE_PRD:
      bc[2] = 0;
      break;  
    case BOUNDARY_TYPE_APRD:
      bc[2] = 1; 
      break;  
    default:
      ERR.NotImplemented(fname, "Boundary condition type not supported");
  }

  for (int x=0; x<x_sites; x++) {
    for (int y=0; y<y_sites; y++) {
      for (int z=0; z<z_sites; z++) {

        i =  z;
        i += y*z_sites;
        i += x*y_sites*z_sites;

        opposite_parity_index[i] =  (z_sites-z-bc[2])%z_sites;
        opposite_parity_index[i] += ( (y_sites-y-bc[1])%y_sites )*z_sites;
        opposite_parity_index[i] += ( (x_sites-x-bc[0])%x_sites )*y_sites*z_sites;
        opposite_parity_index[i] *= 2;

      }
    }
  }

  //---- Construct two body wave function
  if (two_body_arg.wavefunc_type==WAVEFUNC_TYPE_NONE) {
    for (int i=0; i<vol; i++) {
      wavefunc[i] = 0.0;
    }
  }

  if (two_body_arg.wavefunc_type==WAVEFUNC_TYPE_UNIFORM) {
    for (int i=0; i<vol; i++) {
      wavefunc[i] = 1.0/vol;
    }
  }

  if (two_body_arg.wavefunc_type==WAVEFUNC_TYPE_GND) {

    Dispersion dispersion1(two_body_arg.dispersion_type1, two_body_arg.mass1, CUTOFF_TYPE_HARD);
    Dispersion dispersion2(two_body_arg.dispersion_type2, two_body_arg.mass2, CUTOFF_TYPE_HARD);

    double xi1;
    double xi2;
    double lambda = two_body_arg.lambda;

    for (int i=0; i<vol; i++) {
      xi1 = 1.0 + dispersion1.Get(i);
      xi2 = 1.0 + dispersion2.Get(i);
 
      // Here, I use the fermion wave function, as derived in my note.
      // Note that psi(p) = xi(p)/(lambda*xi(p)^2-1) is the eigenstate of the transfer matrix.
      // The correlation function, however, is given by <final| D^{-1/2} T^N D^{-1/2} |initial>
      // Hence <p|final> = xi(p) psi(p), where xi is the same as D^{1/2} in momentum space).
      //
      // wavefunc[i] = xi^2/(lambda*xi^2-1.0);
      wavefunc[i] = 1.0;                     // This way avoids NAN when dispersion type is QUADRATIC, or PERFECT
      wavefunc[i] /= lambda - 1.0/(xi1*xi2); // so do it this way instead:
    }
  }

  if (two_body_arg.wavefunc_type==WAVEFUNC_TYPE_TGND) {

    Dispersion dispersion(DISPERSION_TYPE_QUADRATIC, 1.0, CUTOFF_TYPE_HARD);
    double b = two_body_arg.lambda;
    double psq;
    double xi;

    for (int i=0; i<vol; i++) {
      psq = 2.0*dispersion.Get(i); 

      if ( psq > GJP.CutoffSq() ) {

        wavefunc[i] = 0.0;

      } else {

        xi = sqrt(psq)/(2.0*b);

        if (xi < 1e-15) { // Wave function is non-singular at p=0, but need to treat this case separately
          wavefunc[i] = 1.0;
        } else {
          wavefunc[i] = DawsonF(xi)/xi;
        }

      }

    }

  }

  if (two_body_arg.wavefunc_type==WAVEFUNC_TYPE_PAIR1) {

    Dispersion dispersion(DISPERSION_TYPE_QUADRATIC, 1.0, CUTOFF_TYPE_HARD);
    double bsq = two_body_arg.lambda;
    double psq;

    for (int i=0; i<vol; i++) {
      psq = 2.0*dispersion.Get(i); 
      wavefunc[i] = 1.0/(psq+bsq);
    }

  }

  if (two_body_arg.wavefunc_type==WAVEFUNC_TYPE_PAIR2) {

    Dispersion dispersion(DISPERSION_TYPE_QUADRATIC, 1.0, CUTOFF_TYPE_HARD);
    double b = two_body_arg.lambda;
    double psq;

    for (int i=0; i<vol; i++) {
      psq = 2.0*dispersion.Get(i); 
      if (psq < 1e-15) {
        wavefunc[i] = 0.0; // Omit divergent contribution to wave funcion--must be treated separately
      } else {
        wavefunc[i] = exp(- b*sqrt(psq) )/psq;
      }
    }

  }


}

TwoBody::~TwoBody()
{
  free(wavefunc);
  free(opposite_parity_index);
}

complex<double> TwoBody::Run(double* prop1, double* prop2)
{

//  char* fname = "complex<double> TwoBody::Run(double*, double*)";
//  VRB.Func(fname);

  // Compute sum_p psi(p) G(p) G(-p)
  int j;
  int k;
  double re=0.0;
  double im=0.0;


  for (int i=0; i<vol; i++) {

//    cout << i << " " << wavefunc[i] << endl;

    j = 2*i;
    k = opposite_parity_index[i];

    re += (prop1[j]*prop2[k] - prop1[j+1]*prop2[k+1]) * wavefunc[i];
    im += (prop1[j]*prop2[k+1] + prop1[j+1]*prop2[k]) * wavefunc[i];

  }

  return complex<double>(re,im);

}

void TwoBody::Deform(string vml_filename)
{

  char* fname = "void TwoBody::Deform(string)";
  VRB.Func(fname);

  ifstream file;
  file.open(vml_filename.c_str());
  if (!file) {
    ERR.FileR(fname,vml_filename.c_str());
  }

  int coord1;
  int coord2;
  int coord3;

  double re;

  int j;
  //---- Need to add error checking, etc....
  while (true) {

    file >> coord1;
    if (coord1<0) coord1 += GJP.Xsites();
    file >> coord2;
    if (coord2<0) coord2 += GJP.Ysites();
    file >> coord3;
    if (coord3<0) coord3 += GJP.Zsites();

    j = coord3;
    j+= coord2*GJP.Zsites();
    j+= coord1*GJP.Ysites()*GJP.Zsites();

    file >> re; 

    VRB.Flow(fname,"Appending %f to %d component of wave function.", re, j);
    wavefunc[j] += re;

    if ( file.eof() ) break; // Should be moved two lines up--if EOF, then accessing random memory.

  }

  file.close();

}


