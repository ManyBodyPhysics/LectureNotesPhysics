#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;
#include "arg.h"
#include "one_body.h"
#include "verbose.h"
#include "special_functions.h"
#include "error.h"
#include "global_job_parameter.h"
#include "constants.h"

OneBody::OneBody(OneBodyArg source_arg)
{

  source_type = source_arg.source_type;
  lambda1 = source_arg.lambda1; 
  lambda2 = source_arg.lambda2; 
  lambda3 = source_arg.lambda3; 
  fermion_size = 2*GJP.Vol();
  psi = (double *) malloc(fermion_size*sizeof(double));  

}

OneBody::~OneBody()
{
  free(psi);
}

void OneBody::Set()
{

  char* fname = "void OneBody::Set()";
  VRB.Func(fname);

  int x_sites = GJP.Xsites();
  int y_sites = GJP.Ysites();
  int z_sites = GJP.Zsites();
  int vol = GJP.Vol();

  if (source_type==SOURCE_TYPE_MOM) {

    //---- No need to do anything!
    return;

  }

  if (source_type==SOURCE_TYPE_SHO) {

    if ( GJP.APBCQuery() ) { ERR.NotImplemented(fname,"Must use PBCs when a potential is present."); }

    //---- This is the SHO state in momentum space

    int j;
    int n1_max = 10;
    int n2_max = 10;
    int n3_max = 10;
  
    if ( n1_max > x_sites ) {
      n1_max = x_sites;
      VRB.Warn(fname, "Reset maximum SHO level nx to %d", n1_max);
    }
    if ( n2_max > y_sites ) {
      n2_max = y_sites;
      VRB.Warn(fname, "Reset maximum SHO level ny to %d", n2_max);
    }
    if ( n3_max > z_sites ) {
      n3_max = z_sites;
      VRB.Warn(fname, "Reset maximum SHO level nz to %d", n2_max);
    }

    VRB.Warn(fname, "Maximum SHO levels (nx,ny,nz) are (%d,%d,%d)",n1_max,n2_max,n3_max);

    //---- Copy coeffs from psi to coeffs
    double coeffs[n1_max][n2_max][n3_max][2];
    for (int n1=0; n1<n1_max; n1++) {
      for (int n2=0; n2<n2_max; n2++) {
        for (int n3=0; n3<n3_max; n3++) {
          j = n3;
          j+= n2*GJP.Zsites();
          j+= n1*GJP.Ysites()*GJP.Zsites();
          j*=2;
          coeffs[n1][n2][n3][0] = psi[j];
          coeffs[n1][n2][n3][1] = psi[j+1];
        }
      }
    }
 
    //---- Now lets set the source
    double dx;
    double dy;
    double dz;
    int x;
    int y;
    int z;
    int Y; // Y = y+x*Y
    int Z; // Z = z + (y+x*Y)*Z
    double h1[n1_max];
    double h2[n2_max];
    double h3[n3_max];

    for(int i=0; i<fermion_size;i+=2) {

      //---- Determine coordinates
      Z = i/2;
      z = Z%z_sites;
      Y = Z/z_sites;
      y = Y%y_sites;
      x = Y/y_sites;

      //---- Momentum space wve function is centered around p=0 (PBCS)
      if (x>=x_sites/2) { x -= x_sites; }
      if (y>=y_sites/2) { y -= y_sites; }
      if (z>=z_sites/2) { z -= z_sites; }

      //---- For each set of integers, compute p_j = 2 pi k_j / L_j and multiply by (L_0)_j
      dx = TWOPI * x * lambda1 / (double) x_sites;
      dy = TWOPI * y * lambda2 / (double) y_sites;
      dz = TWOPI * z * lambda2 / (double) z_sites;

      //---- Compute Hermite polynomials
      HermiteH(h1, n1_max, dx);
      HermiteH(h2, n2_max, dy);
      HermiteH(h3, n3_max, dz);

      //---- Compute wave function as a linear combination of SHO basis wave functions
      double tmp;
      psi[i] = 0.0;
      psi[i+1] = 0.0; 
      for (int n1=0; n1<n1_max; n1++) {
        for (int n2=0; n2<n2_max; n2++) {
          for (int n3=0; n3<n3_max; n3++) {
            tmp = h1[n1]*h2[n2]*h3[n3];
            psi[i] += coeffs[n1][n2][n3][0] * tmp;
            psi[i+1] += coeffs[n1][n2][n3][1] * tmp; 
          }
        }
      }
      
      tmp = exp( -(dx*dx+dy*dy+dz*dz)/2.0 );
      tmp *= sqrt(lambda1*lambda2*lambda3);

      //---- Finally, scale result by (2pi/L)^3
      tmp *= sqrt( (TWOPI*TWOPI*TWOPI)/vol );

      psi[i] *= tmp;
      psi[i+1] *= tmp; 

    }


  }


}

void OneBody::Set(string vml_filename)
{

  char* fname = "void OneBody::Set(string)";
  VRB.Func(fname);

  for (int i=0; i<fermion_size; i++) {
    psi[i]=0.0;
  }

  ifstream file;
  file.open(vml_filename.c_str());
  if (!file) {
    ERR.FileR(fname,vml_filename.c_str());
  }

  int coord1;
  int coord2;
  int coord3;

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

    file >> psi[2*j]; 
    file >> psi[2*j+1];

    if ( file.eof() ) break; // Should be moved two lines up?--if EOF, then accessing random memory.
  }

  file.close();

  Set(); 

}

void OneBody::Set(int coord1, int coord2, int coord3)
{

  char* fname = "void OneBody::Set(int, int, int)";
  VRB.Func(fname);

  for (int i=0; i<fermion_size; i++) {
    psi[i]=0.0;
  }

  int j = coord3;
  j+= coord2*GJP.Zsites();
  j+= coord1*GJP.Ysites()*GJP.Zsites();

  psi[2*j] = 1.0; 

  Set(); 

}

void OneBody::Set(double* source) {

  char* fname = "void OneBody::Set(double*)";
  VRB.Func(fname);

  for (int i=0; i<fermion_size; i++) {
    psi[i] = source[i];
  }

}


complex<double> OneBody::Project(double* phi)
{

  double re = 0.0;
  double im = 0.0;

  for(int i=0; i<fermion_size; i+=2) {
    re += psi[i]*phi[i]+psi[i+1]*phi[i+1];
    im += psi[i]*phi[i+1]-psi[i+1]*phi[i];
  }
  return complex<double>(re,im);

}

double* OneBody::Get()
{
  return psi;
}
