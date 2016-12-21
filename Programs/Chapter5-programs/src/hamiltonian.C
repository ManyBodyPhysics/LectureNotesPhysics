#include <iostream>
#include <math.h>
#include <stdlib.h>
using namespace std;
#include "lattice.h"
#include "arg.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"
#include "hamiltonian.h"
#include "dispersion.h"
#include "constants.h"
#include "fourier.h"

Hamiltonian::Hamiltonian()
{

  char* fname = "void Hamiltonian::Hamiltonian()";
  VRB.Func(fname);

  Fourier::Initialize();

}

Hamiltonian::~Hamiltonian()
{

  char* fname = "void Hamiltonian::~Hamiltonian()";
  VRB.Func(fname);

  Fourier::Finalize();

}

Interactions::Interactions(  vector<Hamiltonian*>  interaction_list)
{

  char* fname = "Interactions::Interactions(  vector<Hamiltonian*> )";
  VRB.Func(fname);

  vol = GJP.Vol();
  a = (double *) malloc(2*vol*sizeof(double));  
  c = (double *) malloc(2*vol*sizeof(double));  

  interactions = interaction_list;

}

Interactions::~Interactions()
{
  free(a);
  free(c);
}

void Interactions::Evolve(double* b)
{

  for (int j=0; j<2*vol; j++) { c[j] = 0.0; }
  
  for (int i=0; i<interactions.size(); i++) {
    for (int j=0; j<2*vol; j++) { a[j] = b[j]; }
    interactions[i]->Evolve(a);
    for (int j=0; j<2*vol; j++) { c[j] += a[j]; }
  }

  for (int j=0; j<2*vol; j++) { b[j] *= 1.0-interactions.size(); }

  for (int j=0; j<2*vol; j++) { b[j] += c[j]; }

}

Interaction::Interaction(Lattice* lat, InteractionArg interaction_arg, KineticArg kinetic_arg1, KineticArg kinetic_arg2, double dt)
{

  char* fname = "Interaction::Interaction(Lattice*, InteractionArg, KineticArg, KineticArg, double)";
  VRB.Func(fname);

  lattice = lat;

  //---- Compute the interactions in momentum space
  int vol = GJP.Vol();
  interaction = (double *) malloc(vol*sizeof(double));  

  Dispersion dispersion1(kinetic_arg1.dispersion_type, kinetic_arg1.mass, CUTOFF_TYPE_NONE);
  Dispersion dispersion2(kinetic_arg2.dispersion_type, kinetic_arg2.mass, CUTOFF_TYPE_NONE);

  double xi1;
  double xi2;

  double psq;
  double mass = 2.0 * kinetic_arg1.mass * kinetic_arg2.mass; // reduced mass
  mass /= kinetic_arg1.mass + kinetic_arg2.mass; // reduced mass

  for(int i=0; i<vol; i++) {

    xi1 = 1.0 + dispersion1.Get(i);
    xi2 = 1.0 + dispersion2.Get(i);

    switch (interaction_arg.interaction_type) {
      case INTERACTION_TYPE_NONE:
        ERR.NotImplemented(fname,"Interaction type INTERACTION_TYPE_NONE.");
        break;
      case INTERACTION_TYPE_ONEMINUSXIINVSQ:
        psq = 1.0 - 1.0 /(xi1*xi2); 
        break;
      case INTERACTION_TYPE_XISQMINUSONE:
        psq = xi1*xi2 - 1.0; 
        break;
      default:
        ERR.NotImplemented(fname,"Unrecognized interaction_type.");
    }

    psq *= mass;

    //---- Evaluate the interaction as a Taylor series in p^2: O(p) = sum_n C_n p^(2n)
    interaction[i] = 0.0;
    for (int j=0; j<interaction_arg.num_couplings; j++) {
      interaction[i] += interaction_arg.couplings[j]*pow(psq,j);
    }

    //---- Make sure O(p) is non-negative before taking the square root!
    if (interaction[i]<0.0) {
      ERR.General(fname, "Interaction is less the zero; cannot take square root.");
    } else{
      interaction[i] *= 2.0 * TWOPI / mass;  // Just a convention
      interaction[i] = sqrt(interaction[i]); // Division by vol because FFT is not normalized
      interaction[i] /= vol; // Division by vol because FFT is not normalized
      interaction[i] *= dt; // Controls sign and step size of interaction
    }

  }

}

Interaction::~Interaction()
{

  char* fname = "void Interaction::~Interaction()";
  VRB.Func(fname);

  free(interaction);

}

void Interaction::Evolve(double* b)
{
  // Apply opererator onto in and return as out (i.e., b -> interaction.b ).
  // In and out fermion vectors are in the fermion storage order: [x][y][z][comp]

  char* fname = "void Interaction::Evolve(double*)";
  VRB.Func(fname);

  double* field = lattice->Get();
  int vol = GJP.Vol();
  int j;

  //---- Copy in vector so that it is not detroyed
  for (int i=0; i<vol; i++) {
    j = 2*i;
    Fourier::in[i][0] = b[j]   * interaction[i];
    Fourier::in[i][1] = b[j+1] * interaction[i];
  }

  //---- Go to position space
  Fourier::Forward();

  //---- Add APBC phases...
  if ( GJP.APBCQuery() ) {
    Fourier::phase.Add(Fourier::out);
  }

  // Apply field
  for (int i=0; i<vol; i++) {
    j = 2*i;
    Fourier::in[i][0] = field[j] * Fourier::out[i][0] - field[j+1] * Fourier::out[i][1]; 
    Fourier::in[i][1] = field[j] * Fourier::out[i][1] + field[j+1] * Fourier::out[i][0]; 
  }

  //---- Remove APBC phases...
  if ( GJP.APBCQuery() ) {
    Fourier::phase.Subtract(Fourier::in);
  }

  //---- Go to momentum space
  Fourier::Backward();

  //---- Add interaction to in vector
  for (int i=0; i<vol; i++) {
    j = 2*i;
    b[j]   += Fourier::out[i][0]; 
    b[j+1] += Fourier::out[i][1]; 
  }

}

Interaction2::Interaction2(Lattice* lat, InteractionArg interaction_arg, KineticArg kinetic_arg1, KineticArg kinetic_arg2, double dt)
{

  char* fname = "Interaction::Interaction(Lattice*, InteractionArg, PropagatorArg, PropagatorArg, double)";
  VRB.Func(fname);

  lattice = lat;

  //---- Compute the interactions in momentum space
  int vol = GJP.Vol();
  interaction = (double *) malloc(vol*sizeof(double));  

  Dispersion dispersion1(kinetic_arg1.dispersion_type, kinetic_arg1.mass, CUTOFF_TYPE_NONE);
  Dispersion dispersion2(kinetic_arg2.dispersion_type, kinetic_arg2.mass, CUTOFF_TYPE_NONE);
  double xi1;
  double xi2;

  double psq;
  double mass = 2.0 * kinetic_arg1.mass * kinetic_arg2.mass; // reduced mass
  mass /= kinetic_arg1.mass + kinetic_arg2.mass; // reduced mass

  double PSQ = 2.0*( PI - GJP.Cutoff() );
  PSQ *= PSQ;
  double XI1 = exp( PSQ/(2.0*kinetic_arg1.mass));
  double XI2 = exp( PSQ/(2.0*kinetic_arg2.mass));

  for(int i=0; i<vol; i++) {

    xi1 = 1.0 + dispersion1.Get(i);
    if (xi1 > XI1) { xi1 = XI1; }
    xi2 = 1.0 + dispersion2.Get(i);
    if (xi2 > XI2) { xi2 = XI2; }

    switch (interaction_arg.interaction_type) {
      case INTERACTION_TYPE_NONE:
        ERR.NotImplemented(fname,"Interaction type INTERACTION_TYPE_NONE.");
        break;
      case INTERACTION_TYPE_ONEMINUSXIINVSQ:
        psq = 1.0 - 1.0 /(xi1*xi2); 
        break;
      case INTERACTION_TYPE_XISQMINUSONE:
        psq = xi1*xi2 - 1.0; 
        break;
      default:
        ERR.NotImplemented(fname,"Unrecognized interaction_type.");
    }

    psq *= mass;

    //---- Evaluate the interaction as a Taylor series in p^2: O(p) = sum_n C_n p^(2n)
    interaction[i] = 0.0;
    for (int j=0; j<interaction_arg.num_couplings; j++) {
      interaction[i] += interaction_arg.couplings[j]*pow(psq,j);
    }

    //---- Make sure O(p) is non-negative before taking the square root!
    if (interaction[i]<0.0) {
      ERR.General(fname, "Interaction is less the zero; cannot take square root.");
    }

    //---- Normalize operator, etc..
    interaction[i] *= 2.0 * TWOPI / mass;  // Just a convention
    interaction[i] = sqrt(interaction[i]); // Division by vol because FFT is not normalized
    interaction[i] /= vol; // Division by vol because FFT is not normalized
    interaction[i] *= dt; // Controls overall sign and step size of the interaction

  }

  tmp = (double *) malloc(2*vol*sizeof(double));  

}

Interaction2::~Interaction2()
{

  char* fname = "void Interaction::~Interaction()";
  VRB.Func(fname);

  free(tmp);
  free(interaction);

}

void Interaction2::Evolve(double* b)
{
  // Apply opererator onto in and return as out (i.e., b -> interaction.b ).
  // In and out fermion vectors are in the fermion storage order: [x][y][z][comp]

  char* fname = "void Interaction::Evolve(double*)";
  VRB.Func(fname);

  double* field = lattice->Get();
  int vol = GJP.Vol();
  int j;

  //---- First apply interaction onto lattice
  for (int i=0; i<vol; i++) {
    j = 2*i;
    Fourier::in[i][0] = field[j];
    Fourier::in[i][1] = field[j+1];
  }

  if ( GJP.APBCQuery() ) {
    Fourier::phase.Subtract(Fourier::in);
  }

  Fourier::Backward();

  for (int i=0; i<vol; i++) {
    Fourier::in[i][0] = interaction[i] * Fourier::out[i][0] /vol; 
    Fourier::in[i][1] = interaction[i] * Fourier::out[i][1]/ vol; 
  }

  Fourier::Forward();

  if ( GJP.APBCQuery() ) {
    Fourier::phase.Add(Fourier::out);
  }

  for (int i=0; i<vol; i++) {
    j = 2*i;
    tmp[j] = Fourier::out[i][0]; 
    tmp[j+1] = Fourier::out[i][1]; 
  }

  //---- Copy in vector so that it is not detroyed
  for (int i=0; i<vol; i++) {
    j = 2*i;
    Fourier::in[i][0] = b[j];
    Fourier::in[i][1] = b[j+1];
  }

  //---- Go to position space
  Fourier::Forward();

  //---- Add APBC phases...
  if ( GJP.APBCQuery() ) {
    Fourier::phase.Add(Fourier::out);
  }

  // Apply field
  for (int i=0; i<vol; i++) {
    j = 2*i;
    Fourier::in[i][0] = tmp[j] * Fourier::out[i][0] - tmp[j+1] * Fourier::out[i][1]; 
    Fourier::in[i][1] = tmp[j] * Fourier::out[i][1] + tmp[j+1] * Fourier::out[i][0]; 
  }

  //---- Remove APBC phases...
  if ( GJP.APBCQuery() ) {
    Fourier::phase.Subtract(Fourier::in);
  }

  //---- Go to momentum space
  Fourier::Backward();

  //---- Add interaction to in vector
  for (int i=0; i<vol; i++) {
    j = 2*i;
    b[j]   += Fourier::out[i][0]; 
    b[j+1] += Fourier::out[i][1]; 
  }

}

Interaction3::Interaction3(Lattice* lat, InteractionArg interaction_arg, KineticArg kinetic_arg1, KineticArg kinetic_arg2, double dt)
{

  char* fname = "Interaction::Interaction(Lattice*, InteractionArg, PropagatorArg, PropagatorArg, double)";
  VRB.Func(fname);

  lattice = lat;

  //---- Compute the interactions in momentum space
  int vol = GJP.Vol();
  interaction = (double *) malloc(vol*sizeof(double));  

  Dispersion dispersion1(kinetic_arg1.dispersion_type, kinetic_arg1.mass, CUTOFF_TYPE_NONE);
  Dispersion dispersion2(kinetic_arg2.dispersion_type, kinetic_arg2.mass, CUTOFF_TYPE_NONE);
  double xi1;
  double xi2;

  double psq;
  double mass = 2.0 * kinetic_arg1.mass * kinetic_arg2.mass; // reduced mass
  mass /= kinetic_arg1.mass + kinetic_arg2.mass; // reduced mass

  double PSQ = GJP.Cutoff();
  PSQ *= PSQ;
  double XI1 = exp( PSQ/(2.0*kinetic_arg1.mass));
  double XI2 = exp( PSQ/(2.0*kinetic_arg2.mass));

  for(int i=0; i<vol; i++) {

    xi1 = 1.0 + dispersion1.Get(i);
    if (xi1 > XI1) { xi1 = XI1; }
    xi2 = 1.0 + dispersion2.Get(i);
    if (xi2 > XI2) { xi2 = XI2; }

    switch (interaction_arg.interaction_type) {
      case INTERACTION_TYPE_NONE:
        ERR.NotImplemented(fname,"Interaction type INTERACTION_TYPE_NONE.");
        break;
      case INTERACTION_TYPE_ONEMINUSXIINVSQ:
        psq = 1.0 - 1.0 /(xi1*xi2); 
        break;
      case INTERACTION_TYPE_XISQMINUSONE:
        psq = xi1*xi2 - 1.0; 
        break;
      default:
        ERR.NotImplemented(fname,"Unrecognized interaction_type.");
    }

    psq *= mass;

    //---- Evaluate the interaction as a Taylor series in p^2: O(p) = sum_n C_n p^(2n)
    interaction[i] = 0.0;
    for (int j=0; j<interaction_arg.num_couplings; j++) {
      interaction[i] += interaction_arg.couplings[j]*pow(psq,j);
    }

    //---- Make sure O(p) is non-negative before taking the square root!
    if (interaction[i]<0.0) {
      ERR.General(fname, "Interaction is less the zero; cannot take square root.");
    }

    //---- Normalize operator, etc..
    interaction[i] *= 2.0 * TWOPI / mass;  // Just a convention
    interaction[i] = sqrt(interaction[i]); // Division by vol because FFT is not normalized
    interaction[i] /= vol; // Division by vol because FFT is not normalized
//    cout << i << " " << interaction[i] << "\n";
    interaction[i] *= dt; // Controls overall sign and step size of the interaction

  }

  tmp = (double *) malloc(2*vol*sizeof(double));  

}

Interaction3::~Interaction3()
{

  char* fname = "void Interaction::~Interaction()";
  VRB.Func(fname);

  free(tmp);
  free(interaction);

}

void Interaction3::Evolve(double* b)
{
  // Apply opererator onto in and return as out (i.e., b -> interaction.b ).
  // In and out fermion vectors are in the fermion storage order: [x][y][z][comp]

  char* fname = "void Interaction::Evolve(double*)";
  VRB.Func(fname);

  double* field = lattice->Get();
  int vol = GJP.Vol();
  int j;

  //---- First apply interaction onto lattice
  for (int i=0; i<vol; i++) {
    j = 2*i;
    Fourier::in[i][0] = field[j];
    Fourier::in[i][1] = field[j+1];
  }

  if ( GJP.APBCQuery() ) {
    Fourier::phase.Subtract(Fourier::in);
  }

  Fourier::Backward();

  for (int i=0; i<vol; i++) {
    Fourier::in[i][0] = interaction[i] * Fourier::out[i][0] /vol; 
    Fourier::in[i][1] = interaction[i] * Fourier::out[i][1]/ vol; 
  }

  Fourier::Forward();

  if ( GJP.APBCQuery() ) {
    Fourier::phase.Add(Fourier::out);
  }

  for (int i=0; i<vol; i++) {
    j = 2*i;
    tmp[j] = Fourier::out[i][0]; 
    tmp[j+1] = Fourier::out[i][1]; 
  }

  //---- Copy in vector so that it is not detroyed
  for (int i=0; i<vol; i++) {
    j = 2*i;
    Fourier::in[i][0] = b[j];
    Fourier::in[i][1] = b[j+1];
  }

  //---- Go to position space
  Fourier::Forward();

  //---- Add APBC phases...
  if ( GJP.APBCQuery() ) {
    Fourier::phase.Add(Fourier::out);
  }

  // Apply field
  for (int i=0; i<vol; i++) {
    j = 2*i;
    Fourier::in[i][0] = tmp[j] * Fourier::out[i][0] - tmp[j+1] * Fourier::out[i][1]; 
    Fourier::in[i][1] = tmp[j] * Fourier::out[i][1] + tmp[j+1] * Fourier::out[i][0]; 
  }

  //---- Remove APBC phases...
  if ( GJP.APBCQuery() ) {
    Fourier::phase.Subtract(Fourier::in);
  }

  //---- Go to momentum space
  Fourier::Backward();

  //---- Add interaction to in vector
  for (int i=0; i<vol; i++) {
    j = 2*i;
    b[j]   += Fourier::out[i][0]; 
    b[j+1] += Fourier::out[i][1]; 
  }

}


Potential::Potential(PotentialArg potential_arg, double dt)
{

  char* fname = "Potential::Potential(PotentialArg, double)";
  VRB.Func(fname);

  if ( GJP.APBCQuery() ) { ERR.NotImplemented(fname,"Must use PBCs when a potential is present."); }

  int vol = GJP.Vol();
  potential = (double *) malloc(vol*sizeof(double));  

  //---- Compute the potential
  int x;
  int y;
  int z;
  int Y; // Y = y+x*Y

  int x_sites = GJP.Xsites();
  int y_sites = GJP.Ysites();
  int z_sites = GJP.Zsites();

  double v; 

  for(int i=0; i<vol; i++) {

    z = i%z_sites;
    Y = i/z_sites;
    y = Y%y_sites;
    x = Y/y_sites;

    if (x>=x_sites/2) { x -= x_sites; }
    if (y>=y_sites/2) { y -= y_sites; }
    if (z>=z_sites/2) { z -= z_sites; }

    //----- Potential form
    switch (potential_arg.potential_form) {
      case POTENTIAL_FORM_NONE:
        v = 0.0;
        break;
      case POTENTIAL_FORM_HARMONIC:
        v =  potential_arg.spring_constant1*x*x;
        v += potential_arg.spring_constant2*y*y;
        v += potential_arg.spring_constant3*z*z;
        v /= 2.0;
        break;
      case POTENTIAL_FORM_COULOMB:
        ERR.NotImplemented(fname, "Potential form POTENTIAL_FORM_COULOMB not supported.");
        v =  x*x/potential_arg.spring_constant1;
        v += y*y/potential_arg.spring_constant2;
        v += z*z/potential_arg.spring_constant3;
        v = 1.0/sqrt(v);
        break;
      default:
        ERR.NotImplemented(fname, "Potential form not supported.");
    }

    //---- Potential discretization
    switch (potential_arg.potential_type) {
      case POTENTIAL_TYPE_LIN:
        potential[i] = dt*v;
        potential[i] /= vol;
        break;
      case POTENTIAL_TYPE_EXP:
        potential[i] = 1.0 - exp(-dt*v);
        potential[i] /= vol;
        break;
      default:
        ERR.NotImplemented(fname, "Potential type not supported.");
    }

  }

}

Potential::~Potential()
{

  char* fname = "void SplitPotential::~Potential()";
  VRB.Func(fname);

  free(potential);

}

void Potential::Evolve(double* b)
{
  // Apply opererator onto in and return as out (i.e., b -> interaction.a ).
  // In and out fermion vectors are in the fermion storage order: [x][y][z][comp]

  char* fname = "void Potential::Evolve(double*)";
  VRB.Func(fname);

  int vol = GJP.Vol();
  int j;

  //---- Copy in vector so that it is not detroyed
  for (int i=0; i<vol; i++) {
    j = 2*i;
    Fourier::in[i][0] = b[j];
    Fourier::in[i][1] = b[j+1];
  }

  //---- Go to position space
  Fourier::Forward();

  for(int i=0; i<vol; i++) {
    Fourier::in[i][0] = Fourier::out[i][0] * potential[i];
    Fourier::in[i][1] = Fourier::out[i][1] * potential[i];
  }

  //---- Go to momentum space
  Fourier::Backward();

  //---- Add interaction to in vector
  for (int i=0; i<vol; i++) {
    j = 2*i;
    b[j]   -= Fourier::out[i][0]; 
    b[j+1] -= Fourier::out[i][1]; 
  }

}




Kinetic::Kinetic(KineticArg kinetic_arg, double dt)
{

  char* fname = "Kinetic::Kinetic(KineticArg, double)";
  VRB.Func(fname);

  int vol = GJP.Vol();
  xi = (double *) malloc(2*vol*sizeof(double));

  double d;
  Dispersion dispersion(kinetic_arg.dispersion_type, kinetic_arg.mass/dt, CUTOFF_TYPE_HARD);
  for (int i=0; i<vol; i++) {
    d = 1.0 + dispersion.Get(i);
    xi[2*i]   = d;
    xi[2*i+1] = d;
  }


}

Kinetic::~Kinetic()
{

  char* fname = "void Kinetic::~Kinetic()";
  VRB.Func(fname);

  free(xi);

}

void Kinetic::Evolve(double* b)
{
  // Apply opererator onto in and return as out (i.e., b -> interaction.a ).
  // In and out fermion vectors are in the fermion storage order: [x][y][z][comp]

  char* fname = "void Kinetic::Evolve(double*)";
  VRB.Func(fname);

  //---- Apply D^-1 to source b; b ->  D^-1 b
  for (int i=0; i<2*GJP.Vol(); i++) {
    b[i] /= xi[i];
  }

}



