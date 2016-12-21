#include <fstream>
#include <string>
#include <iostream>
using namespace std;
#include "arg.h"
#include "verbose.h"
#include "error.h"

DoArg::DoArg()
{
  t_sites = 4;
  x_sites = 4;
  y_sites = 4;
  z_sites = 4;
  x_boundary_type = BOUNDARY_TYPE_PRD;
  y_boundary_type = BOUNDARY_TYPE_PRD;
  z_boundary_type = BOUNDARY_TYPE_PRD;
  cutoff = 1.0;
}

void DoArg::Decode(string vml_filename)
{
  char* fname = "void DoArg::Decode(string)";
  VRB.Func(fname);

  string input;
  ifstream file;
  file.open(vml_filename.c_str());
  if (!file) {
    ERR.FileR(fname,vml_filename.c_str());
  }
  file >> t_sites;
  file >> x_sites;
  file >> y_sites;
  file >> z_sites;

  file >> input;
  if (input=="BOUNDARY_TYPE_PRD") {
    x_boundary_type = BOUNDARY_TYPE_PRD;
  } else if (input=="BOUNDARY_TYPE_APRD") {
    x_boundary_type = BOUNDARY_TYPE_APRD;
  } else {
    ERR.General(fname, "Unable to recognize x_boundary_type.");
  }

  file >> input;
  if (input=="BOUNDARY_TYPE_PRD") {
    y_boundary_type = BOUNDARY_TYPE_PRD;
  } else if (input=="BOUNDARY_TYPE_APRD") {
    y_boundary_type = BOUNDARY_TYPE_APRD;
  } else {
    ERR.General(fname, "Unable to recognize y_boundary_type.");
  }

  file >> input;
  if (input=="BOUNDARY_TYPE_PRD") {
    z_boundary_type = BOUNDARY_TYPE_PRD;
  } else if (input=="BOUNDARY_TYPE_APRD") {
    z_boundary_type = BOUNDARY_TYPE_APRD;
  } else {
    ERR.General(fname, "Unable to recognize z_boundary_type.");
  }

  file >> cutoff;

  file.close();
}

LatticeArg::LatticeArg()
{
  field_type = FIELD_TYPE_GAUSS;
}

void LatticeArg::Decode(string vml_filename)
{

  char* fname = "void LatticeArg::Decode(string)";
  VRB.Func(fname);
  string input;
  ifstream file;
  file.open(vml_filename.c_str());
  if (!file) {
    ERR.FileR(fname,vml_filename.c_str());
  }

  file >> input;
  if (input=="FIELD_TYPE_GAUSS") {
    field_type = FIELD_TYPE_GAUSS;
  } else if (input=="FIELD_TYPE_COMPLEXGAUSS") {
    field_type = FIELD_TYPE_COMPLEXGAUSS;
  } else if (input=="FIELD_TYPE_ZTWO") {
    field_type = FIELD_TYPE_ZTWO;
  } else if (input=="FIELD_TYPE_ZTHREE") {
    field_type = FIELD_TYPE_ZTHREE;
  } else {
    ERR.General(fname, "Unable to recognize field_type.");
  }

  file.close();

}

RandomArg::RandomArg()
{
  random_type = RANDOM_TYPE_RAN2;
  seed_type = SEED_TYPE_INPUT;
  seed = 7012;
  file_stem = "rng";
}

void RandomArg::Decode(string vml_filename)
{
  char* fname = "void RandomArg::Decode(string)";
  VRB.Func(fname);

  string input;
  ifstream file;
  file.open(vml_filename.c_str());
  if (!file) {
    ERR.FileR(fname,vml_filename.c_str());
  }

  file >> input;
  if (input=="RANDOM_TYPE_RAN0") {
    random_type = RANDOM_TYPE_RAN0;
  } else if (input=="RANDOM_TYPE_RAN1") {
    random_type = RANDOM_TYPE_RAN1;
  } else if (input=="RANDOM_TYPE_RAN2") {
    random_type = RANDOM_TYPE_RAN2;
  } else if (input=="RANDOM_TYPE_RAN3") {
    random_type = RANDOM_TYPE_RAN3;
  } else if (input=="RANDOM_TYPE_RANLXS0") {
    random_type = RANDOM_TYPE_RANLXS0;
  } else if (input=="RANDOM_TYPE_RANLXS1") {
    random_type = RANDOM_TYPE_RANLXS1;
  } else if (input=="RANDOM_TYPE_RANLXS2") {
    random_type = RANDOM_TYPE_RANLXS2;
  } else if (input=="RANDOM_TYPE_RANLXD1") {
    random_type = RANDOM_TYPE_RANLXD1;
  } else if (input=="RANDOM_TYPE_RANLXD2") {
    random_type = RANDOM_TYPE_RANLXD2;
  } else {
    ERR.General(fname, "Unable to recognize random_type.");
  }


  file >> input;
  if (input=="SEED_TYPE_DEFAULT") {
    seed_type = SEED_TYPE_DEFAULT;
  } else if (input=="SEED_TYPE_INPUT") {
    seed_type = SEED_TYPE_INPUT;
  } else if (input=="SEED_TYPE_TIME") {
    seed_type = SEED_TYPE_TIME;
  } else if (input=="SEED_TYPE_FILE") {
    seed_type = SEED_TYPE_FILE;
  } else {
    ERR.General(fname, "Unable to recognize seed_type.");
  }
  file >> seed;
  file >> file_stem;

  file.close();
}

PropagatorArg::PropagatorArg()
{
  n_species = 1;
  file_stem = "propagator";

}

void PropagatorArg::Decode(string vml_filename)
{

  char* fname = "void PropagatorArg::Decode(string)";
  VRB.Func(fname);

//  string input;
  ifstream file;
  file.open(vml_filename.c_str());
  if (!file) {
    ERR.FileR(fname,vml_filename.c_str());
  }

  file >> n_species;
  file >> file_stem;

  file.close();

}

OneBodyArg::OneBodyArg()
{
  source_type = SOURCE_TYPE_MOM;
  lambda1 = 1.0;
  lambda2 = 1.0;
  lambda3 = 1.0;
}

void OneBodyArg::Decode(string vml_filename)
{

  char* fname = "void OneBodyArg::Decode(string)";
  VRB.Func(fname);

  string input;
  ifstream file;
  file.open(vml_filename.c_str());
  if (!file) {
    ERR.FileR(fname,vml_filename.c_str());
  }

  file >> input;
  if (input=="SOURCE_TYPE_MOM") {
    source_type = SOURCE_TYPE_MOM;
  } else if (input=="SOURCE_TYPE_SHO") {
    source_type = SOURCE_TYPE_SHO;
  } else {
    ERR.General(fname, "Unable to recognize source_type.");
  }

  file >> lambda1;
  file >> lambda2;
  file >> lambda3;

  file.close();
}

TwoBodyArg::TwoBodyArg()
{
  wavefunc_type = WAVEFUNC_TYPE_UNIFORM;
  dispersion_type1 = DISPERSION_TYPE_STANDARD;
  dispersion_type2 = DISPERSION_TYPE_STANDARD;
  mass1 = 1.0;
  mass2 = 1.0;
  lambda = 1.0; 
  file_stem = "two_body";
}

void TwoBodyArg::Decode(string vml_filename)
{

  char* fname = "void TwoArg::Decode(string)";
  VRB.Func(fname);

  string input;
  ifstream file;
  file.open(vml_filename.c_str());
  if (!file) {
    ERR.FileR(fname,vml_filename.c_str());
  }

  file >> input;
  if (input=="WAVEFUNC_TYPE_NONE") {
    wavefunc_type = WAVEFUNC_TYPE_NONE;
  } else if (input=="WAVEFUNC_TYPE_UNIFORM") {
    wavefunc_type = WAVEFUNC_TYPE_UNIFORM;
  } else if (input=="WAVEFUNC_TYPE_GND") {
    wavefunc_type = WAVEFUNC_TYPE_GND;
  } else if (input=="WAVEFUNC_TYPE_TGND") {
    wavefunc_type = WAVEFUNC_TYPE_TGND;
  } else if (input=="WAVEFUNC_TYPE_PAIR1") {
    wavefunc_type = WAVEFUNC_TYPE_PAIR1;
  } else if (input=="WAVEFUNC_TYPE_PAIR2") {
    wavefunc_type = WAVEFUNC_TYPE_PAIR2;
  } else {
    ERR.General(fname, "Unable to recognize wavefunc_type.");
  }

  file >> input;
  if (input=="DISPERSION_TYPE_STANDARD") {
    dispersion_type1 = DISPERSION_TYPE_STANDARD;
  } else if (input=="DISPERSION_TYPE_QUADRATIC") {
    dispersion_type1 = DISPERSION_TYPE_QUADRATIC;
  } else if (input=="DISPERSION_TYPE_PERFECT") {
    dispersion_type1 = DISPERSION_TYPE_PERFECT;
  } else {
    ERR.General(fname, "Unable to recognize dispersion_type.");
  }

  file >> input;
  if (input=="DISPERSION_TYPE_STANDARD") {
    dispersion_type2 = DISPERSION_TYPE_STANDARD;
  } else if (input=="DISPERSION_TYPE_QUADRATIC") {
    dispersion_type2 = DISPERSION_TYPE_QUADRATIC;
  } else if (input=="DISPERSION_TYPE_PERFECT") {
    dispersion_type2 = DISPERSION_TYPE_PERFECT;
  } else {
    ERR.General(fname, "Unable to recognize dispersion_type.");
  }

  file >> mass1;
  file >> mass2;
  file >> lambda;
  file >> file_stem;

  file.close();
}

VerboseArg::VerboseArg()
{
  func_level = true;
  warn_level = true;
  result_level = true;
  flow_level = true;
  debug_level = true;
}

void VerboseArg::Decode(string vml_filename)
{

  char* fname = "void VerboseArg::Decode(string)";
  VRB.Func(fname);

  string input;
  ifstream file;
  file.open(vml_filename.c_str());
  if (!file) {
    ERR.FileR(fname,vml_filename.c_str());
  }

  file >> input;
  if (input == "1") { func_level=true; } else { func_level=false; }
  file >> input;
  if (input == "1") { warn_level=true; } else { warn_level=false; }
  file >> input;
  if (input == "1") { result_level=true; } else { result_level=false; }
  file >> input;
  if (input == "1") { flow_level=true; } else { flow_level=false; }
  file >> input;
  if (input == "1") { debug_level=true; } else { debug_level=false; }

  file.close();

}

EvoArg::EvoArg()
{
  start = 0;
  unload_period = 1;
  configurations = 10;
}

void EvoArg::Decode(string vml_filename)
{

  char* fname = "void EvoArg::Decode(string)";
  VRB.Func(fname);

  ifstream file;
  file.open(vml_filename.c_str());
  if (!file) {
    ERR.FileR(fname,vml_filename.c_str());
  }

  file >> start;
  file >> unload_period;
  file >> configurations;

  file.close();
}


MomentaArg::MomentaArg()
{
   dispersion_type = DISPERSION_TYPE_STANDARD;
   fermi_energy = 1.0;
   mass = 1.0;
}

void MomentaArg::Decode(string vml_filename)
{

  char* fname = "void MomentaArg::Decode(string)";
  VRB.Func(fname);

  string input;
  ifstream file;
  file.open(vml_filename.c_str());
  if (!file) {
    ERR.FileR(fname,vml_filename.c_str());
  }

  file >> input;
  if (input=="DISPERSION_TYPE_STANDARD") {
    dispersion_type = DISPERSION_TYPE_STANDARD;
  } else if (input=="DISPERSION_TYPE_QUADRATIC") {
    dispersion_type = DISPERSION_TYPE_QUADRATIC;
  } else if (input=="DISPERSION_TYPE_PERFECT") {
    dispersion_type = DISPERSION_TYPE_PERFECT;
  } else {
    ERR.General(fname, "Unable to recognize dispersion_type.");
  }

  file >> fermi_energy;
  file >> mass;

  file.close();
}

InteractionArg::InteractionArg()
{
  interaction_type = INTERACTION_TYPE_NONE;
  num_couplings = 0;
}

void InteractionArg::Decode(string vml_filename)
{

  char* fname = "void InteractionArg::Decode(string)";
  VRB.Func(fname);

  string input;
  ifstream file;
  file.open(vml_filename.c_str());
  if (!file) {
    ERR.FileR(fname,vml_filename.c_str());
  }

  file >> input;
  if (input=="INTERACTION_TYPE_NONE") {
    interaction_type = INTERACTION_TYPE_NONE;
  } else if  (input=="INTERACTION_TYPE_ONEMINUSXIINVSQ") {
    interaction_type = INTERACTION_TYPE_ONEMINUSXIINVSQ;
  } else if (input=="INTERACTION_TYPE_XISQMINUSONE") {
    interaction_type = INTERACTION_TYPE_XISQMINUSONE;
  } else {
    ERR.General(fname, "Unable to recognize interaction_type.");
  }

  double coupling;
  while (true) {
    file >> coupling;
    if ( file.eof() ) break;
    couplings.push_back(coupling);
  }

  file.close();

  num_couplings = couplings.size();

  for (int n=0; n<num_couplings; n++) {
    VRB.Flow(fname,"Coupling %d of %d = %e", n+1, num_couplings, couplings[n]);
  }

}

PotentialArg::PotentialArg()
{
  potential_form = POTENTIAL_FORM_NONE;
  potential_type = POTENTIAL_TYPE_LIN;
  spring_constant1 = 0.0;
  spring_constant2 = 0.0;
  spring_constant3 = 0.0;
}

void PotentialArg::Decode(string vml_filename)
{

  char* fname = "void PotentialArg::Decode(string)";
  VRB.Func(fname);

  string input;
  ifstream file;
  file.open(vml_filename.c_str());
  if (!file) {
    ERR.FileR(fname,vml_filename.c_str());
  }

  file >> input;
  if (input=="POTENTIAL_FORM_NONE") {
    potential_form = POTENTIAL_FORM_NONE;
  } else if (input=="POTENTIAL_FORM_HARMONIC") {
    potential_form = POTENTIAL_FORM_HARMONIC;
  } else if (input=="POTENTIAL_FORM_COULOMB") {
    potential_form = POTENTIAL_FORM_COULOMB;
  } else {
    ERR.General(fname, "Unable to recognize potential_type.");
  }

  file >> input;
  if (input=="POTENTIAL_TYPE_LIN") {
    potential_type = POTENTIAL_TYPE_LIN;
  } else if (input=="POTENTIAL_TYPE_EXP") {
    potential_type = POTENTIAL_TYPE_EXP;
  } else {
    ERR.General(fname, "Unable to recognize potential_type.");
  }

  file >> spring_constant1;
  file >> spring_constant2;
  file >> spring_constant3;

  file.close();
}

KineticArg::KineticArg()
{
  dispersion_type = DISPERSION_TYPE_STANDARD;
  mass = 1.0;
}

void KineticArg::Decode(string vml_filename)
{

  char* fname = "void KineticArg::Decode(string)";
  VRB.Func(fname);

  string input;
  ifstream file;
  file.open(vml_filename.c_str());
  if (!file) {
    ERR.FileR(fname,vml_filename.c_str());
  }

  file >> input;
  if (input=="DISPERSION_TYPE_STANDARD") {
    dispersion_type = DISPERSION_TYPE_STANDARD;
  } else if (input=="DISPERSION_TYPE_QUADRATIC") {
    dispersion_type = DISPERSION_TYPE_QUADRATIC;
  } else if (input=="DISPERSION_TYPE_PERFECT") {
    dispersion_type = DISPERSION_TYPE_PERFECT;
  } else {
    ERR.General(fname, "Unable to recognize dispersion_type.");
  }

  file >> mass;

  file.close();

}



