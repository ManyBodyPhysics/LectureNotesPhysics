#include <string>
#include <iostream>
#include <stdarg.h>
using namespace std;
#include "arg.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"
#include "constants.h"

GlobalJobParameter GJP;

GlobalJobParameter::GlobalJobParameter()
{
  t_sites = 4;
  x_sites = 4;
  y_sites = 4;
  z_sites = 4;
  x_boundary_type = BOUNDARY_TYPE_PRD;
  y_boundary_type = BOUNDARY_TYPE_PRD;
  z_boundary_type = BOUNDARY_TYPE_PRD;
  cutoff = 1.0*PI; 
}

GlobalJobParameter::~GlobalJobParameter()
{
}


void GlobalJobParameter::Initialize(DoArg do_arg)
{

  char* fname = "void GlobalJobParameter::initialize(DoArg)";
  VRB.Func(fname);

  t_sites = do_arg.t_sites;
  x_sites = do_arg.x_sites;
  y_sites = do_arg.y_sites;
  z_sites = do_arg.z_sites;

  if ( x_sites%2!=0 ) { 
    ERR.General(fname, "x_sites not divisible by two."); 
  }
  if ( y_sites%2!=0 ) { 
    ERR.General(fname, "y_sites not divisible by two."); 
  }
  if ( z_sites%2!=0 ) { 
    ERR.General(fname, "z_sites not divisible by two."); 
  }

  x_boundary_type = do_arg.x_boundary_type;
  y_boundary_type = do_arg.y_boundary_type;
  z_boundary_type = do_arg.z_boundary_type;

//  if (do_arg.cutoff > 1.0) {  ERR.General(fname, "Cutoff factor is greater than unity.");}
  if (do_arg.cutoff <= 0.0) {  ERR.General(fname, "Cutoff factor is not positive.");}
  cutoff = do_arg.cutoff*PI;
 
}

int GlobalJobParameter::Tsites() { return t_sites; }
int GlobalJobParameter::Xsites() { return x_sites; }
int GlobalJobParameter::Ysites() { return y_sites; }
int GlobalJobParameter::Zsites() { return z_sites; }
int GlobalJobParameter::Vol() { return x_sites*y_sites*z_sites; }

BoundaryType GlobalJobParameter::Xbc() { return x_boundary_type; }
BoundaryType GlobalJobParameter::Ybc() { return y_boundary_type; }
BoundaryType GlobalJobParameter::Zbc() { return z_boundary_type; }

bool GlobalJobParameter::APBCQuery()
{
  if ( (x_boundary_type==BOUNDARY_TYPE_APRD)||(y_boundary_type==BOUNDARY_TYPE_APRD)||(z_boundary_type==BOUNDARY_TYPE_APRD) ) {
    return true;
  } else { 
    return false;
  }
}

void GlobalJobParameter::GetBcShift(double* bc)
{
  if (x_boundary_type == BOUNDARY_TYPE_PRD) {
    bc[0] = 0.0;
  } else if (x_boundary_type == BOUNDARY_TYPE_APRD) {
    bc[0] = 0.5;
  }
  if (y_boundary_type == BOUNDARY_TYPE_PRD) {
    bc[1] = 0.0;
  } else if (y_boundary_type == BOUNDARY_TYPE_APRD) {
    bc[1] = 0.5;
  }
  if (z_boundary_type == BOUNDARY_TYPE_PRD) {
    bc[2] = 0.0;
  } else if (z_boundary_type == BOUNDARY_TYPE_APRD) {
    bc[2] = 0.5;
  }
}

double GlobalJobParameter::Cutoff() { return cutoff; }
double GlobalJobParameter::CutoffSq() { return cutoff*cutoff; }
