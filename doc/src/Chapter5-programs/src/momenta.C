#include <iostream>
#include "math.h"
using namespace std;
#include "momenta.h"
#include "dispersion.h"
#include "global_job_parameter.h"
#include "error.h"
#include "verbose.h"
#include "arg.h"
#include "constants.h"

Momenta::Momenta(MomentaArg momenta_arg)
{

  char* fname = "Momenta::Momenta(MomentaArg)";

  int x_sites = GJP.Xsites();
  int y_sites = GJP.Ysites();
  int z_sites = GJP.Zsites();

  double energy;
  Dispersion dispersion(momenta_arg.dispersion_type, momenta_arg.mass, CUTOFF_TYPE_HARD);

  //---- Perform two passes. First pass determines how many momenta; second pass, allocates enough memory and stores the momenta.
  int count;
  for (int i=0; i<2;i++) {
    if (i==1) {
      if (count==0) {
        //--- Error and quit
        ERR.General(fname,"No energy levels exceed the specified fermi level.");
      }
      //---- Allocate memory to store results
      momenta_count = count;
      omega = new double[count];
      momentum = new vector<int>[count];
      parity_list = new int[count];
    }
    count = 0;
    for (int n1=0; n1<x_sites; n1++) {
      for (int n2=0; n2<y_sites; n2++) {
        for (int n3=0; n3<z_sites; n3++) {
          energy = log(1+dispersion.Get(n1, n2, n3));
          if (energy <= momenta_arg.fermi_energy) {
            if (i==1) {
              //---- Store momentum config and corresponding energy in memory
              omega[count] = energy;
              momentum[count].push_back(n1);
              momentum[count].push_back(n2);
              momentum[count].push_back(n3);
            }
            count += 1;
          }
        }
      }
    }
  }

  //---- Sort momenta according to energy using insertion sort
  //---- I know, I know... not the fastest of sorting methods.
  vector<int> tmp;
  int j;
  for (int i=1; i<momenta_count; i++) {
    energy = omega[i];
    tmp = momentum[i];
    j = i-1;
    while ( (j>0)&&(omega[j]>energy) ) {
      omega[j+1] = omega[j];
      momentum[j+1] = momentum[j];
      j--;
    }
    omega[j+1]=energy;
    momentum[j+1]=tmp;
  }

  //---- Compute the shell count.
  shell_count.push_back(1);
  for (int i=1; i<momenta_count; i++) {
    if (omega[i]==omega[i-1]) {
        shell_count[shell_count.size()-1]++;
    } else {
      shell_count.push_back(1);
    }
  }

  //---- Determine locations of opposite parity vectors
  vector<int> pmomentum;
  int check=0;
  for (int i=0; i<momenta_count; i++) {
    pmomentum = momentum[i];

    if (GJP.Xbc()==BOUNDARY_TYPE_PRD) {
      pmomentum[0] = (x_sites-pmomentum[0])%x_sites;
    } else if (GJP.Xbc()==BOUNDARY_TYPE_APRD) {
      pmomentum[0] = (x_sites-pmomentum[0]-1)%x_sites;
    } else {
      ERR.NotImplemented(fname, "Unknown boundary type.");
    }

    if (GJP.Ybc()==BOUNDARY_TYPE_PRD) {
      pmomentum[1] = (y_sites-pmomentum[1])%y_sites;
    } else if (GJP.Ybc()==BOUNDARY_TYPE_APRD) {
      pmomentum[1] = (y_sites-pmomentum[1]-1)%y_sites;
    } else {
      ERR.NotImplemented(fname, "Unknown boundary type.");
    }

    if (GJP.Zbc()==BOUNDARY_TYPE_PRD) {
      pmomentum[2] = (z_sites-pmomentum[2])%z_sites;
    } else if (GJP.Zbc()==BOUNDARY_TYPE_APRD) {
      pmomentum[2] = (z_sites-pmomentum[2]-1)%z_sites;
    } else {
      ERR.NotImplemented(fname, "Unknown boundary type.");
    }

    for (int j=0; j<momenta_count; j++) {
      if (pmomentum==momentum[j]) {
        parity_list[i]=j;
      check=1;
      }
    }
    if (check==1) {
      check = 0;
    } else {
      ERR.General(fname, "Opposite parity vector is missing; try adjusting the fermi energy.");
    }
  }

  //---- Print out some information
  VRB.Flow(fname, "Fermi energy: %f", momenta_arg.fermi_energy);
  VRB.Flow(fname, "Number of momenta below Fermi surface: %d", MomentaCount() );
  VRB.Flow(fname, "Total number of shells below Fermi surface: %d", NumShells() );
  VRB.Flow(fname, "Number of momenta in each shell:");
  for (int i=0; i<NumShells(); i++) {
    VRB.Flow(fname, "Shell (%d/%d): %d", i, NumShells()-1, ShellCount(i));
  }


}

Momenta::~Momenta()
{
  delete[] omega;
  delete[] momentum;
  delete[] parity_list;
}


int Momenta::MomentaCount()
{
  return momenta_count;
}

vector<int> Momenta::GetMomentum(int n)
{

  char* fname = "vector<int> Momenta::GetMomentum(int)";
  VRB.Func(fname);

  if ( (n >= momenta_count)||(n<0) ) {
        ERR.General(fname, "Specified momentum exceeds momentum count.");
  }

  return momentum[n];

}

int Momenta::OppositeParityIndex(int n) { return parity_list[n]; }

int Momenta::ShellCount(int n)
{
  char* fname = "vector<int> Momenta::ShellCount(int)";
  if (n < 0) { ERR.General(fname, "Invalid argument. Shell count starts at zero."); }
  if (n >= NumShells() ) { ERR.General(fname, "Invalid argument. Exceeds number of shells."); }
  return shell_count[n];
}

int Momenta::AccShellCount(int n)
{
  char* fname = "vector<int> Momenta::AccShellCount(int)";
  if (n < 0) { return 0; }
  if (n >= NumShells() ) { ERR.General(fname, "Invalid argument. Exceeds number of shells."); }
  int acc_shell_count = 0;
  for (int s=0; s<=n; s++) {
    acc_shell_count += ShellCount(s);
  }
  return acc_shell_count;
}

int Momenta::NumShells() { return shell_count.size(); }
