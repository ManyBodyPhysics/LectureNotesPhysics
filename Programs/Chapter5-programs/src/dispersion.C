#include <iostream>
#include <limits>
#include <math.h>
#include <stdlib.h>
using namespace std;
#include "dispersion.h"
#include "error.h"
#include "global_job_parameter.h"
#include "constants.h"

Dispersion::Dispersion(DispersionType dispersion_type, double mass, CutoffType cutoff_type)
{

  char* fname = "Dispersion::Dispersion(DispersionType)";

  dispersion = (double *) malloc(GJP.Vol()*sizeof(double));  

  double bc[3];
  GJP.GetBcShift(bc);

  int x_sites = GJP.Xsites();
  int y_sites = GJP.Ysites();
  int z_sites = GJP.Zsites();

  double sx;
  double sy;
  double sz;

  double psq;

  int z;
  int Y;
  int y;
  int x;

  for(int i=0; i<GJP.Vol(); i++) {
 
    z = i%z_sites;
    Y = i/z_sites;
    y = Y%y_sites;
    x = Y/y_sites;


    if (dispersion_type==DISPERSION_TYPE_STANDARD) {
      //
      // \Delta(p) = 2 \sum_j \sin^2(p_j/2)
      //
      sx = sin( (x+bc[0]) * PI / x_sites );
      sy = sin( (y+bc[1]) * PI / y_sites );
      sz = sin( (z+bc[2]) * PI / z_sites );
      dispersion[i] = 2.0 * ( sx*sx + sy*sy + sz*sz ); 
      dispersion[i] /= mass; 
    }
  
    if ( (dispersion_type==DISPERSION_TYPE_QUADRATIC)||(dispersion_type==DISPERSION_TYPE_PERFECT) ) {

      if (x >= x_sites/2) {  
        switch ( GJP.Xbc() ) {
          case BOUNDARY_TYPE_APRD:
            x = x_sites - x - 1;
            break;
          case BOUNDARY_TYPE_PRD:
            x = x_sites - x;
            break;
          default:
            ERR.General(fname,"Invalid x-BOUNDARY_TYPE.");
        }
      }

      if (y >= y_sites/2) {  
        switch ( GJP.Ybc() ) {
          case BOUNDARY_TYPE_APRD:
            y = y_sites - y - 1;
            break;
          case BOUNDARY_TYPE_PRD:
            y = y_sites - y;
            break;
          default:
            ERR.General(fname,"Invalid y-BOUNDARY_TYPE.");
        }
      }

      if (z >= z_sites/2) {  
        switch ( GJP.Zbc() ) {
          case BOUNDARY_TYPE_APRD:
            z = z_sites - z - 1;
            break;
          case BOUNDARY_TYPE_PRD:
            z = z_sites - z;
            break;
          default:
            ERR.General(fname,"Invalid z-BOUNDARY_TYPE.");
        }
      }

      sx = (x + bc[0]) / (double)x_sites;
      sy = (y + bc[1]) / (double)y_sites;
      sz = (z + bc[2]) / (double)z_sites;
      psq = (sx*sx + sy*sy + sz*sz)*4.0*PISQ;

      //---- Cutoff, beyond which we take psq = infinity
      if ( (cutoff_type==CUTOFF_TYPE_HARD)&&(psq > GJP.CutoffSq() ) ) {
        dispersion[i] = numeric_limits<double>::infinity();
      } else {
        if (dispersion_type==DISPERSION_TYPE_QUADRATIC) { dispersion[i] = psq/(2.0*mass); }
        if (dispersion_type==DISPERSION_TYPE_PERFECT) { dispersion[i] = ( exp( psq/(2.0*mass) ) -1.0 ); }
      }



    }

  }

}

Dispersion::~Dispersion()
{
  free(dispersion);
}


double Dispersion::Get(int n1, int n2, int n3)
{

  char* fname = "double Dispersion::Get(int, int, int)";
  return dispersion[n3+n2*GJP.Zsites()+n1*GJP.Ysites()*GJP.Zsites()];

}

double Dispersion::Get(int i)
{

  char* fname = "double Dispersion::Get(int)";

  return dispersion[i];
}
