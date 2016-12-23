//------------------------------------------------------------------------------
//
//   'MattersK' code fo rcalculation of spectral function
//  in nfinite nucleonic matter
//
//  For more details, refer to:
//
//  https://github.com/ManyBodyPhysics/LectureNotesPhysics/blob/master/doc/src/Chapter11-programs/Inf_Matter/READme.txt
//
//  and
//
//  C. Barbieri and A. Carbone, "Self-consistent Green's function approaches",
//  chapter 11 of Lecture Notes in Physics "An advanced course in computational
//  nuclear physics: Bridging the scales from quarks to neutron stars",
//  Edited by Hjorth-Jensen, M. P. Lombardo, U. van Kolck, Editors
//  ISBN:
//  http://arxiv.org/abs/1611.03923
//
//
//  MtInf-Main.cpp
//  Matters
//
//  (c) C. Barbieri and A. Carbone,   Surrey/Darmstadt,   December 2016.
//
//------------------------------------------------------------------------------
//
//  MtK-Potential_minnesota.cpp  --  to generate the matrix elements of the
//                                  interaction in the momentum basis.
//

#include <cmath>
#include <iostream>

using namespace std;


#include "MtK-Physics_constants.h"
#include "MtK-Bases.h"
#include "MtK-Potential_minnesota.h"


#define PI  3.14159265358979323846264338327950288  // first 36 digits (from Num. Recipes)


//double qt2_transf(SpBasisK *bas, int ia, int ib, int ig, int id, double *qu2/*=NULL*/) {
//  
//  int nq_x,nq_y,nq_z;
//  
//  double qt2;
//  
//  qt2  = PI / bas->Lbox;
//  qt2 *= qt2 / 4.0;
//
//  if (NULL != qu2) {
//    nq_x = bas->nx[ia] - bas->nx[id] + bas->nx[ig] - bas->nx[ib];
//    nq_y = bas->ny[ia] - bas->ny[id] + bas->ny[ig] - bas->ny[ib];
//    nq_z = bas->nz[ia] - bas->nz[id] + bas->nz[ig] - bas->nz[ib];
//    
//    (*qu2)  = qt2 * (nq_x*nq_x + nq_y*nq_y + nq_z*nq_z);
//  }
//
//  nq_x = bas->nx[ia] - bas->nx[ig] + bas->nx[id] - bas->nx[ib];
//  nq_y = bas->ny[ia] - bas->ny[ig] + bas->ny[id] - bas->ny[ib];
//  nq_z = bas->nz[ia] - bas->nz[ig] + bas->nz[id] - bas->nz[ib];
//
//  qt2  = qt2 * (nq_x*nq_x + nq_y*nq_y + nq_z*nq_z);
//  
//  return qt2;} //end function qt2_transf


//
//  Constant for the Minnesota potential
//
const double v0r = 200.0;    // MeV
const double v0t = 178.0;    // MeV
const double v0s =  91.85;   // MeV
const double kr  =   1.487;  // fm**-2
const double kt  =   0.639;  // fm**-2
const double ks  =   0.465;  // fm**-2

//
// Matrix elements for the Minnesota potential
//
double V_Minnesota(SpBasisK *bas, int ia, int ib, int ig, int id, double *Vdir/*=NULL*/) {

  int nx1,nx2,nx3,nx4;
  int ny1,ny2,ny3,ny4;
  int nz1,nz2,nz3,nz4;
  int spin1,spin2,spin3,spin4;
  int isospin1,isospin2,isospin3,isospin4;
  
  int nqs_x, nqs_y, nqs_z;//, nqt_x, nqt_y, nqt_z, nqu_x, nqu_y, nqu_z;
  

  nx1 = bas->nx[ia];
  ny1 = bas->ny[ia];
  nz1 = bas->nz[ia];
  spin1 = bas->spin[ia];
  isospin1 = bas->chrg[ia];
  
  nx2 = bas->nx[ib];
  ny2 = bas->ny[ib];
  nz2 = bas->nz[ib];
  spin2 = bas->spin[ib];
  isospin2 = bas->chrg[ib];
  
  nqs_x = nx1 + nx2;
  nqs_y = ny1 + ny2;
  nqs_z = nz1 + nz2;
  
  nx3 = bas->nx[ig];
  ny3 = bas->ny[ig];
  nz3 = bas->nz[ig];
  spin3 = bas->spin[ig];
  isospin3 = bas->chrg[ig];

  nx4 = bas->nx[id];
  ny4 = bas->ny[id];
  nz4 = bas->nz[id];
  spin4 = bas->spin[id];
  isospin4 = bas->chrg[id];

  if ( (nqs_x != nx3 + nx4) || (nqs_y != ny3 + ny4) || (nqs_z != nz3 + nz4) ) return 0.0;

  int nqt_x = nx1 - nx2 + nx4 - nx3;
  int nqt_y = ny1 - ny2 + ny4 - ny3;
  int nqt_z = nz1 - nz2 + nz4 - nz3;
  
  int nqu_x = nx1 - nx2 + nx3 - nx4;
  int nqu_y = ny1 - ny2 + ny3 - ny4;
  int nqu_z = nz1 - nz2 + nz3 - nz4;

  double qt2, qu2;
  double xV1,xV2,xV3,xV4;
  double vr_dir, vs_dir, vt_dir, vr_exc, vs_exc, vt_exc;

  bool Spin13 = (spin1 == spin3) && (spin2 == spin4);
  bool Spin14 = (spin1 == spin4) && (spin2 == spin3);

  bool IsoSp13 = (isospin1 == isospin3) && (isospin2 == isospin4);
  bool IsoSp14 = (isospin1 == isospin4) && (isospin2 == isospin3);

  double x1 = pow( (PI / bas->Lbox) , 2);

  qt2 = x1 * ( nqt_x*nqt_x + nqt_y*nqt_y + nqt_z*nqt_z );
  qu2 = x1 * ( nqu_x*nqu_x + nqu_y*nqu_y + nqu_z*nqu_z );

  

  vr_dir =  ( v0r/pow(bas->Lbox,3) ) * pow(PI/kr, 1.5) * exp(-qt2/kr/4.0);

  vt_dir = -( v0t/pow(bas->Lbox,3) ) * pow(PI/kt, 1.5) * exp(-qt2/kt/4.0);

  vs_dir = -( v0s/pow(bas->Lbox,3) ) * pow(PI/ks, 1.5) * exp(-qt2/ks/4.0);

  vr_exc =  ( v0r/pow(bas->Lbox,3) ) * pow(PI/kr, 1.5) * exp(-qu2/kr/4.0);
  
  vt_exc = -( v0t/pow(bas->Lbox,3) ) * pow(PI/kt, 1.5) * exp(-qu2/kt/4.0);
  
  vs_exc = -( v0s/pow(bas->Lbox,3) ) * pow(PI/ks, 1.5) * exp(-qu2/ks/4.0);

  xV1 = xV2 = xV3 = xV4 = 0.0;
  
  if (Spin13 && IsoSp13) xV1 = 0.5 *(vr_dir + 0.5*vt_dir + 0.5*vs_dir);
  
  if (Spin14 && IsoSp13) xV2 = 0.25*(vt_dir - vs_dir);
  
  if (Spin14 && IsoSp14) xV3 = 0.5 *(vr_dir + 0.5*vt_dir + 0.5*vs_dir);
  
  if (Spin13 && IsoSp14) xV4 = 0.25*(vt_dir - vs_dir);
  
  x1 = xV1+xV2-xV3-xV4;

  if (NULL != Vdir) (*Vdir) = x1;

  xV1 = xV2 = xV3 = xV4 = 0.0;

  if (Spin14 && IsoSp14) xV1 = 0.5 *(vr_exc + 0.5*vt_exc + 0.5*vs_exc);
  
  if (Spin13 && IsoSp14) xV2 = 0.25*(vt_exc - vs_exc);
  
  if (Spin13 && IsoSp13) xV3 = 0.5 *(vr_exc + 0.5*vt_exc + 0.5*vs_exc);
  
  if (Spin14 && IsoSp13) xV4 = 0.25*(vt_exc - vs_exc);

  
  x1 -= (xV1+xV2-xV3-xV4);


  return x1;
}  // END FUNCTION V_Minnesota


//
//  Matrix elements of the HF **Hamiltonian** (i.e. Tkin+Vhf), assuming the Minnesota potential
double HF_Minn_MtxEls(SpBasisK *bas, int ia, int ib, int Nocc) {
  
  double VHF;
  
  VHF = 0.0;
  
  for(int ig=0; ig<Nocc; ++ig) VHF += V_Minnesota(bas,ia,ig,ib,ig);

  //cout << "   HF:  "<< VHF << "   ";

  if (ia == ib) VHF += bas->e_kin[ia];
  
  //cout << "   "<< VHF << "   \n";
  

  return VHF;} // end function HF_Minn_MtxEls
