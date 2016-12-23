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
//  MtK-Lanczos.cpp  --  Lanczos algorithm to reduce the 'C' and 'D' submatrices
//                      in the Dyson eigenvalue problem.
//


#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

#include "MtK-Bases.h"
#include "MtK-Potential_minnesota.h"
#include "MtK-Global_data.h"


// BLAS routines
//---------------

extern "C" void dgemv_(char*, int*, int*,    double*, double*, int*,
                       double*, int*, double*,
                       double*, int*);
//DGEMV  performs one of the matrix-vector operations
//
//y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
//
//where alpha and beta are scalars, x and y are vectors and A is an
//m by n matrix.
//
//subroutine dgemv (character TRANS, integer	M, integer N,
//                  double precision ALPHA, double precision, dimension(lda,*) A, integer LDA,
//                  double precision, dimension(*) X, integer INCX, double precision  	BETA,
//                  double precision, dimension(*) Y, integer INCY)


extern "C" void dgemm_(char*, char*,   int*, int*, int*,    double*, double*, int*,
                                          double*, int*,    double*, double*, int*);
//DGEMM  performs one of the matrix-matrix operations
//
//C := alpha*op( A )*op( B ) + beta*C,
//
//where  op( X ) is one of
//
//op( X ) = X   or   op( X ) = X**T,
//
//alpha and beta are scalars, and A, B and C are matrices, with op( A )
//an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
//
//subroutine dgemm 	(character TRANSA, character TRANSB,   integer M, integer N, integer K,
//                   double precision ALPHA, double precision, dimension(lda,*)	A, integer LDA,
//                                           double precision, dimension(ldb,*) B, integer LDB,
//                   double precision	BETA,  double precision, dimension(ldc,*) C, integer LDC)



int Build_DysonADC_LancMatrix(ADC3BasisK *ADC3Bas_in, int i_fwbk, int NLanc,
                                             double *Mst, double *Cst, int LDC ) {
  
  SpBasisK *WK_bas = ADC3Bas_in->SpBasLoc;
  
  int isp = ADC3Bas_in->iSpLoc;
  int Ndim;

  double *MNmtx = NULL;
  double *CDmtx = NULL;
  double *PVmtx = NULL;
  
  if (i_fwbk > 0) {
    //
    Ndim = ADC3Bas_in->Nbas_2p1h;
    //
    MNmtx = new double[Ndim];
    CDmtx = new double[Ndim*Ndim];
    //
    Build_DysonADC_Matrix_2p1h(ADC3Bas_in, MNmtx, CDmtx, Ndim);
    //
  } else if (i_fwbk < 0) {
    //
    Ndim = ADC3Bas_in->Nbas_2h1p;
    //
    MNmtx = new double[Ndim];
    CDmtx = new double[Ndim*Ndim];
    //
    Build_DysonADC_Matrix_2h1p(ADC3Bas_in, MNmtx, CDmtx, Ndim);
    //
  } else {
    //
    cerr << "\n\n WARNING Build_DysonADC_LancMatrix is returning without doing anything, i_fwbk = "<<i_fwbk<<endl<<endl<<flush;
  }
  
  if (NLanc > LDC) {
    cerr << "\n\n WARNING Build_DysonADC_LancMatrix is lowering NLanc from "<<NLanc
         << " to  " << LDC <<" (==LDC); called with i_fwbk = "<<i_fwbk<<endl<<endl<<flush;
    NLanc = LDC;
  }

  if (NLanc > Ndim) {
    cerr << "\n\n WARNING Build_DysonADC_LancMatrix is lowering NLanc from "<<NLanc
         << " to  " << int(Ndim * .98) <<" (Ndim="<<Ndim<<"); called with i_fwbk = "
         << i_fwbk<<endl<<endl<<flush;
    NLanc = int(Ndim * .98);
  }

//  for (int i1=0; i1<NLanc; ++i1) {
//    for (int i2=0; i2<NLanc; ++i2) Cst[i1*LDC + i2] =CDmtx[i1*Ndim+i2];;
//    Mst[i1] = MNmtx[i1];
//  }
//
//  if (NULL != MNmtx) delete [] MNmtx;
//  if (NULL != CDmtx) delete [] CDmtx;
//  if (NULL != PVmtx) delete [] PVmtx;
//
//return 0;
  
  double x1, x2, x3;
  double alpha, beta;

  char No_trans = 'N';
  double   Zero = 0.0, One =1.0;
  int    IntOne = 1;

  
  PVmtx = new double[(NLanc+1)*Ndim];
  for (int i1=0; i1<(NLanc+1)*Ndim; ++i1) PVmtx[i1] = 0.0;
  
  for (int i1=0; i1<NLanc; ++i1) {
    for (int i2=0; i2<NLanc; ++i2) Cst[i1*LDC + i2] =0.0;
    Cst[i1*LDC + i1] =1.e6 * i_fwbk;
    Mst[i1] =0.0;
  }

  int     i_lncz  = 0;
  double *dp_vec  = PVmtx + Ndim*i_lncz;
  double *dp_vec2;


//x1 = 0.0;
  for (int i1=0; i1<Ndim; ++i1) {dp_vec[i1] = MNmtx[i1];}// x1 += (dp_vec[i1]*dp_vec[i1]);}
//x1 = sqrt(x1);
//for (int i1=0; i1<Ndim; ++i1)  dp_vec[i1] /= x1;


  for(int i_lncz=0; i_lncz<NLanc; ++i_lncz) {

    dp_vec2 = PVmtx + Ndim*(i_lncz+1);
    dp_vec  = PVmtx + Ndim*(i_lncz);


    //
    // Perform matrix x vector multiplication:
    //
    //    dp_vec2 = CDmtx x dp_vec
    //
    dgemv_(&No_trans, &Ndim,  &Ndim,  &One,   CDmtx,   &Ndim,
                                              dp_vec,  &IntOne,
                                      &Zero,  dp_vec2, &IntOne);
    //
    // Alternative way with BLAS3:
    //dgemm_(&No_trans, &No_trans, &Ndim, &IntOne, &Ndim, &One,   CDmtx,  &Ndim,
    //                                                           dp_vec,  &Ndim,
    //                                                   &Zero,  dp_vec2, &Ndim);


//    for (int i1=0; i1<Ndim; ++i1) {
//      x1 = 0.0;
//      for (int i2=0; i2<Ndim; ++i2) x1 += CDmtx[i2*Ndim+i1] * dp_vec[i2];
//      dp_vec2[i1] = x1;
//    }

    x1 = 0.0; x2 = 0.0; x3 = 0.0;
    for(int i1=0; i1<Ndim; ++i1) {
      x3 += dp_vec2[i1]*dp_vec[i1];
      x1 += dp_vec[i1]*dp_vec[i1];
      x2 += dp_vec[i1]*MNmtx[i1];
    }
  //alpha[i_lncz]  = x3 / x1;
  //beta[i_lncz]   = sqrt(x1);
    alpha = x3 / x1;
    beta  = sqrt(x1);
    Mst[i_lncz]  = x2 / beta;
    Cst[(LDC+1)*i_lncz] = alpha;
    if (0 < i_lncz) {
      Cst[(LDC+1)*i_lncz - 1 ]   = beta;
      Cst[(LDC+1)*i_lncz - LDC ] = beta;
    }
    //cout << "\n x's, a, b = "<<x1 <<" ,  "<<x2 <<" ,  "<<x3 <<" ,  "<<alpha <<" ,  "<<beta;
    for (int i1=0; i1<Ndim; ++i1) {
      dp_vec[i1]  /= beta;
      dp_vec2[i1] /= beta;
      dp_vec2[i1] -= alpha*dp_vec[i1];
    }
    if (0 < i_lncz) for (int i1=0; i1<Ndim; ++i1) dp_vec2[i1] -= beta*dp_vec[i1-Ndim];


  
  //G-S
    for(int i2=i_lncz-2; i2 >= 0; --i2) {
      x1=0.0; x2=0.0;
      dp_vec =  PVmtx + Ndim*i2;
      for (int i1=0; i1<Ndim; ++i1) {x1 += dp_vec2[i1]*dp_vec[i1]; x2 += dp_vec[i1]*dp_vec[i1];}
      x1 /= x2;
      for (int i1=0; i1<Ndim; ++i1) {dp_vec2[i1] -= x1*dp_vec[i1];}
    }

  }

  
  
  
//  for(int i_lncz=0; i_lncz<NLanc; ++i_lncz) {
//    
//    dp_vec2 = PVmtx + Ndim*(i_lncz+1);
//    dp_vec  = PVmtx + Ndim*(i_lncz);
//    
// 
//    x1 = 0.0;
//    for(int i1=0; i1<Ndim; ++i1) x1 += dp_vec[i1]*dp_vec[i1];
//    x1 = sqrt(x1);
//    for(int i1=0; i1<Ndim; ++i1) dp_vec[i1] /= x1;
//
//    //
//    // Perform matrix x vector multiplication:
//    //
//    //    dp_vec2 = CDmtx x dp_vec
//    //
//    dgemv_(&No_trans, &Ndim,  &Ndim,  &One,   CDmtx,   &Ndim,
//           dp_vec,  &IntOne,
//           &Zero,  dp_vec2, &IntOne);
//    //
//    // Alternative way with BLAS3:
//    //dgemm_(&No_trans, &No_trans, &Ndim, &IntOne, &Ndim, &One,   CDmtx,  &Ndim,
//    //                                                           dp_vec,  &Ndim,
//    //                                                   &Zero,  dp_vec2, &Ndim);
//    
//    
//    //    for (int i1=0; i1<Ndim; ++i1) {
//    //      x1 = 0.0;
//    //      for (int i2=0; i2<Ndim; ++i2) x1 += CDmtx[i2*Ndim+i1] * dp_vec[i2];
//    //      dp_vec2[i1] = x1;
//    //    }
//
//    alpha = 0.0;
//    for(int i1=0; i1<Ndim; ++i1) alpha += dp_vec2[i1]*dp_vec[i1];
//    
//    beta = 0.0;
//    if (0 < i_lncz) for(int i1=0; i1<Ndim; ++i1) beta +=  dp_vec2[i1]*dp_vec[i1-Ndim];
//
//    x2 = 0.0;
//    for(int i1=0; i1<Ndim; ++i1) x2 += MNmtx[i1]*dp_vec[i1];
//
//    Mst[i_lncz]  = x2;
//    Cst[(LDC+1)*i_lncz] = alpha;
//    if (0 < i_lncz) {
//      Cst[(LDC+1)*i_lncz - 1 ]   = beta;
//      Cst[(LDC+1)*i_lncz - LDC ] = beta;
//    }
//    
//    //G-S
//    for(int i2=i_lncz; i2 >= 0; --i2) {
//      x1=0.0; x2=0.0;
//      dp_vec =  PVmtx + Ndim*i2;
//      for (int i1=0; i1<Ndim; ++i1) {x1 += dp_vec2[i1]*dp_vec[i1]; x2 += dp_vec[i1]*dp_vec[i1];}
//      x1 /= x2;
//      for (int i1=0; i1<Ndim; ++i1) {dp_vec2[i1] -= x1*dp_vec[i1];}
//    }
//    
//  }

  if (NULL != MNmtx) delete [] MNmtx;
  if (NULL != CDmtx) delete [] CDmtx;
  if (NULL != PVmtx) delete [] PVmtx;

  return 0;}

