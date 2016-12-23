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
//  MtK-Potential_minnesota.h  --  to generate the matrix elements of the
//                                  interaction in the momentum basis.
//

#ifndef ____MtK_Potential_minnesota__
#define ____MtK_Potential_minnesota__

#include <stdio.h>

//double qt2_transf(SpBasisK*, int, int, int, int, double *qu2=NULL); //end function qt2_transf

double V_Minnesota(SpBasisK*, int, int, int, int, double *Vdir=NULL);



double HF_Minn_MtxEls(SpBasisK*, int , int , int );


#endif /* defined(____MtK_Potential_minnesota__) */
