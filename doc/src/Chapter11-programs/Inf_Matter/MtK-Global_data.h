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
//  MtK-Physics_constants.h  --  global variables to be used throught the code.
//

#ifndef ____MtK_Global_data__
#define ____MtK_Global_data__


#define PI  3.14159265358979323846264338327950288  // first 36 digits (from Num. Recipes)

#include <stdio.h>

#include "MtK-Bases.h"
#include "MtK-SpProppagator.h"
//#include "MtK-Potential_minnesota.h"

extern int NLanczos; // = -100;

extern double GammaPlot_1, GammaPlot_2;



extern SpctDist  *g_sp_out, *Sigma_irred; // Global pointers to the dressed prop. and Self energy

//
//  Dyson/Lanczos stuff:
//
int Solve_Dyson(    ADC3BasisK*,      double*, double*, double*);
int Solve_Dyson_sp( ADC3BasisK*, int, double*, double*, double*, SpctDist *gsp_in=NULL);
int Solve_Dyson_cHF(double,      int, double*, double*, double*, SpctDist*, SpctDist*, SpctDist*);

int Build_DysonADC_Matrix_2p1h(ADC3BasisK*, double*, double*, int);
int Build_DysonADC_Matrix_2h1p(ADC3BasisK*, double*, double*, int);

int Build_DysonADC_LancMatrix(ADC3BasisK*,  int, int, double*, double*, int);


// Generate datafile with the spectral distribution in the HF approximation
int Plot_HF_dist(SpBasisK* );


//
// Sorting routines (in MtK-Utilities.cpp):
//
extern void Sort_u_double2dim(int , int , int ,double []);

extern void Sort_d_double2dim(int , int , int ,double []);


#endif /* defined(____MtK_Global_data__) */
