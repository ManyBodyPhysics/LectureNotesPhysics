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
//  MtK-Physics_constants.h  --  global definitions of relvent physics constants.
//



#ifndef __Matters__Phys_consts__
#define __Matters__Phys_consts__

#define PI  3.14159265358979323846264338327950288  // first 36 digits (from Num. Recipes)

static double MeVfm = 197.326968; // hbar*c

//
//
//  Values also used in arXiv:1502.04682 [nucl-th] (supp. material):
//    r_p^2 =  0.8775 fm    [Rev. Mod. Phys. 84, 1527 (2012)]
//    r_n^2 = −0.1149 fm^2  [At. Data Nucl. Data Tables 99, 69 (2013)]
#define PROTON_CH_RADIUS        0.8775
#define NEUTRON_ABS_CH_RADIUS   0.338969025133566
//


//
// Usual BoccaDorata  values
//
//#define MeVfm       197.326968
#define ELECTRONmass   0.511
//#define NUCLEONmass  938.9182491198266689
#define PROTONmass   938.9182491198266689
#define NEUTRONmass  938.9182491198266689

#define NUCLEONmass  939.565


//
//  CODATA Internationally recommended 2010 values, taken
// from NIST website [averaged nucleon mass is calculated
// as 2*mp*mn/(mp+mn)].
//
//#define MeVfm       197.3269718121
//#define ELECTRONmass   0.510998928
//#define PROTONmass   938.272046
//#define NEUTRONmass  939.565379
//#define NUCLEONmass  938.91826711792735
//#define PROTON_CH_RADIUS    0.8775
//#define NEUTRON_ABS_CH_RADIUS   ??
//

//
// Used in the FaddO16 code
//
//#define MeVfm 197.327053
//#define NUCLEONmass 938.65

//
// Values used in the Vucom package
//
//#define MeVfm       197.326968
//#define NUCLEONmass 938.9182491198266689

//
// Values used in Hjorth-Jensen codes
//
//#define MeVfm       197.3269680000000000
//#define NUCLEONmass 938.9260000000000000

#endif /* defined(__Matters__Phys_consts__) */
