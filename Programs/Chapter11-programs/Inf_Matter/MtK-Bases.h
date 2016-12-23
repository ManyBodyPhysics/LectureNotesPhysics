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
//  MtK-Bases.h  --  classes for the s.-p. and 2p1h/2h1p configurations.
//

#ifndef __Matters__MtK_Bases__
#define __Matters__MtK_Bases__



class SpBasisK {
  
public:
  /// Integration bounduaries and number of mesh points (NOT intervals)
  int SpNmax, SpNAlloc;
  int *nx, *ny, *nz, *spin, *chrg, *nsq;
  int *group;
  double *k, *e_kin, *e_HF, *e_sp;
  
  double Lbox, k_Fermi, E_Fermi;
  int ch_min, ch_max, N_holes;
  
  int N_grps;
  int *gr_mlt, *gr_rep;


public:
  
  // functions:
  SpBasisK();
  ~SpBasisK();
  int  Count_sp_basis(int, int, int );
  void Build_sp_basis(int, int, int, int );
  
  
  
  int Build_groups_tables();

  
  
};

class ADC3BasisK {
  
public:
  int *Bas_2p1h, *Bas_2h1p; // pointers to 2p1h/2h1p bases
  int Nbas_2p1h, Nbas_2h1p; // dimensions of the 2p1h/2h1p bases

  int iSpLoc; // to keep the overall k-channel in s.p. basis
  SpBasisK *SpBasLoc;

public:
  
  // functions:
  ADC3BasisK(SpBasisK* );
  ~ADC3BasisK();
  void Count_ALL_2p1h_2h1p_bases( );
  void Count_2p1h_2h1p_basis(int, int*, int* );
  void Build_2p1h_2h1p_basis(int  );
  
  
};


#endif /* defined(__Matters__MtK_Bases__) */
