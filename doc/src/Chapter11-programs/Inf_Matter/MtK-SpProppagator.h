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
//  MtK-SpProppagator.h  --  class to store the spectral (or Lehmann) representation; this
//                            is meant for *both* the sp proppagator and the self-energy.
//

#ifndef __Matters__MtK_SpctDist__
#define __Matters__MtK_SpctDist__



class SpctDist {
  
public:
  SpBasisK *SpBasLoc;  // associated s.p. basis

  int Ndim_Leh, N_LEH_ALLOC; // number of k values from the s.. basis
  
  int *i_sp, *i_grp, *i_grp_mult;
  int *nsq;
  double *t_kin, *e_sp, *Sig_inf;
  double *Tot_fw_str, *Tot_bk_str;
  
  int     *N_fw_pls,     *N_bk_pls;
  int     *N_PLS_ALLOC;
  double **ek_fw,       **ek_bk; // poles of the propagator/self-energy
  double **Ak_fw,       **Ak_bk; // this is the FULL residue (i.e. the SF, not the amplitude X,Y)
  double **Bk_fw,       **Bk_bk; // Bk_xx  are not there for the propagator but are  used to store the tridiagonal
                                 // line from the Lanczos when this class is used to store C/D  and M/N matrices.

  double  EFermi, Atot, EKoltun, Tkin;  // expectation values


public:
  
  // functions:
  SpctDist(SpBasisK* );
  ~SpctDist();
  void Initialize(void );

  int free_mem(); // needed by the destructor and to 'allocate_k_kchannel' to store a new propagator
  int add_k_channel(int, int, double*, double*, int, double*, double*,
                    double in_Sig_inf=0.0, double *B_fw_in=NULL, double *B_bk_in=NULL);
  int allocate_k_kchannel(int, int, int);
  
  int Make_HF_propagator(void ); // Initialises to a HF propagator based on the current basis 'SpBasLoc'


  double Get_EFermi(double *xeF_h=NULL, double *xeF_p=NULL,
                    double *xe_mn=NULL, double *xe_mx=NULL  );

  double Set_EFermi(double );
  double Seek_Fermi_level(double );

  double Koltun(void );
                              
  int Get_Leh_index(int );

  double Make_LorSF_vs_E(int, double**, double**, double**, int*, bool InvIm=false);
  int Plot_SpectFnct_vs_E(int );
  int Plot_SelfEn_vs_E(int );
  int Plot_SpectFnct_3D(void );
  int DiagTridiagSE(int );

};


#endif /* defined(__Matters__MtK_SpctDist__) */
