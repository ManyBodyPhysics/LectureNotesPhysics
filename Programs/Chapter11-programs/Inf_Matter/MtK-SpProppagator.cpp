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
//  MtK-SpProppagator.cpp  --  class to store the spectral (or Lehmann) representation; this
//                            is meant for *both* the sp proppagator and the self-energy.
//


#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

#include "MtK-Global_data.h"
#include "MtK-Bases.h"
#include "MtK-SpProppagator.h"
#include "MtK-Physics_constants.h"


// LAPACK driver for diagonalizaton
extern "C" void dsyevd_(char*, char*, int*, double*, int*,
                        double*, double*, int*, int*, int*,  int* );

static char Vectors = 'V';
static char Upper   = 'U';



///
///   Creators and destructors
///

// The rationale of Init_NULL and Free_mem si that
// they are written only once and called when needed
// so changes are always done in one place only.

SpctDist::SpctDist(SpBasisK *bas_in) {
  
  SpBasLoc = bas_in;
  this->Initialize();

  return;}

void SpctDist::Initialize(void ) {
    
  //
  //this->Init_NULL();
  Ndim_Leh    = -100;
  N_LEH_ALLOC = -100;
  //
  i_sp       = NULL;
  i_grp      = NULL;
  i_grp_mult = NULL;
  //
  nsq     = NULL;
  t_kin   = NULL;
  e_sp    = NULL;
  Sig_inf = NULL;
  EFermi = 0.0;
  Atot = 0.0;
  EKoltun = 0.0;
  Tkin = 0.0;
  //
  Tot_fw_str = NULL;    Tot_bk_str = NULL;
  //
  N_fw_pls   = NULL;    N_bk_pls   = NULL;
  N_PLS_ALLOC = NULL;
  ek_fw      = NULL;    ek_bk      = NULL;
  Ak_fw      = NULL;    Ak_bk      = NULL;
  Bk_fw      = NULL;    Bk_bk      = NULL;
  //

  //
  // For now, always initialise as many slots as there are groups:
  //
  Ndim_Leh    = 0;
  N_LEH_ALLOC = this->SpBasLoc->N_grps;
  i_sp       = new int[N_LEH_ALLOC];
  i_grp      = new int[N_LEH_ALLOC];
  i_grp_mult = new int[N_LEH_ALLOC];
  //
  nsq     = new int[N_LEH_ALLOC];
  t_kin   = new double[N_LEH_ALLOC];
  e_sp    = new double[N_LEH_ALLOC];
  Sig_inf = new double[N_LEH_ALLOC];
  //
  Tot_fw_str = new double[N_LEH_ALLOC];     Tot_bk_str = new double[N_LEH_ALLOC];
  //
  N_fw_pls   = new int[N_LEH_ALLOC];        N_bk_pls   = new int[N_LEH_ALLOC];
  N_PLS_ALLOC = new int[N_LEH_ALLOC];
  ek_fw      = new double*[N_LEH_ALLOC];    ek_bk      = new double*[N_LEH_ALLOC];
  Ak_fw      = new double*[N_LEH_ALLOC];    Ak_bk      = new double*[N_LEH_ALLOC];
  Bk_fw      = new double*[N_LEH_ALLOC];    Bk_bk      = new double*[N_LEH_ALLOC];

  for (int isp=0; isp<N_LEH_ALLOC; ++isp) {
    i_sp      [isp] = -100;
    i_grp     [isp] = -100;
    i_grp_mult[isp] = -100;
    //
    nsq    [isp] = -100;
    t_kin  [isp] = 0.0;
    e_sp   [isp] = 0.0;
    Sig_inf[isp] = 0.0;
    //
    Tot_fw_str[isp] = 0.0;     Tot_bk_str[isp] = 0.0;
    //
    N_fw_pls  [isp] = -100;    N_bk_pls  [isp] = -100;
    N_PLS_ALLOC[isp] = -100;
    ek_fw[isp]      = NULL;    ek_bk[isp]      = NULL;
    Ak_fw[isp]      = NULL;    Ak_bk[isp]      = NULL;
    Bk_fw[isp]      = NULL;    Bk_bk[isp]      = NULL;
 }
  
  return;}

SpctDist::~SpctDist(void ) {
  // Must free all the allocated memory
  //  this->Free_mem();
  for (int i_all=0; i_all<N_LEH_ALLOC; ++i_all) {
    if (NULL != ek_bk[i_all]) delete [] ek_bk[i_all];  ek_bk[i_all] = NULL;  ek_fw[i_all] = NULL;
    if (NULL != Ak_bk[i_all]) delete [] Ak_bk[i_all];  Ak_bk[i_all] = NULL;  Ak_fw[i_all] = NULL;
    if (NULL != Bk_bk[i_all]) delete [] Bk_bk[i_all];  Bk_bk[i_all] = NULL;  Bk_fw[i_all] = NULL;
  }
  //
  delete [] i_sp;
  delete [] i_grp;
  delete [] i_grp_mult;
  //
  delete [] nsq;
  delete [] t_kin;
  delete [] e_sp;
  delete [] Sig_inf;
  //
  delete [] Tot_fw_str;    delete [] Tot_bk_str;
  //
  delete [] N_fw_pls;      delete [] N_bk_pls;
  delete [] N_PLS_ALLOC;
  delete [] ek_fw;         delete [] ek_bk;
  delete [] Ak_fw;         delete [] Ak_bk;
  delete [] Bk_fw;         delete [] Bk_bk;


  //  this->Init_NULL();
  Ndim_Leh    = -100;
  N_LEH_ALLOC = -100;
  //
  i_sp       = NULL;
  i_grp      = NULL;
  i_grp_mult = NULL;
  //
  nsq     = NULL;
  t_kin   = NULL;
  e_sp    = NULL;
  Sig_inf = NULL;
  EFermi = 0.0;
  Atot = 0.0;
  EKoltun = 0.0;
  Tkin = 0.0;
  //
  Tot_fw_str = NULL;    Tot_bk_str = NULL;
  //
  N_fw_pls   = NULL;    N_bk_pls   = NULL;
  N_PLS_ALLOC = NULL;
  ek_fw      = NULL;    ek_bk      = NULL;
  Ak_fw      = NULL;    Ak_bk      = NULL;
  Bk_fw      = NULL;    Bk_bk      = NULL;

  return;}




//
//
//

int SpctDist::add_k_channel(int isp_in,
                                 int N_fw_in, double *A_fw_in, double *E_fw_in,
                                 int N_bk_in, double *A_bk_in, double *E_bk_in,
                                 double in_Sig_inf /*=0.0*/,
                                 double *B_fw_in/*=NULL*/, double *B_bk_in/*=NULL*/){

  int i_Leh_loc = -100;
  for (int i1=0; i1<Ndim_Leh; ++i1)
    if (isp_in == i_sp[i1]) {i_Leh_loc = i1; break;}

  //cout << " i_Leh_loc = "<<i_Leh_loc<<"   "<<Ndim_Leh<<"  *   "<<isp_in <<"   "<<N_fw_in<<"   "<<N_bk_in;
  
  bool allocate_new = false;
       if  (0 > i_Leh_loc) {allocate_new = true;}
  else if ( N_fw_in + N_bk_in > N_PLS_ALLOC[i_Leh_loc] ) {allocate_new = true;}
  //
  if (allocate_new) allocate_k_kchannel(isp_in, N_fw_in+10, N_bk_in);

  i_Leh_loc = -100;
  for (int i1=0; i1<Ndim_Leh; ++i1)
    if (isp_in == i_sp[i1]) {i_Leh_loc = i1; break;}

  //cout << "  ---   i_Leh_loc = "<<i_Leh_loc<<"   "<<Ndim_Leh<<"  *   "<<isp_in<<endl<< flush;

  if  (0 > i_Leh_loc) {cerr << "\n Something really wrong happened here !!! \n\n"; exit(100);}

  double x1 = 0.0;
  N_bk_pls[i_Leh_loc] = N_bk_in;
  for (int ibk=0; ibk<N_bk_in; ++ibk) {
    ek_bk[i_Leh_loc][ibk] = E_bk_in[ibk];
    Ak_bk[i_Leh_loc][ibk] = A_bk_in[ibk];
    Bk_bk[i_Leh_loc][ibk] = 0.0;
    if (NULL != B_bk_in) Bk_bk[i_Leh_loc][ibk] = B_bk_in[ibk];
    x1 += A_bk_in[ibk]*A_bk_in[ibk];
  }
  Tot_bk_str[i_Leh_loc] = x1;
  // sort ek_bk ,  Ak_bk  and  Bk_bk
  
  ek_fw[i_Leh_loc] = ek_bk[i_Leh_loc] + N_bk_in;
  Ak_fw[i_Leh_loc] = Ak_bk[i_Leh_loc] + N_bk_in;
  Bk_fw[i_Leh_loc] = Bk_bk[i_Leh_loc] + N_bk_in;
  //
  x1 = 0.0;
  N_fw_pls[i_Leh_loc] = N_fw_in;
  for (int ifw=0; ifw<N_fw_in; ++ifw) {
    ek_fw[i_Leh_loc][ifw] = E_fw_in[ifw];
    Ak_fw[i_Leh_loc][ifw] = A_fw_in[ifw];
    Bk_fw[i_Leh_loc][ifw] = 0.0;
    if (NULL != B_fw_in) Bk_fw[i_Leh_loc][ifw] = B_fw_in[ifw];
    x1 += A_fw_in[ifw]*A_fw_in[ifw];
  }
  Tot_fw_str[i_Leh_loc] = x1;
  // sort ek_fw ,  Ak_fw  and  Bk_fw

  Sig_inf[i_Leh_loc] = in_Sig_inf;

  return 0;}


int SpctDist::allocate_k_kchannel(int isp_in, int N_fw_max, int N_bk_max) {

  
  if (N_fw_max < 1) N_fw_max = 1;
  if (N_bk_max < 1) N_bk_max = 1;


  int i_Leh = Ndim_Leh;
  for (int i1=0; i1<Ndim_Leh; ++i1)
    if (isp_in == i_sp[i1]) {
      //
      cerr << "\n\n WARNING (SpctDist::allocate_k_kchannel): channel i_sp = "<< isp_in
           <<" will be deleted and reallocated\n\n"<<flush;
      i_Leh = i1;
      --Ndim_Leh;  // because it will be increased at the end of the function
      delete [] ek_bk[i_Leh];
      delete [] Ak_bk[i_Leh];
      delete [] Bk_bk[i_Leh];
      break;
  }

  if (Ndim_Leh >= N_LEH_ALLOC) {
    cerr << "\n\n ERROR (SpctDist::allocate_k_kchannel): cannot allocate any new channels!!! \n\n"<<flush;
  }
  
  i_sp      [i_Leh] = isp_in;
  i_grp     [i_Leh] = this->SpBasLoc->group[isp_in];
  i_grp_mult[i_Leh] = this->SpBasLoc->gr_mlt[this->SpBasLoc->group[isp_in]];
  //
  nsq    [i_Leh] = this->SpBasLoc->nsq[isp_in];
  t_kin  [i_Leh] = this->SpBasLoc->e_kin[isp_in];
  e_sp   [i_Leh] = this->SpBasLoc->e_sp[isp_in];
  Sig_inf[i_Leh] = 0.0;
  //
  Tot_fw_str[i_Leh] = 0.0;         Tot_bk_str[i_Leh] = 0.0;
  //
  N_fw_pls  [i_Leh] = -100;        N_bk_pls  [i_Leh] = -100;
  N_PLS_ALLOC[i_Leh] = N_bk_max + N_fw_max;
  ek_bk[i_Leh]      = new double[N_PLS_ALLOC[i_Leh]];  ek_fw[i_Leh] = ek_bk[i_Leh] + N_bk_max;
  Ak_bk[i_Leh]      = new double[N_PLS_ALLOC[i_Leh]];  Ak_fw[i_Leh] = Ak_bk[i_Leh] + N_bk_max;
  Bk_bk[i_Leh]      = new double[N_PLS_ALLOC[i_Leh]];  Bk_fw[i_Leh] = Bk_bk[i_Leh] + N_bk_max;

  for (int ibk=0; ibk<N_PLS_ALLOC[i_Leh]; ++ibk) {
    ek_bk[i_Leh][ibk] = -1.e6;
    Ak_bk[i_Leh][ibk] = 0.0;
    Bk_bk[i_Leh][ibk] = 0.0;
  }

  ++Ndim_Leh;
  
  return 0;}



int SpctDist::Get_Leh_index(int i_ksp) {

  for (int i_Leh=0; i_Leh<Ndim_Leh; ++i_Leh)
    if (this->i_grp[i_Leh] == this->SpBasLoc->group[i_ksp]) return i_Leh;

  // If not found, the return fail:
  return -100;}


double SpctDist::Get_EFermi(double *xeF_h/*=NULL*/, double *xeF_p/*=NULL*/,
                            double *xe_mn/*=NULL*/, double *xe_mx/*=NULL*/ ) {
 
  double ye_mn=1.e8, yeF_h=-1.e8, yeF_p=1.e8, ye_mx=-1.e8;

  if (Ndim_Leh > 0) {
    ye_mn = yeF_h = ek_bk[0][0];
    ye_mx = yeF_p = ek_fw[0][0];
  }

  for (int iksp=0; iksp<Ndim_Leh; ++iksp) {
    for (int i=0; i<N_bk_pls[iksp]; ++i) {
      ye_mn = (ye_mn < ek_bk[iksp][i]) ? ye_mn : ek_bk[iksp][i];
      yeF_h = (yeF_h > ek_bk[iksp][i]) ? yeF_h : ek_bk[iksp][i];
    }
    for (int i=0; i<N_fw_pls[iksp]; ++i) {
      ye_mx = (ye_mx > ek_fw[iksp][i]) ? ye_mx : ek_fw[iksp][i];
      yeF_p = (yeF_p < ek_fw[iksp][i]) ? yeF_p : ek_fw[iksp][i];
    }
  }

  if (NULL != xe_mn) *xe_mn = ye_mn;
  if (NULL != xe_mx) *xe_mx = ye_mx;
  if (NULL != xeF_h) *xeF_h = yeF_h;
  if (NULL != xeF_p) *xeF_p = yeF_p;
  
  return (yeF_h + yeF_p)/2.0;
  }


int SpctDist::DiagTridiagSE(int i_loc_Leh) {
  
  int LDA = N_fw_pls[i_loc_Leh]; if (LDA < N_bk_pls[i_loc_Leh]) LDA = N_bk_pls[i_loc_Leh];
  ++LDA; // safety
  
  // The stuff needed by the LAPACK eigenvalue pakage:
  // For 'DSYEVD':
  int INFO;
  int LWORK  = 2+(6+2*LDA)*(LDA);   // For 'DSYEVD' only...
  int LIWORK = 3+5*LDA;          // For 'DSYEVD' only...
  int LWopt=0, LIWopt=0;
  
  int IWORK[LIWORK];
  double MTX[LDA*LDA], W[LDA], WORK[LWORK], MN_new[LDA];
  
  
  
  int nc, ne;
  double x1,x2;
  //
  //  Do hole part:
  //
  for (nc=0; nc<(LDA*LDA); ++nc) MTX[nc] = 0.0;
  
  int N_dim = N_bk_pls[i_loc_Leh];
  for (nc=0; nc<N_dim; ++nc) {
    MTX[(LDA+1)*nc       ] = ek_bk[i_loc_Leh][nc];
    MTX[(LDA+1)*nc +  1  ] = Bk_bk[i_loc_Leh][nc];
    MTX[(LDA+1)*nc + LDA ] = Bk_bk[i_loc_Leh][nc];
  }
  
  
  // LAPACK library (w/ DSYEVD):
  INFO = 0;
  /*408*/
  dsyevd_(&Vectors,&Upper,&N_dim,MTX,&LDA,W,WORK, &LWORK,IWORK,&LIWORK,&INFO);
  if (0 != INFO) cout<< "\nDyson (DSYEV), Wrong value of IFAIL: INFO= "<<INFO<<endl;
  if (0 == INFO) {
    ne = int(WORK[0] + 0.01);
    LWopt = (ne > LWopt) ? ne : LWopt;
    LIWopt = (IWORK[0] > LIWopt) ? IWORK[0] : LIWopt;
    /*410*/} else {
      cout << " 'DEEGV' gave INFO = "   << INFO << endl;
      if (INFO < 0) cout << "\n The " <<-INFO <<"-th aggument in line 408 was illegal:\n";
      cout << "\n\n Program has been aborted since IERR != 0."
      <<   "\nAborted at the line 410, Dyson.f!"
      <<   "\n   --> stop.\n\n";
      exit(1);
    }
  
  x1 = 0.0;
  for (ne=0; ne<N_dim; ++ne) {
    x2 = 0.0;
    for (nc=0; nc<N_dim; ++nc) x2 += MTX[LDA*ne + nc] * Ak_bk[i_loc_Leh][nc];
    MN_new[ne] = x2;
    x1 += x2*x2;
  }
  
  Tot_bk_str[i_loc_Leh] = x1;
  for (ne=0; ne<N_dim; ++ne) {
    ek_bk[i_loc_Leh][ne] = W[ne];
    Ak_bk[i_loc_Leh][ne] = MN_new[ne];
    Bk_bk[i_loc_Leh][ne] = 0.0;
  }
  
  
  //
  //  Do particle part:
  //
  for (nc=0; nc<(LDA*LDA); ++nc) MTX[nc] = 0.0;
  
  N_dim = N_fw_pls[i_loc_Leh];
  for (nc=0; nc<N_dim; ++nc) {
    MTX[(LDA+1)*nc       ] = ek_fw[i_loc_Leh][nc];
    MTX[(LDA+1)*nc +  1  ] = Bk_fw[i_loc_Leh][nc];
    MTX[(LDA+1)*nc + LDA ] = Bk_fw[i_loc_Leh][nc];
  }
  
  // LAPACK library (w/ DSYEVD):
  INFO = 0;
  /*408*/
  dsyevd_(&Vectors,&Upper,&N_dim,MTX,&LDA,W,WORK, &LWORK,IWORK,&LIWORK,&INFO);
  if (0 != INFO) cout<< "\nDyson (DSYEV), Wrong value of IFAIL: INFO= "<<INFO<<endl;
  if (0 == INFO) {
    ne = int(WORK[0] + 0.01);
    LWopt = (ne > LWopt) ? ne : LWopt;
    LIWopt = (IWORK[0] > LIWopt) ? IWORK[0] : LIWopt;
    /*410*/} else {
      cout << " 'DEEGV' gave INFO = "   << INFO << endl;
      if (INFO < 0) cout << "\n The " <<-INFO <<"-th aggument in line 408 was illegal:\n";
      cout << "\n\n Program has been aborted since IERR != 0."
      <<   "\nAborted at the line 410, Dyson.f!"
      <<   "\n   --> stop.\n\n";
      exit(1);
    }
  
  x1 = 0.0;
  for (ne=0; ne<N_dim; ++ne) {
    x2 = 0.0;
    for (nc=0; nc<N_dim; ++nc) x2 += MTX[LDA*ne + nc] * Ak_fw[i_loc_Leh][nc];
    MN_new[ne] = x2;
    x1 += x2*x2;
  }
  Tot_fw_str[i_loc_Leh] = x1;
  for (ne=0; ne<N_dim; ++ne) {
    ek_fw[i_loc_Leh][ne] = W[ne];
    Ak_fw[i_loc_Leh][ne] = MN_new[ne];
    Bk_fw[i_loc_Leh][ne] = 0.0;
  }
  
  
  return 0;}


static int I_p_h_DIVIDE = -100;

double SpctDist::Make_LorSF_vs_E(int iLeh_ksp, double **Emesh, double **ReSF, double **ImSF, int *N_mesh, bool InvIm/*=false*/) {
  
  
  double Gamma_1 = GammaPlot_1;
  double Gamma_2 = GammaPlot_2;
  
  double xe_mn, xeF_h, xeF_p, xe_mx;

  double Ef = Get_EFermi(&xeF_h, &xeF_p, &xe_mn, &xe_mx);

  double xe_rad = 2.0 * (xeF_p - xeF_h );
  xe_mn -= 5.0;
  xe_mx += 5.0;
  
  *N_mesh = int(20.0*xe_rad/Gamma_1) + int(10.0*(xe_mx - xe_mn - 2.0*xe_rad)/Gamma_2) + 5;
  
  if (NULL != *Emesh) delete [] *Emesh;   *Emesh = new double[*N_mesh];
  if (NULL != *ReSF ) delete [] *ReSF;    *ReSF  = new double[*N_mesh];
  if (NULL != *ImSF ) delete [] *ImSF;    *ImSF  = new double[*N_mesh];
  
  double x1,x2,x3;
  
  I_p_h_DIVIDE = -10;
  x1 = xe_mn - 2*Gamma_2;
  for (int i=0; i<(*N_mesh); ++i) {
    (*Emesh)[i] = x1;
    (*ReSF)[i] = 0.0;
    (*ImSF)[i] = 0.0;
    
    if (abs(x1-Ef) > xe_rad)   {x1 += Gamma_2/10.0;}  else {x1 += Gamma_1/10.0;}
    if ((x1 > xeF_h+5*Gamma_1) and (I_p_h_DIVIDE < 0)) {x1 = xeF_p-5*Gamma_1; I_p_h_DIVIDE = i;}
  }
  
  for (int i=0; i<(*N_mesh); ++i) {
    
    x1 = Gamma_1; if (abs(x1-Ef) > xe_rad) x2 = Gamma_2;
    
    for (int ih=0; ih<N_bk_pls[iLeh_ksp]; ++ih) {
      x2 = (*Emesh)[i] - ek_bk[iLeh_ksp][ih];
      x3 = x1*x1 + x2*x2;
      (*ReSF)[i] += x2 * Ak_bk[iLeh_ksp][ih] * Ak_bk[iLeh_ksp][ih] / x3;
      (*ImSF)[i] += x1 * Ak_bk[iLeh_ksp][ih] * Ak_bk[iLeh_ksp][ih] / x3;
    }
    
    if (InvIm) x1 = -x1;
    for (int ip=0; ip<N_fw_pls[iLeh_ksp]; ++ip) {
      x2 = (*Emesh)[i] - ek_fw[iLeh_ksp][ip];
      x3 = x1*x1 + x2*x2;
      (*ReSF)[i] += x2 * Ak_fw[iLeh_ksp][ip] * Ak_fw[iLeh_ksp][ip] / x3;
      (*ImSF)[i] -= x1 * Ak_fw[iLeh_ksp][ip] * Ak_fw[iLeh_ksp][ip] / x3;
    }
    
  }
  
  return Ef;}


int SpctDist::Plot_SpectFnct_3D(void ) {

  double *Emesh = NULL;
  double *ReSF  = NULL;
  double *ImSF  = NULL;
  int    N_mesh;

  double xe_mn, ef_lo, ef_up, xe_mx;

  char filename[100];
  
  double E_F = this->Get_EFermi(&ef_lo, &ef_up, &xe_mn, &xe_mx); // Get the Fermi energy and the energy bounds
  
  sprintf(filename, "SpectFunct3D_parts.dat");
  ofstream file_qp(filename, ios::trunc|ios::out);

  file_qp << "#\n#\n#   Lonertzian-smoothed PARTICLE spectral function (both aprticle and hole) as\n";
  file_qp << "#  function of momentum and energy\n";
  file_qp << "#\n#      Range [Ef+ -:- E_max] = "<< ef_up << "  -:-  " << xe_mx <<" MeV\n";
  file_qp << "#\n#      E_F = "<<E_F<<" MeV\n";
  file_qp << "#\n#\n#    k              E               A(k,E)\n";


  sprintf(filename, "SpectFunct3D_holes.dat");
  ofstream file_qh(filename, ios::trunc|ios::out);

  file_qh << "#\n#\n#   Lonertzian-smoothed HOLE spectral function (both aprticle and hole) as\n";
  file_qh << "#  function of momentum and energy\n";
  file_qh << "#\n#      Range [Emin -:- Ef-] = "<< xe_mn << "  -:-  " << ef_lo <<" MeV\n";
  file_qh << "#\n#      E_F = "<<E_F<<" MeV\n";
  file_qh << "#\n#\n#    k              E               A(k,E)\n";

  
  sprintf(filename, "SpectFunct3D_fence_parts.dat");
  ofstream fnce_qp(filename, ios::trunc|ios::out);
  
  fnce_qp << "#\n#\n#   Lonertzian-smoothed PARTICLE spectral function (both aprticle and hole) as\n";
  fnce_qp << "#  function of momentum and energy\n";
  fnce_qp << "#\n#      Range [Ef+ -:- E_max] = "<< ef_up << "  -:-  " << xe_mx <<" MeV\n";
  fnce_qp << "#\n#      E_F = "<<E_F<<" MeV\n";
  fnce_qp << "#\n#\n#    k              E               A(k,E)\n";
  
  
  sprintf(filename, "SpectFunct3D_fence_holes.dat");
  ofstream fnce_qh(filename, ios::trunc|ios::out);
  
  fnce_qh << "#\n#\n#   Lonertzian-smoothed HOLE spectral function (both aprticle and hole) as\n";
  fnce_qh << "#  function of momentum and energy\n";
  fnce_qh << "#\n#      Range [Emin -:- Ef-] = "<< xe_mn << "  -:-  " << ef_lo <<" MeV\n";
  fnce_qh << "#\n#      E_F = "<<E_F<<" MeV\n";
  fnce_qh << "#\n#\n#    k              E               A(k,E)\n";

  for (int i_loc_Leh=0; i_loc_Leh<Ndim_Leh; ++i_loc_Leh) {
    
    //
    //' true' means that the imaginary part in the particle side is taken with a
    // minus sign (as needed for the spectral function).
    double E_F = Make_LorSF_vs_E(i_loc_Leh, &Emesh, &ReSF, &ImSF, &N_mesh, true);

    file_qp << endl;
    file_qh << endl;
    fnce_qp << endl;
    fnce_qh << endl;

    for (int i=0; i<=I_p_h_DIVIDE; ++i) {
      file_qh << SpBasLoc->k[this->i_sp[i_loc_Leh]] <<"    "<< Emesh[i] << "       "  << ImSF[i]/PI <<endl;
      fnce_qh << SpBasLoc->k[this->i_sp[i_loc_Leh]] <<"    "<< Emesh[i] << "       "  << 0.0 <<endl;
      fnce_qh << SpBasLoc->k[this->i_sp[i_loc_Leh]] <<"    "<< Emesh[i] << "       "  << ImSF[i]/PI <<endl;
      fnce_qh << SpBasLoc->k[this->i_sp[i_loc_Leh]] <<"    "<< Emesh[i] << "       "  << 0.0 <<endl;
    }

    for (int i=I_p_h_DIVIDE+1; i<N_mesh; ++i) {
      file_qp << SpBasLoc->k[this->i_sp[i_loc_Leh]] <<"    "<< Emesh[i] << "       "  << ImSF[i]/PI <<endl;
      fnce_qp << SpBasLoc->k[this->i_sp[i_loc_Leh]] <<"    "<< Emesh[i] << "       "  << 0.0 <<endl;
      fnce_qp << SpBasLoc->k[this->i_sp[i_loc_Leh]] <<"    "<< Emesh[i] << "       "  << ImSF[i]/PI <<endl;
      fnce_qp << SpBasLoc->k[this->i_sp[i_loc_Leh]] <<"    "<< Emesh[i] << "       "  << 0.0 <<endl;
    }

    file_qp << endl;
    file_qh << endl;
    fnce_qp << endl;
    fnce_qh << endl;
  }

  file_qp.close();
  file_qh.close();
  fnce_qp.close();
  fnce_qh.close();

  if (NULL != Emesh) delete [] Emesh;
  if (NULL != ReSF ) delete [] ReSF;
  if (NULL != ImSF ) delete [] ImSF;

  
  return 0;}



int SpctDist::Plot_SelfEn_vs_E(int i_ksp) {
  
  double *Emesh = NULL;
  double *ReSF  = NULL;
  double *ImSF  = NULL;
  int    N_mesh;
  
//  int i_loc_Leh = Get_Leh_index(iLeh_ksp);
  int i_loc_Leh = -100;
  for (int i=0; i<Ndim_Leh; ++i)
    if (i_ksp == this->i_sp[i]) {i_loc_Leh=i; break;}
  
  if (i_loc_Leh != Get_Leh_index(i_ksp)) {
    cerr << "AAAAAAAA1\n\n;"; exit(1);
  }

  this->DiagTridiagSE(i_loc_Leh);
  
  double E_F = Make_LorSF_vs_E(i_loc_Leh, &Emesh, &ReSF, &ImSF, &N_mesh, false);
  
  char filename[100];
  
  //cout << "plotting Ileh="<<i_ksp<<"     Nmesh"<<N_mesh<<"     nsq"<<nsq[i_loc_Leh]<<" \n"<<flush;
  
  
  sprintf(filename, "SelfEn_k%.3lf_fm-1__grp%i__poles.dat" , SpBasLoc->k[i_ksp] , this->i_grp[i_loc_Leh] );
  
  ofstream file(filename, ios::trunc|ios::out);
  
  file << "#\n#\n#   poles of the irreducible self-energy (both particle and hole)\n";
  file << "#\n#   k = "<<SpBasLoc->k[i_ksp]<<" fm^-1      and E_F = "<<E_F<<" MeV\n";
  file << "#\n#\n#    E              |M_kE|^2\n\n";
  
  for (int ih=0; ih<N_bk_pls[i_loc_Leh]; ++ih) file << ek_bk[i_loc_Leh][ih] << "       " << Ak_bk[i_loc_Leh][ih] << endl;
  file << endl << endl;
  for (int ip=0; ip<N_fw_pls[i_loc_Leh]; ++ip) file << ek_fw[i_loc_Leh][ip] << "       " << Ak_fw[i_loc_Leh][ip] << endl;
  
  file.close();
  
  //
  //
  //
  
  sprintf(filename, "SelfEn_k%.3lf_fm-1__grp%i__Lorntz.dat" , SpBasLoc->k[i_ksp] , this->i_grp[i_loc_Leh] );
  
  
  file.open(filename, ios::trunc|ios::out);
  
  //cout << "plotting Ileh ="<<i_ksp<<"       Nmesh = "<<N_mesh<<" \n"<<flush;
  
  
  file << "#\n#\n#   Lonertzian-smoothed self-energy (both particle and hole)\n";
  file << "#\n#   k = "<<SpBasLoc->k[i_ksp]<<" fm^-1 (group = "<< this->i_grp[i_loc_Leh]<<")     and E_F = "<<E_F<<" MeV\n";
  file << "#\n#   Sigma^{infty} = "<<Sig_inf[i_loc_Leh]<<" MeV\n";
  file << "#\n#\n#    E               Re \\Sig(k,E)       Im \\Sig(k,E)             k \n\n";
  
  
  for (int i=0; i<N_mesh; ++i) {
    if (i-1 == I_p_h_DIVIDE) file  << Emesh[i] << "       " << 0.0 << "       " << 0.0 << endl;
    
    file << Emesh[i] << "      " << ReSF[i] + Sig_inf[i_loc_Leh] << "       " << ImSF[i] << "            "<< SpBasLoc->k[i_ksp] <<endl;
    
    if (i == I_p_h_DIVIDE) file  << Emesh[i] << "       " << 0.0 << "       " << 0.0 << "\n\n\n";

  }
  
  
  file.close();
  
  
  if (NULL != Emesh) delete [] Emesh;
  if (NULL != ReSF ) delete [] ReSF;
  if (NULL != ImSF ) delete [] ImSF;
  
  
  return 0;}



int SpctDist::Plot_SpectFnct_vs_E(int i_ksp) {
  
  double *Emesh = NULL;
  double *ReSF  = NULL;
  double *ImSF  = NULL;
  int    N_mesh;
  
  //  int i_loc_Leh = Get_Leh_index(iLeh_ksp);
  int i_loc_Leh = -100;
  for (int i=0; i<Ndim_Leh; ++i)
    if (i_ksp == this->i_sp[i]) {i_loc_Leh=i; break;}
  
  if (i_loc_Leh != Get_Leh_index(i_ksp)) {
    cerr << "AAAAAAAA2\n\n;"; exit(1);
  }
  
  
  //
  //' true' means that the imaginary part in the particle side is takend with a
  // minus sign (as needed for the spectral function).
  double E_F = Make_LorSF_vs_E(i_loc_Leh, &Emesh, &ReSF, &ImSF, &N_mesh, true);
  
  char filename[100];
  
  //cout << "plotting Ileh="<<i_ksp<<"     Nmesh"<<N_mesh<<"     nsq"<<nsq[i_loc_Leh]<<" \n"<<flush;
  
  
  sprintf(filename, "SpectFunct_k%.3lf_fm-1__grp%i__poles.dat" , SpBasLoc->k[i_ksp] , this->i_grp[i_loc_Leh] );
  
  ofstream file(filename, ios::trunc|ios::out);
  
  file << "#\n#\n#   poles of the spectral functions (both aprticle and hole)\n";
  file << "#\n#   k = "<<SpBasLoc->k[i_ksp]<<" fm^-1 (group = "<< this->i_grp[i_loc_Leh]<<")     and E_F = "<<E_F<<" MeV\n";
  file << "#\n#\n#    E               A(k,E)\n\n";
  
  for (int ih=0; ih<N_bk_pls[i_loc_Leh]; ++ih) file << ek_bk[i_loc_Leh][ih] << "       " << Ak_bk[i_loc_Leh][ih] << endl;
  file << endl << endl;
  for (int ip=0; ip<N_fw_pls[i_loc_Leh]; ++ip) file << ek_fw[i_loc_Leh][ip] << "       " << Ak_fw[i_loc_Leh][ip] << endl;
  
  file.close();
  
  //
  //
  //
  
  sprintf(filename, "SpectFunct_k%.3lf_fm-1__grp%i__Lorntz.dat" , SpBasLoc->k[i_ksp] , this->i_grp[i_loc_Leh] );
  
  
  file.open(filename, ios::trunc|ios::out);
  
  //cout << "plotting Ileh ="<<i_ksp<<"       Nmesh = "<<N_mesh<<" \n"<<flush;
  
  
  file << "#\n#\n#   Lonertzian-smoothed spectral function (both aprticle and hole)\n";
  file << "#\n#   k = "<<SpBasLoc->k[i_ksp]<<" fm^-1      and E_F = "<<E_F<<" MeV\n";
  file << "#\n#\n#    E               A(k,E)           k \n\n";
  
  for (int i=0; i<N_mesh; ++i) {
    if (i-1 == I_p_h_DIVIDE) file  << Emesh[i] << "       " << 0.0 << "       " << 0.0 << endl;
    
    file << Emesh[i] << "       " << ImSF[i]/PI << "       "<< SpBasLoc->k[i_ksp] <<endl;
    
    if (i == I_p_h_DIVIDE) file  << Emesh[i] << "       " << 0.0 << "       " << 0.0 << "\n\n\n";
  }
  
  file.close();
  
  
  if (NULL != Emesh) delete [] Emesh;
  if (NULL != ReSF ) delete [] ReSF;
  if (NULL != ImSF ) delete [] ImSF;
  
  
  return 0;}




double SpctDist::Koltun(void ) {
  
  double xkin, xA, xKlt, xSF;
  
  
  this->Atot = 0.0;
  this->EKoltun = 0.0;
  this->Tkin = 0.0;
  for (int iLeh_ksp=0; iLeh_ksp<Ndim_Leh; ++iLeh_ksp) {
    xkin = this->t_kin[iLeh_ksp];
    xA = 0.0;
    xKlt = 0.0;
    for (int ifgr=0; ifgr<N_bk_pls[iLeh_ksp]; ++ifgr) {

      xSF =  pow(Ak_bk[iLeh_ksp][ifgr], 2);
    
      xKlt += (ek_bk[iLeh_ksp][ifgr] ) * xSF;
    
      xA += xSF;
    }
    
    xKlt += xkin * xA;
    xKlt *= i_grp_mult[iLeh_ksp] / 2.0;
    xA   *= i_grp_mult[iLeh_ksp];
    
    this->Atot += xA;
    this->EKoltun += xKlt;
    this->Tkin += xkin * xA;
  }

  return this->EKoltun;}



double SpctDist::Set_EFermi(double EF_in ) {

  int i_max, ifgr;
  double x1;
  
  this->EFermi = EF_in;

  this->Atot = 0.0;
  for (int iLeh_ksp=0; iLeh_ksp<Ndim_Leh; ++iLeh_ksp) {
    x1 = 0.0;
    i_max = N_fw_pls[iLeh_ksp]+N_bk_pls[iLeh_ksp];
    //cout << " imax="<< i_max <<"      Zb, Ztot:" << Tot_bk_str[iLeh_ksp] <<"    "
    //<<Tot_fw_str[iLeh_ksp] +Tot_bk_str[iLeh_ksp]  << flush;
//    bool found = false;
    for (ifgr=0; ifgr<i_max; ++ifgr)
      if (ek_bk[iLeh_ksp][ifgr] <= this->EFermi) {x1+=pow(Ak_bk[iLeh_ksp][ifgr], 2);} else {break;}
    
    N_bk_pls[iLeh_ksp] = ifgr;
    N_fw_pls[iLeh_ksp] = i_max - ifgr;
    Tot_fw_str[iLeh_ksp] += Tot_bk_str[iLeh_ksp] - x1;
    Tot_bk_str[iLeh_ksp] = x1;
    //
    ek_fw[iLeh_ksp] = ek_bk[iLeh_ksp] + N_bk_pls[iLeh_ksp];
    Ak_fw[iLeh_ksp] = Ak_bk[iLeh_ksp] + N_bk_pls[iLeh_ksp];
    Bk_fw[iLeh_ksp] = Bk_bk[iLeh_ksp] + N_bk_pls[iLeh_ksp];
    
    this->Atot += x1 * i_grp_mult[iLeh_ksp];
    //cout << " -->  Zb, Ztot:" << Tot_bk_str[iLeh_ksp] <<"    "
    //<<Tot_fw_str[iLeh_ksp] +Tot_bk_str[iLeh_ksp] <<"     Atot="<<Atot <<endl << flush;
  }
  
  return this->Atot;}



double SpctDist::Seek_Fermi_level(double N_trgt ) {
  
  double x1, E_lo, E_up, Ef_st, N_lo, N_up;
  
  int i_pls, i_max, n1, n2;

  Ef_st = Get_EFermi();
  x1 = 0.0;
  
  N_lo = N_up = N_trgt + 1.0;
  while ((N_lo >= N_trgt) || (N_up <= N_trgt)) {
    x1 += 1.0;
    E_lo = Ef_st - x1;
    E_up = Ef_st + x1;
    N_lo = this->Set_EFermi(E_lo);
    N_up = this->Set_EFermi(E_up);
  }
  x1 += 1.0;

  n1 = 0;
  for (int iLeh_ksp=0; iLeh_ksp<Ndim_Leh; ++iLeh_ksp) {
    i_max = N_fw_pls[iLeh_ksp]+N_bk_pls[iLeh_ksp];
    for (int ifgr=0; ifgr<i_max; ++ifgr) {
      x1 = ek_bk[iLeh_ksp][ifgr];
      if ( (E_lo < x1) && (x1 < E_up) ) ++n1;
  } }
  i_pls = n1;

  double E_loc[i_pls+1], Z_loc[i_pls+1];

  n1 = 1;
  for (int iLeh_ksp=0; iLeh_ksp<Ndim_Leh; ++iLeh_ksp) {
    i_max = N_fw_pls[iLeh_ksp]+N_bk_pls[iLeh_ksp];
    for (int ifgr=0; ifgr<i_max; ++ifgr) {
      x1 = ek_bk[iLeh_ksp][ifgr];
      if ( (E_lo < x1) && (x1 < E_up) ) {
        E_loc[n1] = x1;
        //Z_loc[n1] = Ak_bk[iLeh_ksp][i_max] * Ak_bk[iLeh_ksp][i_max];
        ++n1;}
    } }

  Sort_u_double2dim(i_pls, 1, 0, E_loc+1);

  x1 = E_lo;
  for (n1=0; n1<i_pls; ++n1) {E_loc[n1] = (x1+E_loc[n1+1])/2.0; x1 = E_loc[n1+1];}
  E_loc[i_pls] = (x1+E_up)/2.0;


  x1 = abs(N_trgt) + 1.0;
  n2 = -100;
  for (n1=0; n1<=i_pls; ++n1) {
    Z_loc[n1] = this->Set_EFermi(E_loc[n1]);
    if ( abs(Z_loc[n1]-N_trgt) < x1 ) {x1 = abs(Z_loc[n1]-N_trgt); n2=n1;}
 }
 
 
  this->Set_EFermi(E_loc[n2]);
 
  return this->Atot;}



  


