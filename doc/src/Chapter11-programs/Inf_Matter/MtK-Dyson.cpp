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
//  MtK-Dyson  --  Diagonalisaton of Dyson matrix; here is where the action is ;-).
//


#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

#include "MtK-Bases.h"
#include "MtK-SpProppagator.h"
#include "MtK-Potential_minnesota.h"
#include "MtK-Global_data.h"


// LAPACK routines
//----------------
extern "C" void dsyevd_(char*, char*, int*, double*, int*,
                        double*, double*, int*, int*, int*,  int* );

static char Vectors = 'V';
static char Upper   = 'U';


//
//  Gerenral driver to solve the Dyson equation recursively for each group of single particle states.
//
int Solve_Dyson(ADC3BasisK *ADC3Bas_in, double *N_tot, double *E_Kolt, double *E_kin) {

  double xN,xE,xK, x1;
  bool GrpsFlag = true; // set == false to diagonalize ALL s.p. channels never really needed)

  ADC3Bas_in->SpBasLoc->E_Fermi = (ADC3Bas_in->SpBasLoc->e_sp[ADC3Bas_in->SpBasLoc->N_holes-1]
                                + ADC3Bas_in->SpBasLoc->e_sp[ADC3Bas_in->SpBasLoc->N_holes])/2.0;

  cout << "\n\n Assuming E_F = "<<ADC3Bas_in->SpBasLoc->E_Fermi<< " MeV for the Fermi energy.\n";
  
  *N_tot  = 0.0;
  *E_Kolt = 0.0;
  *E_kin  = 0.0;
  for (int isp=0; isp<ADC3Bas_in->SpBasLoc->SpNmax; ++isp) {

    x1 = 0.0;
    if (GrpsFlag) {
      if (isp != ADC3Bas_in->SpBasLoc->gr_rep[ADC3Bas_in->SpBasLoc->group[isp]]) continue;
      //
      x1 = ADC3Bas_in->SpBasLoc->gr_mlt[ADC3Bas_in->SpBasLoc->group[isp]]; // # of equivalent channels
    }

    Solve_Dyson_sp(ADC3Bas_in, isp, &xN, &xE, &xK, g_sp_out);

    (*N_tot ) += xN *x1;
    (*E_Kolt) += xE *x1;
    (*E_kin ) += xK *x1;
    
  }
  
  
  //
  //  At this point *N_tot, *E_Kolt, *E_kin ocntain solutions from the Dyson equation, but based on
  // the above assumption for E_F. If that is not good, the total numebr fo particle will deviate from
  // the target valuse (ADC3Bas_in->SpBasLoc->N_holes) and all sum rules from the spectral function
  // will be incorrect.
  //
  //  Thus we first reset the Fermi energy to get the number of particle as closest as possible to
  // our intended one and the realculate all sum rules. the function SpctDist::Seek_Fermi_level(double A)
  // does this for us:
  //
  
  (*N_tot ) = g_sp_out->Seek_Fermi_level(double(ADC3Bas_in->SpBasLoc->N_holes));
  
  cout << "\n New Fermi level reset to E_F = "  << g_sp_out->EFermi  << " ,   with   A_tot = "  << g_sp_out->Atot << flush;
  
  g_sp_out->Koltun();
  (*E_Kolt) = g_sp_out->EKoltun;
  (*E_kin ) = g_sp_out->Tkin;
  (*N_tot ) = g_sp_out->Atot;
  
  cout << " ==>  E_Kolt = "  << g_sp_out->EKoltun  << " MeV ,   E_kin = "  << g_sp_out->Tkin << " MeV,   A_tot = "  << g_sp_out->Atot << endl << endl  << flush;
  
  //
  //  Since Dyson is not self consistent, it will break sligly the number of particles. Here we fine
  // tune our result by recalculation the density and rescaling the energies w.r.t. the calculated
  // value of A_tot = "  << g_sp_out->Atot;
  
  double rho_loc = (*N_tot)/ pow(ADC3Bas_in->SpBasLoc->Lbox , 3.0);
  double kF_loc =  6.0 * PI * PI * rho_loc /double(2*(ADC3Bas_in->SpBasLoc->ch_max-ADC3Bas_in->SpBasLoc->ch_min+1));
  kF_loc = pow( kF_loc , (1.0/3.0) );
  
  cout << "\n\n The final density, Fermi momentum and energy per particle are:";
  cout << "\n Solution from the Dyson diagonalization.:";
  cout << "\n rho  = "  << rho_loc << " fm^-3 ,    f_F  = "  << kF_loc << " fm^-1";
  cout << "\n A_tot  = "  << *N_tot;
  cout << "\n E_tot  = "  << *E_Kolt << " ,    E_tot/A_tot = "  << (*E_Kolt)/(*N_tot) << "  [MeV]";
  cout << "\n K_tot  = "  << *E_kin  << " ,    K_tot/A_tot = "  << (*E_kin)/(*N_tot)  << "  [MeV]\n\n";

  return 0;}



//
//  Solve the Dyson eq for the first iteration, here the static self energy is simply the HF potential.
//
int Solve_Dyson_sp(ADC3BasisK *ADC3Bas_in, int isp, double *N_out, double *E_out, double *K_out, SpctDist *gsp_out/*=NULL*/) {
  
  
  ADC3Bas_in->Build_2p1h_2h1p_basis(isp);
  
  SpBasisK *WK_bas = ADC3Bas_in->SpBasLoc;
  
  int Ndim_sig_fw = ADC3Bas_in->Nbas_2p1h;
  int Ndim_sig_bk = ADC3Bas_in->Nbas_2h1p;
  
  bool DoLanc_fw = false;
  //if ( (NLanczos > 0) && (NLanczos < (9*Ndim_sig_fw/10)) ) {
  if  (NLanczos > 0)  {
    Ndim_sig_fw = NLanczos;
    DoLanc_fw = true;
    int i1 = int(0.98 * ADC3Bas_in->Nbas_2p1h); if (Ndim_sig_fw>i1) Ndim_sig_fw=i1;
  }
  
  bool DoLanc_bk = false;
  //if ( (NLanczos > 0) && (NLanczos < (9*Ndim_sig_bk/10)) ) {
  if   (NLanczos > 0)  {
    Ndim_sig_bk = NLanczos;
    DoLanc_bk = true;
    int i1 = int(0.98 * ADC3Bas_in->Nbas_2h1p); if (Ndim_sig_bk>i1) Ndim_sig_bk=i1;
  }
  
  cout << isp << "  Nfw=" << Ndim_sig_fw << "  Nbk=" << Ndim_sig_bk << flush;
  
  int ndim = 1 + Ndim_sig_fw + Ndim_sig_bk;
  
  
  // The stuff needed by the LAPACK eigenvalue pakage:
  // For 'DSYEV':
  int INFO, ne;
  int LDA = ndim+1;
  //int LWORK  = 3*LDA+2;          // For 'DSYEV'  only...
  int LWORK  = 2+(6+2*ndim)*(ndim);   // For 'DSYEVD' only...
  int LIWORK = 3+5*ndim;          // For 'DSYEVD' only...
  
  //cout << " ,  LDA    = " << LDA   ;
  //cout << " ,  LWORK  = " << LWORK ;
  //cout << " ,  LIWORK = " << LIWORK;
  
  double *MTX, *W, *WORK;
  int    *IWORK;
  W     = new double[LDA];
  WORK  = new double[LWORK];
  IWORK = new    int[LIWORK];
  
  int LWopt=0, LIWopt=0;
  
  MTX = new double[LDA*LDA];
  for (int i=0; i<(LDA*LDA); ++i) MTX[i]=0.0;
  
  
  
  //
  //  the structure of the Dyson Matrix is as follow:
  //
  //         /     |                 |                \
  //         | HF  |   M2p1h^dag     |    M2h1p^dag   |
  //         |     |                 |                |
  //         |-----+-----------------+----------------|
  //         |     |                 |                |
  //         | M   |                 |                |
  //         |2p1h |     diag{       |                |
  //         |     |     e_2p1h      |                |
  //   MTX = {     |       }         |      0         |
  //         |     |                 |                |
  //         |     |                 |                |
  //         |-----+-----------------+----------------|
  //         |     |                 |                |
  //         |     |                 |                |
  //         |  M  |                 |     diag{      |
  //         |2h1p |        0        |    e_2h1p      |
  //         |     |                 |       }        |
  //         |     |                 |                |
  //         \     |                 |                /
  //
  //
  //
  
  cout << " , ndim=" << ndim << flush;
  
  
  MTX[0] = HF_Minn_MtxEls(WK_bas, isp, isp, WK_bas->N_holes);
  
  
  double *dptr_1 = MTX        + 1;
  double *dptr_2 = MTX + LDA    ;
  double *dptr_3 = MTX + LDA + 1;
  
  cout << " -- "<< dptr_1-MTX << "...  " << flush;
  
  if (DoLanc_fw) {Build_DysonADC_LancMatrix( ADC3Bas_in, +1, Ndim_sig_fw, dptr_1, dptr_3, LDA);}
  else {Build_DysonADC_Matrix_2p1h(ADC3Bas_in,                  dptr_1, dptr_3, LDA);}
  for (int i1=0; i1 < Ndim_sig_fw; ++i1) {
    dptr_2[i1*LDA]=dptr_1[i1];
  }
  //
  dptr_1 += Ndim_sig_fw;
  dptr_2 += LDA*Ndim_sig_fw;
  dptr_3 += (LDA+1)*Ndim_sig_fw; // == MTX + (LDA+1)*(1+Ndim_sig_fw)
  //
  cout << dptr_1-MTX << "...  " << flush;
  
  if (DoLanc_bk) {Build_DysonADC_LancMatrix( ADC3Bas_in, -1, Ndim_sig_bk, dptr_1, dptr_3, LDA);}
  else {Build_DysonADC_Matrix_2h1p(ADC3Bas_in,                  dptr_1, dptr_3, LDA);}
  for (int i1=0; i1 < Ndim_sig_bk; ++i1) {
    dptr_2[i1*LDA]=dptr_1[i1];
  }
  //
  dptr_1 += Ndim_sig_bk;
  dptr_2 += LDA*Ndim_sig_bk;
  //
  cout << dptr_1-MTX;
  
  
  cout << " --   group=" << WK_bas->group[isp];
  cout << " " << WK_bas->gr_rep[WK_bas->group[isp]]<< flush;
  cout << " " << WK_bas->nsq[WK_bas->gr_rep[WK_bas->group[isp]]]<< flush;
  cout << " " << WK_bas->gr_mlt[WK_bas->group[isp]]<< flush;
  cout << " (e_sp=" << WK_bas->e_sp[isp]<<") "<< flush;
  
  {
    // Store the self-energy:
    double  M_fw[LDA], CA_sol_fw[LDA], CB_sol_fw[LDA];
    double  N_bk[LDA], DA_sol_bk[LDA], DB_sol_bk[LDA];
    
    for (int i_sig_pls=0; i_sig_pls<Ndim_sig_fw; ++i_sig_pls) {
      M_fw[i_sig_pls] = MTX[ 1  + i_sig_pls];
      CA_sol_fw[i_sig_pls] = MTX[ (LDA+1) * (1+i_sig_pls) ];
      if (i_sig_pls < Ndim_sig_fw-1) CB_sol_fw[i_sig_pls] = MTX[ (LDA+1) * (1+i_sig_pls) +1 ];
      else CB_sol_fw[i_sig_pls] =0.0;
    }
    
    for (int i_sig_pls=0; i_sig_pls<Ndim_sig_bk; ++i_sig_pls) {
      N_bk[i_sig_pls] = MTX[ 1  + Ndim_sig_fw + i_sig_pls];
      DA_sol_bk[i_sig_pls] = MTX[ (LDA+1) * (1+Ndim_sig_fw+i_sig_pls) ];
      if (i_sig_pls < Ndim_sig_fw-1) DB_sol_bk[i_sig_pls] = MTX[ (LDA+1) * (1+Ndim_sig_fw+i_sig_pls) +1 ];
      else DB_sol_bk[i_sig_pls] =0.0;
    }
    
    // Save the irreducilble self energy
    Sigma_irred->add_k_channel(isp,  Ndim_sig_fw,  M_fw,  CA_sol_fw,  Ndim_sig_bk,  N_bk,  DA_sol_bk,
                               MTX[0], CB_sol_fw, DB_sol_bk);
  }
  
  
  
  
  // LAPACK library (w/ DSYEVD):
  INFO = 0;
  /*408*/
  dsyevd_(&Vectors,&Upper,&ndim,MTX,&LDA,W,WORK, &LWORK,IWORK,&LIWORK,&INFO);
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
  
  bool qparticle;
  int ig;
  double x1;
  
  int     n_sol_fw = 0, n_sol_bk = 0;
  double  A_sol_fw[LDA], e_sol_fw[LDA];
  double  A_sol_bk[LDA], e_sol_bk[LDA];
  
  (*N_out) = 0.0;
  (*E_out) = 0.0;
  for(ne=0; ne<ndim; ++ne) {
    dptr_1=MTX+ne*LDA;
    x1=0.0;
    for (ig=0; ig<ndim; ++ig) x1 += dptr_1[ig]*dptr_1[ig];
    x1 = sqrt(x1);
    for (ig=0; ig<ndim; ++ig) dptr_1[ig] /= x1;
    x1 = dptr_1[0]*dptr_1[0];
    
    qparticle = false; if (W[ne] > WK_bas->E_Fermi) qparticle = true;
    
    if (qparticle) {
      //
      //  Save qp spectral function
      //
      A_sol_fw[n_sol_fw] = sqrt(x1); // the overall phase of spect amplitudes does not matter
      e_sol_fw[n_sol_fw] = W[ne];
      ++n_sol_fw;
      
    } else {
      //
      //  Save qh spectral function
      //
      A_sol_bk[n_sol_bk] = sqrt(x1); // the overall phase of spect amplitudes does not matter
      e_sol_bk[n_sol_bk] = W[ne];
      ++n_sol_bk;
      
      //
      //  Calculate contribution for the Koltun SR
      //
      (*N_out) += x1;
      (*E_out) += x1 * W[ne];
    }
    
  } // en loop `ne'
  (*K_out) = (*N_out) * WK_bas->e_kin[isp];
  (*E_out) = ( (*K_out) + (*E_out) ) / 2.0;
  
  cout << "  N="<<(*N_out)<< "  K="<<(*K_out)<< "  V="<<(*E_out)-(*K_out)<< "  E="<<(*E_out);
  
  cout << endl;
  
  //
  //  Finally, store the spectral function in the output propagator:
  gsp_out->add_k_channel(isp,  n_sol_fw,  A_sol_fw,  e_sol_fw,  n_sol_bk,  A_sol_bk,  e_sol_bk);
  
  delete [] W;
  delete [] WORK;
  delete [] IWORK;
  
  delete [] MTX;
  
  return 0;}


//
//  Builds the interaction (E^> + C) and coupling (M) matrices for the forward ADC(n) calculation (2p1h ISCs)
//
int Build_DysonADC_Matrix_2p1h(ADC3BasisK *ADC3Bas_in, double *Mst, double *Cst, int LDC ) {

  SpBasisK *WK_bas = ADC3Bas_in->SpBasLoc;
  
  int isp = ADC3Bas_in->iSpLoc;
  
  int im,iv, ig;
  double x1, x2;
  int  kin_x, kin_y, kin_z, sp_in, ch_in;
  int   kt_x,  kt_y,  kt_z, sp_t,  ch_t;
  int   h5_x,  h5_y,  h5_z, sp_5,  ch_5,  h5,  h4;
  int   n7_x,  n7_y,  n7_z, sp_7,  ch_7,  n7;//,  n6;
  
  
  int    *iptr_b = ADC3Bas_in->Bas_2p1h;
  double *dptr_1 = Mst;
  
  for (int i1=0; i1 < ADC3Bas_in->Nbas_2p1h; ++i1) {
    im = iptr_b[0];
    iv = iptr_b[1];
    ig = iptr_b[2];
    //   cout << " \n    bbb  " << i1 << flush;
    
    x1 = V_Minnesota(WK_bas,isp,ig,im,iv);
    *dptr_1 = x1;
    
    //
    //  Start with ADC(3) coupling vertices
    //
    kin_x = WK_bas->nx[isp];
    kin_y = WK_bas->ny[isp];
    kin_z = WK_bas->nz[isp];
    sp_in = WK_bas->spin[isp];
    ch_in = WK_bas->chrg[isp];
    
    kt_x = WK_bas->nx[ig]   + kin_x;
    kt_y = WK_bas->ny[ig]   + kin_y;
    kt_z = WK_bas->nz[ig]   + kin_z;
    sp_t = WK_bas->spin[ig] + sp_in;
    ch_t = WK_bas->chrg[ig] + ch_in;
    //
    x1 = 0.0;
    x2 = WK_bas->e_sp[im] + WK_bas->e_sp[iv];
    for (h4=0; h4<WK_bas->N_holes; ++h4) {
      h5_x = kt_x - WK_bas->nx[h4];
      h5_y = kt_y - WK_bas->ny[h4];
      h5_z = kt_z - WK_bas->nz[h4];
      sp_5 = sp_t - WK_bas->spin[h4];
      ch_5 = ch_t - WK_bas->chrg[h4];
      
      for (h5=h4+1; h5<WK_bas->N_holes; ++h5) {
        
        if ( (h5_x == WK_bas->nx[h5]  ) && (h5_y == WK_bas->ny[h5]  ) && (h5_z == WK_bas->nz[h5]) &&
            (sp_5 == WK_bas->spin[h5]) && (ch_5 == WK_bas->chrg[h5]) ) {
          
          x1 += V_Minnesota(WK_bas,im,iv,h4,h5) * V_Minnesota(WK_bas,h4,h5,isp,ig)
          / (WK_bas->e_sp[h4] + WK_bas->e_sp[h5] - x2);
        }
        
      }
    }
    //
    *dptr_1 += x1;
    
    kt_x = WK_bas->nx[im]   - kin_x;
    kt_y = WK_bas->ny[im]   - kin_y;
    kt_z = WK_bas->nz[im]   - kin_z;
    sp_t = WK_bas->spin[im] - sp_in;
    ch_t = WK_bas->chrg[im] - ch_in;
    //
    x1 = 0.0;
    x2 = WK_bas->e_sp[ig] - WK_bas->e_sp[iv];
    for (h5=0; h5<WK_bas->N_holes; ++h5) {
      n7_x = kt_x + WK_bas->nx[h5];
      n7_y = kt_y + WK_bas->ny[h5];
      n7_z = kt_z + WK_bas->nz[h5];
      sp_7 = sp_t + WK_bas->spin[h5];
      ch_7 = ch_t + WK_bas->chrg[h5];
      
      for (n7=WK_bas->N_holes; n7<WK_bas->SpNmax; ++n7) {
        
        if ( (n7_x == WK_bas->nx[n7]  ) && (n7_y == WK_bas->ny[n7]  ) && (n7_z == WK_bas->nz[n7]) &&
            (sp_7 == WK_bas->spin[n7]) && (ch_7 == WK_bas->chrg[n7]) ) {
          
          x1 += V_Minnesota(WK_bas,iv,n7,ig,h5) * V_Minnesota(WK_bas,im,h5,isp,n7)
          / (WK_bas->e_sp[h5] - WK_bas->e_sp[n7] + x2);
        }
        
      }
    }
    //
    *dptr_1 += x1;
    
    kt_x = WK_bas->nx[iv]   - kin_x;
    kt_y = WK_bas->ny[iv]   - kin_y;
    kt_z = WK_bas->nz[iv]   - kin_z;
    sp_t = WK_bas->spin[iv] - sp_in;
    ch_t = WK_bas->chrg[iv] - ch_in;
    //
    x1 = 0.0;
    x2 = WK_bas->e_sp[ig] - WK_bas->e_sp[im];
    for (h5=0; h5<WK_bas->N_holes; ++h5) {
      n7_x = kt_x + WK_bas->nx[h5];
      n7_y = kt_y + WK_bas->ny[h5];
      n7_z = kt_z + WK_bas->nz[h5];
      sp_7 = sp_t + WK_bas->spin[h5];
      ch_7 = ch_t + WK_bas->chrg[h5];
      
      for (n7=WK_bas->N_holes; n7<WK_bas->SpNmax; ++n7) {
        
        if ( (n7_x == WK_bas->nx[n7]  ) && (n7_y == WK_bas->ny[n7]  ) && (n7_z == WK_bas->nz[n7]) &&
            (sp_7 == WK_bas->spin[n7]) && (ch_7 == WK_bas->chrg[n7]) ) {
          
          x1 += V_Minnesota(WK_bas,im,n7,ig,h5) * V_Minnesota(WK_bas,iv,h5,isp,n7)
          / (WK_bas->e_sp[h5] - WK_bas->e_sp[n7] + x2);
        }
        
      }
    }
    //
    *dptr_1 -= x1;
    //
    //  end of ADC(3) coupling vertices
    //
    
    iptr_b += 3;
    dptr_1 += 1;
  }


  
  int ia, ib, ir;
  int *iptr_b2;
  //
  dptr_1 = Cst;
  iptr_b = ADC3Bas_in->Bas_2p1h;
  for (int i1=0; i1 < ADC3Bas_in->Nbas_2p1h; ++i1) {
    im = iptr_b[0];
    iv = iptr_b[1];
    ig = iptr_b[2];
    
    iptr_b2 = ADC3Bas_in->Bas_2p1h;
    for (int i2=0; i2 < ADC3Bas_in->Nbas_2p1h; ++i2) {
      ia = iptr_b2[0];
      ib = iptr_b2[1];
      ir = iptr_b2[2];
      
      x1 = 0.0;
      if (ig==ir) x1 += V_Minnesota(WK_bas,im,iv,ia,ib);
      //
      if (im==ia) x1 += V_Minnesota(WK_bas,iv,ir,ig,ib);
      if (iv==ia) x1 -= V_Minnesota(WK_bas,im,ir,ig,ib);
      if (im==ib) x1 -= V_Minnesota(WK_bas,iv,ir,ig,ia);
      if (iv==ib) x1 += V_Minnesota(WK_bas,im,ir,ig,ia);
      
      dptr_1[i2] = x1;
      
      
      iptr_b2 += 3;
    }
    
    dptr_1[i1] += WK_bas->e_sp[im] + WK_bas->e_sp[iv] - WK_bas->e_sp[ig];
    
    iptr_b += 3;
    dptr_1 += LDC;
  }

  return 0;}


//
//  Builds the interaction (E^< + D) and coupling (N) matrices for the forward ADC(n) calculation (2h1p ISCs)
//
int Build_DysonADC_Matrix_2h1p(ADC3BasisK *ADC3Bas_in, double *Mst, double *Cst, int LDC ) {
  
  SpBasisK *WK_bas = ADC3Bas_in->SpBasLoc;
  
  int isp = ADC3Bas_in->iSpLoc;

  int im,iv, ig;
  double x1, x2;
  int  kin_x, kin_y, kin_z, sp_in, ch_in;
  int   kt_x,  kt_y,  kt_z, sp_t,  ch_t;
//  int   h5_x,  h5_y,  h5_z, sp_5,  ch_5,  h5,  h4;
  int   n7_x,  n7_y,  n7_z, sp_7,  ch_7,  n7,  n6, h5;
  
  
  int    *iptr_b = ADC3Bas_in->Bas_2h1p;
  double *dptr_1 = Mst;

  for (int i1=0; i1 < ADC3Bas_in->Nbas_2h1p; ++i1) {
    im = iptr_b[0];
    iv = iptr_b[1];
    ig = iptr_b[2];
    //   cout << " \n    bbb  " << i1 << flush;
    
    x1 = V_Minnesota(WK_bas,isp,ig,im,iv);
    *dptr_1 = x1;
    
    //
    //  Start with ADC(3) coupling vertices
    //
    kin_x = WK_bas->nx[isp];
    kin_y = WK_bas->ny[isp];
    kin_z = WK_bas->nz[isp];
    sp_in = WK_bas->spin[isp];
    ch_in = WK_bas->chrg[isp];
    
    kt_x = WK_bas->nx[ig]   + kin_x;
    kt_y = WK_bas->ny[ig]   + kin_y;
    kt_z = WK_bas->nz[ig]   + kin_z;
    sp_t = WK_bas->spin[ig] + sp_in;
    ch_t = WK_bas->chrg[ig] + ch_in;
    //
    x1 = 0.0;
    x2 = WK_bas->e_sp[im] + WK_bas->e_sp[iv];
    for (n6=WK_bas->N_holes; n6<WK_bas->SpNmax; ++n6) {
      n7_x = kt_x - WK_bas->nx[n6];
      n7_y = kt_y - WK_bas->ny[n6];
      n7_z = kt_z - WK_bas->nz[n6];
      sp_7 = sp_t - WK_bas->spin[n6];
      ch_7 = ch_t - WK_bas->chrg[n6];
      
      for (n7=n6+1; n7<WK_bas->SpNmax; ++n7) {
        
        if ( (n7_x == WK_bas->nx[n7]  ) && (n7_y == WK_bas->ny[n7]  ) && (n7_z == WK_bas->nz[n7]) &&
            (sp_7 == WK_bas->spin[n7]) && (ch_7 == WK_bas->chrg[n7]) ) {
          
          x1 += V_Minnesota(WK_bas,im,iv,n6,n7) * V_Minnesota(WK_bas,n6,n7,isp,ig)
          / (x2 - WK_bas->e_sp[n6] - WK_bas->e_sp[n7]);
        }
        
      }
    }
    //
    *dptr_1 += x1;
    
    kt_x = kin_x - WK_bas->nx[im];
    kt_y = kin_y - WK_bas->ny[im];
    kt_z = kin_z - WK_bas->nz[im];
    sp_t = sp_in - WK_bas->spin[im];
    ch_t = ch_in - WK_bas->chrg[im];
    //
    x1 = 0.0;
    x2 = WK_bas->e_sp[ig] - WK_bas->e_sp[iv];
    for (h5=0; h5<WK_bas->N_holes; ++h5) {
      n7_x = kt_x + WK_bas->nx[h5];
      n7_y = kt_y + WK_bas->ny[h5];
      n7_z = kt_z + WK_bas->nz[h5];
      sp_7 = sp_t + WK_bas->spin[h5];
      ch_7 = ch_t + WK_bas->chrg[h5];
      
      for (n7=WK_bas->N_holes; n7<WK_bas->SpNmax; ++n7) {
        
        if ( (n7_x == WK_bas->nx[n7]  ) && (n7_y == WK_bas->ny[n7]  ) && (n7_z == WK_bas->nz[n7]) &&
            (sp_7 == WK_bas->spin[n7]) && (ch_7 == WK_bas->chrg[n7]) ) {
          
          x1 += V_Minnesota(WK_bas,iv,h5,ig,n7) * V_Minnesota(WK_bas,im,n7,isp,h5)
          / (WK_bas->e_sp[h5] - WK_bas->e_sp[n7] - x2);
        }
        
      }
    }
    //
    *dptr_1 += x1;
    
    kt_x = kin_x - WK_bas->nx[iv];
    kt_y = kin_y - WK_bas->ny[iv];
    kt_z = kin_z - WK_bas->nz[iv];
    sp_t = sp_in - WK_bas->spin[iv];
    ch_t = ch_in - WK_bas->chrg[iv];
    //
    x1 = 0.0;
    x2 = WK_bas->e_sp[ig] - WK_bas->e_sp[im];
    for (h5=0; h5<WK_bas->N_holes; ++h5) {
      n7_x = kt_x + WK_bas->nx[h5];
      n7_y = kt_y + WK_bas->ny[h5];
      n7_z = kt_z + WK_bas->nz[h5];
      sp_7 = sp_t + WK_bas->spin[h5];
      ch_7 = ch_t + WK_bas->chrg[h5];
      
      for (n7=WK_bas->N_holes; n7<WK_bas->SpNmax; ++n7) {
        
        if ( (n7_x == WK_bas->nx[n7]  ) && (n7_y == WK_bas->ny[n7]  ) && (n7_z == WK_bas->nz[n7]) &&
            (sp_7 == WK_bas->spin[n7]) && (ch_7 == WK_bas->chrg[n7]) ) {
          
          x1 += V_Minnesota(WK_bas,im,h5,ig,n7) * V_Minnesota(WK_bas,iv,n7,isp,h5)
          / (WK_bas->e_sp[h5] - WK_bas->e_sp[n7] - x2);
        }
        
      }
    }
    //
    *dptr_1 -= x1;
    //
    //  end of ADC(3) coupling vertices
    //
    
    iptr_b += 3;
    dptr_1 += 1;
  }

  
  int ia, ib, ir;
  int *iptr_b2;
  //
  iptr_b = ADC3Bas_in->Bas_2h1p;
  dptr_1 = Cst;
  for (int i1=0; i1 < ADC3Bas_in->Nbas_2h1p; ++i1) {
    im = iptr_b[0];
    iv = iptr_b[1];
    ig = iptr_b[2];
    
    iptr_b2 = ADC3Bas_in->Bas_2h1p;
    for (int i2=0; i2 < ADC3Bas_in->Nbas_2h1p; ++i2) {
      ia = iptr_b2[0];
      ib = iptr_b2[1];
      ir = iptr_b2[2];
      
      x1 = 0.0;
      if (ig==ir) x1 -= V_Minnesota(WK_bas,im,iv,ia,ib);
      //
      if (im==ia) x1 -= V_Minnesota(WK_bas,iv,ir,ig,ib);
      if (iv==ia) x1 += V_Minnesota(WK_bas,im,ir,ig,ib);
      if (im==ib) x1 += V_Minnesota(WK_bas,iv,ir,ig,ia);
      if (iv==ib) x1 -= V_Minnesota(WK_bas,im,ir,ig,ia);
      
      dptr_1[i2] = x1;
      
      
      iptr_b2 += 3;
    }
    
    dptr_1[i1] += WK_bas->e_sp[im] + WK_bas->e_sp[iv] - WK_bas->e_sp[ig];
    
    iptr_b += 3;
    dptr_1 += LDC;
  }


  return 0;}

