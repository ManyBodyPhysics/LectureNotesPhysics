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
//  MtK-Bases.cpp  --  classes for the s.-p. and 2p1h/2h1p configurations.
//


#include <cmath>
#include <iostream>

using namespace std;

#include "MtK-Physics_constants.h"
#include "MtK-Bases.h"


//
//   Creators and destructors
//   ------------------------
//

SpBasisK::SpBasisK(void ) {
  //this->Init_NULL();
  Lbox = 1.0;
  return;}


SpBasisK::~SpBasisK(void ) {
  // Must free allthe allocated memory
//  this->Free_mem();
//  this->Init_NULL();
  return;}

// The rationale of Init_NULL and Free_mem si that
// they are written only once and called when needed
// so changes are always done in one place only.

// Initialise parameters to an 'empty object'

int SpBasisK::Count_sp_basis(int nsq_mn, int nsq_mx, int imax) {
    
  int i_count, itest;
  //int imax = int(sqrt(double(nsq_mx))) + 1;
  //
  // first count the k configurations:
    i_count = 0;
    for (int ix=-imax; ix<=imax; ++ix) {
      for (int iy=-imax; iy<=imax; ++iy) {
        for (int iz=-imax; iz<=imax; ++iz) {
          itest = ix*ix + iy*iy + iz*iz;
          if ((nsq_mn <= itest) && (itest <= nsq_mx)) ++i_count;
        }

      }

    //cout << "   " << isq << "    ( " << sqrt(double(isq))<<" ) -    ";
    //if (i_count> 0) cout  << i_count;
    //cout << endl;
  }

  return i_count;}


void SpBasisK::Build_sp_basis(int nsq_mx, int imax, int ch_min_in, int ch_max_in) {

  ch_min = ch_min_in;
  ch_max = ch_max_in;
  
  // First count the k configurations:
  //
  int i1 = this->Count_sp_basis(0, nsq_mx, imax);
 
  SpNAlloc = (ch_max - ch_min + 1) * 2 * i1 + 10; // +10 for safety
  
  cout << "\n allocating space for "<< SpNAlloc << " sp states... \n";
  
  nx    = new int[SpNAlloc];
  ny    = new int[SpNAlloc];
  nz    = new int[SpNAlloc];
  nsq   = new int[SpNAlloc];
  spin  = new int[SpNAlloc];
  chrg  = new int[SpNAlloc];
  k     = new double[SpNAlloc];
  e_kin = new double[SpNAlloc];
  e_HF  = new double[SpNAlloc];
  e_sp  = new double[SpNAlloc];
  group = new int[SpNAlloc];
  
  int ich, is;

  double xk = 0.0, xek = 0.0;
  
  i1 = 0;
  for (int isq=0; isq<=nsq_mx; ++isq) {
    for (int ix=-imax; ix<=imax; ++ix) {
      for (int iy=-imax; iy<=imax; ++iy) {
        for (int iz=-imax; iz<=imax; ++iz) {
          if ((ix*ix + iy*iy + iz*iz) != isq) continue;

          xek = double(isq);
          xk  = sqrt(xek) * 2.0 * PI / Lbox;
          xek = xk * xk * MeVfm * MeVfm / 2.0 / NUCLEONmass;  //nuc_mass_ave;
          cout << i1 << "  " << ix << "  " << iy << "  " << iz << "  ";
          cout << isq << "  " << xk << "  " << xek << endl;

          for (ich=ch_min; ich<=ch_max; ++ich)
            for (is=-1; is<2; is+=2) {
              nx[i1] = ix;
              ny[i1] = iy;
              nz[i1] = iz;
             nsq[i1] = isq;
            spin[i1] = is;
            chrg[i1] = ich;
               k[i1] = xk;
           e_kin[i1] = xek;
            e_HF[i1] = 0.0;
            e_sp[i1] = 0.0;
           group[i1] = -100;
            ++i1;
            }
        
        }
      }
    }
  }
  SpNmax = i1;
  cout << "\n\n The total dimension of the sp basis is " << SpNmax << endl;

  return;}


int SpBasisK::Build_groups_tables(void ) {
  
  //
  // Maximum value for nmax
  //
  int AbsN_mx = 0;
  for (int isp=0; isp<this->SpNmax; ++isp) {
    AbsN_mx = ( AbsN_mx > abs(nx[isp]) ) ? AbsN_mx : abs(nx[isp]);
    AbsN_mx = ( AbsN_mx > abs(ny[isp]) ) ? AbsN_mx : abs(ny[isp]);
    AbsN_mx = ( AbsN_mx > abs(nz[isp]) ) ? AbsN_mx : abs(nz[isp]);
  }
  
  //
  // Get the upper limit for the number of grroups
  //
  N_grps = (AbsN_mx+1)*(AbsN_mx+2)*(AbsN_mx+3)/6;
  
  gr_mlt = new int[N_grps];
  gr_rep = new int[N_grps];

  
  //
  //  Unset the group index in the sp basis, to check later that everything
  // has been associated to a group
  //
  for (int isp=0; isp<this->SpNmax; ++isp) group[isp] = -100;

  
  int chrg_loc, spin_loc, i1, i2, i3, i_mult, i_rep, n1, n2, n3, isp, itmp;
  int count=0;
  //for (chrg_loc=ch_min; chrg_loc<=ch_max; ++chrg_loc)
    //for (spin_loc=-1; spin_loc<2; spin_loc+=2)
      for (i1=0; i1<=AbsN_mx; ++i1)
        for (i2=i1; i2<=AbsN_mx; ++i2)
          for (i3=i2; i3<=AbsN_mx; ++i3) {
            
            i_mult = 0;
            i_rep = -100;
            for (isp=0; isp<this->SpNmax; ++isp) {
              //if ((chrg_loc != chrg[isp]) || (spin_loc != spin[isp])) continue;
              
              n1 = abs(nx[isp]);
              n2 = abs(ny[isp]);
              n3 = abs(nz[isp]);
              
              if (n1 > n2) {itmp=n1; n1=n2; n2=itmp;}
              if (n1 > n3) {itmp=n1; n1=n3; n3=itmp;}
              if (n2 > n3) {itmp=n2; n2=n3; n3=itmp;}
              
              if ((n1==i1) && (n2==i2) && (n3==i3)) {
                ++i_mult;
                if (i_rep < 0) i_rep = isp;
                group[isp] = count;
              }
              
            }
            
            gr_mlt[count] = i_mult;
            gr_rep[count] = i_rep;
            
            ++count;
            
          }
  
  //
  //  Now checks that multiplicities are correct and the each orbit has a group
  //
  if (count != N_grps) {cout << " ERROR (count != N_grps) "<<count<<"   "<<N_grps<<"\n"; exit(100);}
  //

  cout << "\n\n Symmetry groups:\n n.  --  n_tot   mult   i_rep\n";
  count = 0;
  for(i1=0; i1<N_grps; ++i1) {
    count += gr_mlt[i1];
    cout << " " << i1 << "   --    "<< count << "      " << gr_mlt[i1] << "      " << gr_rep[i1] <<endl;
    if ( (0 < gr_mlt[i1]) && ( (0 > gr_rep[i1]) || (SpNmax <= gr_rep[i1]) ) ) {cout << "BBB\n"; exit(100);}
  }
  cout << "\n";
  if (count != SpNmax) {cout << " ERROR (count != SpNmax) "<<count<<"   "<<SpNmax<<"\n"; exit(100);}
  //
  for (int isp=0; isp<this->SpNmax; ++isp) if (0 > group[isp]) {cout << "DDD\n"; exit(100);}
  
  return 0;}


//
//  2p1h and 2h1p intermediate state configurations (ISC)
//

ADC3BasisK::ADC3BasisK(SpBasisK *bas_in) {

  SpBasLoc = bas_in;

  //
  //this->Init_NULL();
  iSpLoc = -100;
  //
  Bas_2p1h = NULL;
  Bas_2h1p = NULL;
  //
  Nbas_2p1h = -100;
  Nbas_2h1p = -100;
  //
  return;
}

ADC3BasisK::~ADC3BasisK(void ) {
  // Must free allthe allocated memory
  //  this->Free_mem();
  if (NULL != Bas_2p1h) delete [] Bas_2p1h;  Bas_2p1h = NULL;
  if (NULL != Bas_2h1p) delete [] Bas_2h1p;  Bas_2h1p = NULL;

  //  this->Init_NULL();
  iSpLoc = -100;
  //
  Bas_2p1h = NULL;
  Bas_2h1p = NULL;
  //
  Nbas_2p1h = -100;
  Nbas_2h1p = -100;
  return;}


void ADC3BasisK::Count_ALL_2p1h_2h1p_bases( ) {
  
  int  i_pph, i_hhp, i_pph_mx, i_hhp_mx, i_tot, i_tot_mx;
  int sum_pph, sum_hhp;

  bool GrpsFlag = true;

  i_pph_mx = 0;
  i_hhp_mx = 0;
  i_tot_mx = 0;
  sum_pph  = 0;
  sum_hhp  = 0;
  
  for (int ia=0; ia<SpBasLoc->SpNmax; ++ia) {

    if (GrpsFlag) if (ia != this->SpBasLoc->gr_rep[this->SpBasLoc->group[ia]]) continue;
    
    Count_2p1h_2h1p_basis(ia, &i_pph, &i_hhp);

    sum_pph += i_pph;
    sum_hhp += i_hhp;
    
    i_tot = 1 + i_pph + i_hhp;
    
    i_pph_mx = (i_pph_mx > i_pph) ? i_pph_mx : i_pph;
    i_hhp_mx = (i_hhp_mx > i_hhp) ? i_hhp_mx : i_hhp;
    i_tot_mx = (i_tot_mx > i_tot) ? i_tot_mx : i_tot;
    
    cout << ia << "          " << i_pph << "          " << i_hhp << "          " << i_tot << "\n";
    
  }
  
  cout << "\n\n Max block sizes:     " << i_pph_mx << " (2p1h)         " << i_hhp_mx << " (2h1p)         " << i_tot_mx << " (full Dys mtx)\n";
  
  cout << "\n\n Total # of 2p1h and  2h1p states:  " << sum_pph << "     " << sum_hhp << endl;


  return;}

void ADC3BasisK::Count_2p1h_2h1p_basis(int ia, int *i_pph, int *i_hhp) {
  
  int    ka_x, ka_y, ka_z, sp_a, ch_a;
  int    p1_x, p1_y, p1_z, sp_1, ch_1, p1;
  int    p2_x, p2_y, p2_z, sp_2, ch_2, p2, h3;
  
  
    ka_x = SpBasLoc->nx[ia];
    ka_y = SpBasLoc->ny[ia];
    ka_z = SpBasLoc->nz[ia];
    sp_a = SpBasLoc->spin[ia];
    ch_a = SpBasLoc->chrg[ia];
    //cout << "ia=" << ia << endl<<flush;
    
    *i_pph = 0;
    for (p1=SpBasLoc->N_holes; p1<SpBasLoc->SpNmax; ++p1) {
      p1_x = SpBasLoc->nx[p1]   - ka_x;
      p1_y = SpBasLoc->ny[p1]   - ka_y;
      p1_z = SpBasLoc->nz[p1]   - ka_z;
      sp_1 = SpBasLoc->spin[p1] - sp_a;
      ch_1 = SpBasLoc->chrg[p1] - ch_a;
      //cout << "ia, p1=" << ia <<"   " << p1<< endl<<flush;
      
      for (p2=p1+1; p2<SpBasLoc->SpNmax; ++p2) {
        p2_x = SpBasLoc->nx[p2]   + p1_x;
        p2_y = SpBasLoc->ny[p2]   + p1_y;
        p2_z = SpBasLoc->nz[p2]   + p1_z;
        sp_2 = SpBasLoc->spin[p2] + sp_1;
        ch_2 = SpBasLoc->chrg[p2] + ch_1;
        //cout << "ia, p1, p2=" << ia <<"   " << p1 << "   " << p2 << endl<<flush;
        
        for (h3=0; h3<SpBasLoc->N_holes; ++h3) {
          if ( (p2_x != SpBasLoc->nx[h3]  ) || (p2_y != SpBasLoc->ny[h3]  ) || (p2_z != SpBasLoc->nz[h3]) ||
               (sp_2 != SpBasLoc->spin[h3]) || (ch_2 != SpBasLoc->chrg[h3]) ) continue;
          
          //cout << "ia, p1, p2, h3=" << ia <<"   " << p1 << "   " << p2 << "  " << h3<< endl<<flush;
          
          ++(*i_pph);
          
        }
      }
    }
    
    *i_hhp = 0;
    for (p1=0; p1<SpBasLoc->N_holes; ++p1) {
      p1_x = SpBasLoc->nx[p1]   - ka_x;
      p1_y = SpBasLoc->ny[p1]   - ka_y;
      p1_z = SpBasLoc->nz[p1]   - ka_z;
      sp_1 = SpBasLoc->spin[p1] - sp_a;
      ch_1 = SpBasLoc->chrg[p1] - ch_a;
      
      for (p2=p1+1; p2<SpBasLoc->N_holes; ++p2) {
        p2_x = SpBasLoc->nx[p2]   + p1_x;
        p2_y = SpBasLoc->ny[p2]   + p1_y;
        p2_z = SpBasLoc->nz[p2]   + p1_z;
        sp_2 = SpBasLoc->spin[p2] + sp_1;
        ch_2 = SpBasLoc->chrg[p2] + ch_1;
        
        for (h3=SpBasLoc->N_holes; h3<SpBasLoc->SpNmax; ++h3) {
          if ( (p2_x != SpBasLoc->nx[h3]  ) || (p2_y != SpBasLoc->ny[h3  ]) || (p2_z != SpBasLoc->nz[h3]) ||
               (sp_2 != SpBasLoc->spin[h3]) || (ch_2 != SpBasLoc->chrg[h3]) ) continue;
          
          ++(*i_hhp);
          
        }
      }
    }
  
  return;}

void ADC3BasisK::Build_2p1h_2h1p_basis(int ia) {
  
  int i1, i2, *iptr_b;

  this->iSpLoc = ia; // for use by other functions that exploit this bases

  this->Count_2p1h_2h1p_basis( ia,  &i1, &i2);

  if (NULL != Bas_2p1h) delete [] Bas_2p1h;
  if (NULL != Bas_2h1p) delete [] Bas_2h1p;

  Nbas_2p1h = i1;
  Nbas_2h1p = i2;

  Bas_2p1h = new int[3*(Nbas_2p1h+4)];  // +4 is just for safety...
  Bas_2h1p = new int[3*(Nbas_2h1p+4)];  // +4 is just for safety...

  int    ka_x, ka_y, ka_z, sp_a, ch_a;
  int    p1_x, p1_y, p1_z, sp_1, ch_1, p1;
  int    p2_x, p2_y, p2_z, sp_2, ch_2, p2, h3;
  
  
  ka_x = SpBasLoc->nx[ia];
  ka_y = SpBasLoc->ny[ia];
  ka_z = SpBasLoc->nz[ia];
  sp_a = SpBasLoc->spin[ia];
  ch_a = SpBasLoc->chrg[ia];



  iptr_b = Bas_2p1h;
  //
  for (p1=SpBasLoc->N_holes; p1<SpBasLoc->SpNmax; ++p1) {
    p1_x = SpBasLoc->nx[p1]   - ka_x;
    p1_y = SpBasLoc->ny[p1]   - ka_y;
    p1_z = SpBasLoc->nz[p1]   - ka_z;
    sp_1 = SpBasLoc->spin[p1] - sp_a;
    ch_1 = SpBasLoc->chrg[p1] - ch_a;
    
    for (p2=p1+1; p2<SpBasLoc->SpNmax; ++p2) {
      p2_x = SpBasLoc->nx[p2]   + p1_x;
      p2_y = SpBasLoc->ny[p2]   + p1_y;
      p2_z = SpBasLoc->nz[p2]   + p1_z;
      sp_2 = SpBasLoc->spin[p2] + sp_1;
      ch_2 = SpBasLoc->chrg[p2] + ch_1;
      
      for (h3=0; h3<SpBasLoc->N_holes; ++h3) {
        if ( (p2_x != SpBasLoc->nx[h3]  ) || (p2_y != SpBasLoc->ny[h3]  ) || (p2_z != SpBasLoc->nz[h3]) ||
            (sp_2 != SpBasLoc->spin[h3]) || (ch_2 != SpBasLoc->chrg[h3]) ) continue;
        
        iptr_b[0] = p1;
        iptr_b[1] = p2;
        iptr_b[2] = h3;
        iptr_b += 3;
        
      }
    }
  }
  if ( (iptr_b - Bas_2p1h) != 3*Nbas_2p1h) {cout << " \n\n ERROR 2p1h xxx!!!!!"; exit(100);}


  iptr_b = Bas_2h1p;
  for (p1=0; p1<SpBasLoc->N_holes; ++p1) {
    p1_x = SpBasLoc->nx[p1]   - ka_x;
    p1_y = SpBasLoc->ny[p1]   - ka_y;
    p1_z = SpBasLoc->nz[p1]   - ka_z;
    sp_1 = SpBasLoc->spin[p1] - sp_a;
    ch_1 = SpBasLoc->chrg[p1] - ch_a;
    
    for (p2=p1+1; p2<SpBasLoc->N_holes; ++p2) {
      p2_x = SpBasLoc->nx[p2]   + p1_x;
      p2_y = SpBasLoc->ny[p2]   + p1_y;
      p2_z = SpBasLoc->nz[p2]   + p1_z;
      sp_2 = SpBasLoc->spin[p2] + sp_1;
      ch_2 = SpBasLoc->chrg[p2] + ch_1;
      
      for (h3=SpBasLoc->N_holes; h3<SpBasLoc->SpNmax; ++h3) {
        if ( (p2_x != SpBasLoc->nx[h3]  ) || (p2_y != SpBasLoc->ny[h3  ]) || (p2_z != SpBasLoc->nz[h3]) ||
            (sp_2 != SpBasLoc->spin[h3]) || (ch_2 != SpBasLoc->chrg[h3]) ) continue;
        
        iptr_b[0] = p1;
        iptr_b[1] = p2;
        iptr_b[2] = h3;
        iptr_b += 3;
        
      }
    }
  }
  if ( (iptr_b - Bas_2h1p) != 3*Nbas_2h1p) {cout << " \n\n ERROR 2h1p xxx!!!!!"; exit(100);}

  return;}

