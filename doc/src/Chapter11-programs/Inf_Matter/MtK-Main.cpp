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
//  MtK-Main.cpp  --  This is the Main function.
//


#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
using namespace std;

#include "MtK-Global_data.h"
#include "MtK-Bases.h"
#include "MtK-SpProppagator.h"
#include "MtK-Potential_minnesota.h"


int NLanczos = -100;


double GammaPlot_1 = 1.2; //  Width of lorentzians used to fold the spectral distribution in
double GammaPlot_2 = 7.0; // class SpctDist::  These are available through 'MtK-Global_data.h' and
                          // are initialised here.


SpctDist  *g_sp_out=NULL, *Sigma_irred=NULL;  // Initialise global pointers (from 'MtK-Global_data.h').

int main(int argc, char **argv) {

  //
  //  Reads the relevant parametes for the calculation and builds
  // the single-particle basis.
  //
  
  cout << " Welcome to MattersK v1.0 --  \n\n";
  
  int nsq_max, n_max, ch_mx, N_holes;
  double xe_hf,  rho_mt,  k_F, x1;
  
  cout << "\n\n Enter in order: max for nsq = sum_i ni^2 ,  n_max, and chrg_max (0 for PNM and 1for SNM) ?\n";
  
  cin >> nsq_max >> n_max >> ch_mx;


  cout << "\n\n Number of nucleons in the Fermi Sea? ";
  cin >> N_holes;
  
  cout << "\n\n Density [in fm^-3]? ";
  cin >> rho_mt;

  SpBasisK SpBas;
  
  k_F = 6.0 * PI * PI * rho_mt / double(2 + 2*ch_mx) ;
  k_F = pow( k_F , (1.0/3.0) );

  SpBas.Lbox = pow( (N_holes/rho_mt) , (1.0/3.0) );
  
  cout << "\n\n density = " << rho_mt << "/fm^3 ,     K_F = " << k_F << "/fm ,   L_box = " << SpBas.Lbox<< " fm\n";

  cout << "\n\n Now building the single particle basis:\n ---------------------------------------\n";

  SpBas.N_holes = N_holes;
  
  SpBas.Build_sp_basis(nsq_max,n_max,0,ch_mx);
  
  SpBas.Build_groups_tables();


  //
  //  Create the ADC3 basis and count its dimensions:
  //
  cout << "\n\n Statistics of the 2p1h and 2h1p configurations:\n -----------------------------------------------\n";
 
  ADC3BasisK ADC3Bas(&SpBas);
  
  ADC3Bas.Count_ALL_2p1h_2h1p_bases();

  
  

  //
  //  Run the HF bit:
  //
  cout << "\n\n Hartree-Fock solution:\n ----------------------\n";

  xe_hf = 0.0;
  for (int isp=0; isp<N_holes; ++isp) xe_hf += SpBas.e_kin[isp];
  cout << "\n    E_kin/A = "<< xe_hf/N_holes <<  " MeV   (kinetic energy of free Fermi gas)";

  for (int itr=0; itr<2; ++itr) { // In reality, only one iteration is necessary since k-basis is already HF
    for (int isp=0; isp<SpBas.SpNmax; ++isp) SpBas.e_HF[isp] = HF_Minn_MtxEls(&SpBas,isp,isp,N_holes);
    xe_hf = 0.0;
    x1 = 0.0;
    for (int isp=0; isp<N_holes; ++isp) {
      xe_hf += (SpBas.e_kin[isp] + SpBas.e_HF[isp]); // Koltun contribuiton
      x1 +=  SpBas.e_kin[isp];  // Kinetic en.
    }
    //cout << "\n   Itr.  " << itr << "    --> E_HF = "<< xe_hf/2.0/N_holes;
    //cout << " MeV        E_kin = " << x1/N_holes << " MeV";
  }
  cout << "\n    E_HF/A  = "<< xe_hf/2.0/N_holes <<  " MeV   (uncorrelated HF energy)";

  cout << endl;

  char SFplot;
  cout << flush;
  cout << "\n\n Would you like to generate data files for plotting the spectral function and self energy [y/Y or n/N] ? ";
  cin >> SFplot;
  
  if (('y' == SFplot) || ('Y' == SFplot)) {
    //
    // Creates data files with the spectral distribution in the HF approximation:
    //
    Plot_HF_dist(&SpBas);
  }


  // For saving the calculate self-energy and propagators, see 'MtK-Global_data.h'.
  Sigma_irred = new SpctDist(&SpBas);
  g_sp_out    = new SpctDist(&SpBas);


  cout << "\n\n\n START RUNNING THE DYSON DIAGONALIZATIONS:\n -------------------------------------------- \n\n";

  NLanczos = -100;
  cout << "\n\n How many Lanczos iterations [< 0 for no Lanczos] ? ";
  cin >> NLanczos;
  

  double xA, xE, xK; // Sum rule contributions to # of particles, kin. energy and Koltun SR
  
  // This sets the uncorrelated s.p. energies (we use HF here, this is not sc0...)
  for (int isp=0; isp<SpBas.SpNmax; ++isp) SpBas.e_sp[isp] = SpBas.e_HF[isp];

  // Solve the Dyson equation; the resulting propagator and self-energy are saved
  // into object 'SpctDist *g_sp_out, *Sigma_irred;' which are global to the code
  // and made available thorugh 'MtK-Global_data.h'.
  Solve_Dyson(&ADC3Bas, &xA, &xE, &xK);


  cout << "\n Density,  E_HF,  corr en   (corr. en rescaled to A="<<xA<<" particles) :"
       << rho_mt << "    " << xe_hf << "    " << (xE/double(N_holes)) - xe_hf << "    ("
       <<  (xE/xA) - xe_hf << ")\n" << flush;
  
  if (('y' == SFplot) || ('Y' == SFplot)) {
    //
    // Write results for spectral function and self-energy to data file.
    //
    cout << "\n\n Generating datafiles with spectral fucnction and self-energy...\n\n";
    //
    for (int isp=0; isp<SpBas.SpNmax; ++isp)
    if (isp == SpBas.gr_rep[SpBas.group[isp]]) {
      g_sp_out->Plot_SpectFnct_vs_E(isp);
      Sigma_irred->Plot_SelfEn_vs_E(isp);
    }
    //
    // data file for the 3D plot
    g_sp_out->Plot_SpectFnct_3D();
  }

  return 0;}


//
// Generate a datafile with the spectral distribution in the HF approximation:
//
int Plot_HF_dist(SpBasisK* Bas) {


  char filename[100];

  sprintf(filename, "SpectFunctHF_3D_parts.dat");
  ofstream file_qp(filename, ios::trunc|ios::out);

  file_qp << "#\n#\n#   Lonertzian-smoothed PARTICLE spectral function (both aprticle and hole) as\n";
  file_qp << "#  function of momentum and energy\n";
  file_qp << "#\n#\n#    k              E               A(k,E)\n";


  sprintf(filename, "SpectFunctHF_3D_holes.dat");
  ofstream file_qh(filename, ios::trunc|ios::out);

  file_qh << "#\n#\n#   Lonertzian-smoothed HOLE spectral function (both aprticle and hole) as\n";
  file_qh << "#  function of momentum and energy\n";
  file_qh << "#\n#\n#    k              E               A(k,E)\n";

  double GmHF = .4 , xe, str;

  for (int isp=0; isp<Bas->SpNmax; ++isp) {

    if (isp != Bas->gr_rep[Bas->group[isp]]) continue;

    if (isp < Bas->N_holes) {
      file_qh << endl;
      file_qh << Bas->k[isp] <<"    "<< xe+Bas->e_HF[isp] << "       "  << 0.0 <<endl;
      file_qh << Bas->k[isp] <<"    "<< xe+Bas->e_HF[isp] << "       "  << 1.0/PI/GammaPlot_1 <<endl;
      file_qh << Bas->k[isp] <<"    "<< xe+Bas->e_HF[isp] << "       "  << 0.0 <<endl;
      file_qh << endl;
    } else {
      file_qp << endl;
      file_qp << Bas->k[isp] <<"    "<< xe+Bas->e_HF[isp] << "       "  << 0.0 <<endl;
      file_qp << Bas->k[isp] <<"    "<< xe+Bas->e_HF[isp] << "       "  << 1.0/PI/GammaPlot_1 <<endl;
      file_qp << Bas->k[isp] <<"    "<< xe+Bas->e_HF[isp] << "       "  << 0.0 <<endl;
      file_qp << endl;
    }

  }

  file_qp.close();
  file_qh.close();


  return 0;
}




