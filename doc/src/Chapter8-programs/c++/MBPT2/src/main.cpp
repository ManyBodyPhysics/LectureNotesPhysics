// CCD code for arbitrary basis
// Written by Justin Lietz starting in June 2015


#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <cblas.h>
#include <omp.h>
#include "SymBlock.hpp"
#include "abstractSPbasis.hpp"
#include "infMatterSPBasis.hpp"

using namespace std;

// Two Body operator using channels
struct TB_OpChannels{ 
  SymmetryBlock * hhpp;   
};

int hhpp_index(int i, int j, int a, int b, abstractSPbasis * SPbasis);

// Pass vnn by ref to not create an image of the whole thing.
// vnn is a struct, so pass with ampersand
// SPbasis is a pointer to a class (abstract SPbasis)
void load_hhpp_blocks(TB_OpChannels &vnn, abstractSPbasis * SPbasis);
void load_full_hhpp(double * vnn_full_hhpp, abstractSPbasis * SPbasis);
double MBPT2corr(TB_OpChannels &vnn, int numChannels, abstractSPbasis * SPbasis);
double mbpt2_full_hhpp(double * vnn_full_hhpp, abstractSPbasis * SPbasis);
double referenceEnergy(abstractSPbasis * SPbasis);
double calcVpqrs(int p, int q, int r, int s, abstractSPbasis * SPbasis);

///////////////////////////////////////////////////
// begin main program
///////////////////////////////////////////////////


int main(int argc, char * argv[])
{  
  double density;
  int tzMax,shellMax,Nparticles;
  
  abstractSPbasis * SPbasis;
  
  // select basis and basis parameters at runtime
  // move this to a config file at some point?
  if( argc != 6){    
    cout << "or : 1 density tzMax shellMax Nparticles: values at command line" << endl;   
    return 0;  
  } else if ( atoi(argv[1]) == 1 ){    
    density = atof(argv[2]);
    tzMax = atoi(argv[3]);
    shellMax = atoi(argv[4]);
    Nparticles = atoi(argv[5]);   
    SPbasis = new infMatterSPBasis (density,tzMax,shellMax,Nparticles);  
  } else {
    cout << "bases: 0 for pairing, 1 for infMatter, 2 for electronGas, 3 for 2d qdots." << endl;
    return 0;
  }
  
  ///////////////////////////////////////////////////////////////
  // Set Up Basis
  ///////////////////////////////////////////////////////////////


  int fermiLevel;
  fermiLevel = Nparticles;
  // Initialize single particle basis  
  SPbasis->generateIndexMap();
  SPbasis->generateBasis();  
  int Nspstates = SPbasis->Nspstates; 
  //Check build
  SPbasis->printBasis(); 
  // calculate reference energy BEFORE
  // rotating into new basis
  double Eref;
  Eref = referenceEnergy(SPbasis); 
  cout << "Reference Energy: " << Eref << endl;
  // rotate sp energies to HF energies
  double ei,ea;
  for(int i=0; i<fermiLevel; i++){
    ei = SPbasis->spEnergy[i];
    for(int j=0; j<fermiLevel; j++){
      ei += calcVpqrs(i,j,i,j,SPbasis);
    }
    SPbasis->spEnergy[i] = ei;
  } 
  for(int a=fermiLevel; a<Nspstates; a++){
    ea = SPbasis->spEnergy[a];
    for(int j=0; j<fermiLevel; j++){
      ea += calcVpqrs(a,j,a,j,SPbasis);
    }
    SPbasis->spEnergy[a] = ea;
  }

  //Check rotation 
  SPbasis->printBasis();


  ///////////////////////////////////////////////////////
  // MBPT2 calculation using block matrices
  ///////////////////////////////////////////////////////
  

  // |i,j> are the two body states
  // set up the indexing scheme for these
  cout << endl << "CHANNELS" << endl;  
  double basisStart = omp_get_wtime();
  SPbasis->setUpTwoStateChannels();  
  double basisEnd = omp_get_wtime();    
  int numChannels = SPbasis->Nchannels;   
  cout << "Two Body Channels setup time: " << basisEnd - basisStart << endl;

  double loadBlockStart = omp_get_wtime();
  TB_OpChannels vnn; 
  vnn.hhpp = new SymmetryBlock[numChannels];   
  load_hhpp_blocks(vnn,SPbasis);  
  double loadBlockEnd = omp_get_wtime();  
  cout << "Block matrix Load time: " << loadBlockEnd - loadBlockStart <<  endl;
  double corrMBPT2; 
  double multTime = omp_get_wtime();
  corrMBPT2 = MBPT2corr(vnn, numChannels, SPbasis);
  double wallTimeEnd = omp_get_wtime();
  cout << "Block matrix mult time: " << wallTimeEnd - multTime << endl;
  cout << "Total Block calc time: " << wallTimeEnd - basisStart << endl;
  cout << "MBPT2corr/A = " << setprecision(16) << corrMBPT2/Nparticles << endl << endl;
  
  /////////////////////////////////////////////////////////////////////////
  // MBPT2 Calculation using the full 4-index structure
  /////////////////////////////////////////////////////////////////////////

  double fullLoadStart = omp_get_wtime();       
  double * vnn_hhpp_full = new double[Nparticles*Nparticles*(Nspstates-Nparticles)*(Nspstates-Nparticles)];
  load_full_hhpp(vnn_hhpp_full,SPbasis);
  double fullLoadEnd = omp_get_wtime();  
  cout << "full matrix load time: " << fullLoadEnd - fullLoadStart << endl;

  double fullCalcStart = omp_get_wtime();
  double corrMBPT2_full = mbpt2_full_hhpp(vnn_hhpp_full, SPbasis);
  double fullCalcEnd = omp_get_wtime();
  cout << "full struct mult time: " << fullCalcEnd - fullCalcStart << endl;  
  cout << "Total full struct calc time: " << fullCalcEnd - fullLoadStart << endl;
  cout << "corrMBPT2_full/a = " <<  corrMBPT2_full/Nparticles << endl << endl;

  // free memory
  for(int ichan = 0; ichan < numChannels; ichan++){   
    vnn.hhpp[ichan].deallocate();       
  }
 
  delete [] vnn.hhpp;
    
  SPbasis->deallocate();   

  return 0;  
}


// calculate the 2nd Order Many-Body Perturbation Theory correlation energy 
double MBPT2corr(TB_OpChannels &vnn, int numChannels, abstractSPbasis * SPbasis)
{
  double corr = 0.;  
  double energyDenom;
  int h1,h2,p1,p2;
  for(int ichan = 0; ichan < numChannels; ichan++){    
    for(int hhIndex = 0; hhIndex < vnn.hhpp[ichan].getRowNum(); hhIndex++){
      for(int ppIndex = 0; ppIndex < vnn.hhpp[ichan].getColNum(); ppIndex++){
	h1 = vnn.hhpp[ichan].rowMap[hhIndex][0];
	h2 = vnn.hhpp[ichan].rowMap[hhIndex][1];
	p1 = vnn.hhpp[ichan].colMap[ppIndex][0];
	p2 = vnn.hhpp[ichan].colMap[ppIndex][1];
	energyDenom = SPbasis->spEnergy[h1] + SPbasis->spEnergy[h2] - SPbasis->spEnergy[p1] - SPbasis->spEnergy[p2];
	corr += vnn.hhpp[ichan].getElement(hhIndex,ppIndex)*vnn.hhpp[ichan].getElement(hhIndex,ppIndex)/energyDenom;	
      }
    }
  }
  corr = 0.25*corr;
  return corr;
} // end MBPT2corr

double mbpt2_full_hhpp(double * vnn_full_hhpp, abstractSPbasis * SPbasis)
{
  double energyDenom,corr;
  int fermiLevel = SPbasis->Nparticles;
  int NparticleStates = SPbasis->Nspstates - SPbasis->Nparticles;
  for(int i=0; i<fermiLevel; i++){
    for(int j=0; j<fermiLevel; j++){
      for(int a=0; a<NparticleStates; a++){
	for(int b=0; b<NparticleStates; b++){	
	  energyDenom = SPbasis->spEnergy[i] + SPbasis->spEnergy[j] - SPbasis->spEnergy[fermiLevel+a] - SPbasis->spEnergy[fermiLevel+b];
	  corr += vnn_full_hhpp[hhpp_index(i,j,a,b,SPbasis)]*vnn_full_hhpp[hhpp_index(i,j,a,b,SPbasis)]/energyDenom;
	}
      }
    }
  }
  corr = 0.25*corr;
  return corr;
} // end load_full_hhpp


double referenceEnergy(abstractSPbasis * SPbasis)
{  
  int Nparticles = SPbasis->Nparticles;
  int fermiLevel = Nparticles;  
  double vijij = 0.;
  double Eref = 0.;
  
  for(int i = 0; i < fermiLevel-1; i++){
    for(int j = i+1; j < fermiLevel; j++){
      if( SPbasis->checkSympqrs(i,j,i,j) == 1 ){
  	vijij = calcVpqrs(i,j,i,j,SPbasis);
	Eref += vijij;
      } // end if  
    } // end j    
  } // end i
  
  for(int i = 0; i < fermiLevel; i++){
    Eref += SPbasis->spEnergy[i];
  }
  
  return Eref;
}

// Could possibly redefine SPbais->calc_TBME to already be
// antisymmetric
double calcVpqrs(int p, int q, int r, int s, abstractSPbasis * SPbasis)
{
  double vpqrs = 0.;
  vpqrs = SPbasis->calc_TBME(p,q,r,s) - SPbasis->calc_TBME(p,q,s,r);
  return vpqrs;
} // end calcVpqrs

int hhpp_index(int i, int j, int a, int b, abstractSPbasis * SPbasis){
  int N_h = SPbasis->Nparticles;
  int N_p = SPbasis->Nspstates-SPbasis->Nparticles;
  int index;
  index = i*N_h*N_p*N_p + j*N_p*N_p + a*N_p + b;
  return index;
}

void load_full_hhpp(double * vnn_full_hhpp, abstractSPbasis * SPbasis)
{
  int fermiLevel = SPbasis->Nparticles;
  int NparticleStates = SPbasis->Nspstates - SPbasis->Nparticles;
  for(int i=0; i<fermiLevel; i++){
    for(int j=0; j<fermiLevel; j++){
      for(int a=0; a<NparticleStates; a++){
	for(int b=0; b<NparticleStates; b++){
	  vnn_full_hhpp[hhpp_index(i,j,a,b,SPbasis)] = calcVpqrs(i,j,fermiLevel+a,fermiLevel+b,SPbasis);	 
	}
      }
    }
  }  
} // end load_full_hhpp

void load_hhpp_blocks(TB_OpChannels &vnn, abstractSPbasis * SPbasis)
{  
  int numChannels = SPbasis->Nchannels;
        
  ////////////////////////////////////////////////////////////////////////
  //cout << "Setting up vnn.hhpp" << endl;
  ////////////////////////////////////////////////////////////////////////
   
   
   
  //#pragma omp parallel
  //  {
  //#pragma omp for
    for(int ichan = 0; ichan < numChannels; ichan++){
      double vpqrs;
      int ppDim,hhDim;
      int p,q,r,s;
      
      hhDim = SPbasis->chanValue[ichan].hhDim;
      ppDim = SPbasis->chanValue[ichan].ppDim;
                     

      ////////////////////////////////////////////////////////////////////////
      //cout << "Setting up vnn.hhpp" << endl;
      ////////////////////////////////////////////////////////////////////////
      
      // actually need rowmap for hhpp
      vnn.hhpp[ichan].allocateWithLabels(hhDim,ppDim,2,2);
      vnn.hhpp[ichan].zeros();

      for(int ihh = 0; ihh < hhDim; ihh++){
      	for(int ipp = 0; ipp < ppDim; ipp++){
      	  p = SPbasis->chanValue[ichan].hhMap[ihh][0];
      	  q = SPbasis->chanValue[ichan].hhMap[ihh][1];
      	  r = SPbasis->chanValue[ichan].ppMap[ipp][0];
      	  s = SPbasis->chanValue[ichan].ppMap[ipp][1];
	  vnn.hhpp[ichan].rowMap[ihh][0] = p;
	  vnn.hhpp[ichan].rowMap[ihh][1] = q;
	  vnn.hhpp[ichan].colMap[ipp][0] = r;
	  vnn.hhpp[ichan].colMap[ipp][1] = s;
      	  vpqrs = calcVpqrs(p,q,r,s,SPbasis);
      	  vnn.hhpp[ichan].setElement(ihh,ipp,vpqrs);
      	} // end ihh loop
      } // end ipp loop
    } // end ichan
     
    //  } // end #pragma omp parallel
} // end load_hhpp_blocks
