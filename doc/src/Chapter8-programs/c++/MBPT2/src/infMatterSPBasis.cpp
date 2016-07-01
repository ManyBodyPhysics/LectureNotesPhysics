// Copyright (c) 2015-2016, Justin Gage Lietz
// All rights reserved.

// Need to write a PDF about this biz.

#include "infMatterSPBasis.hpp"
#include "util.hpp"
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cmath>
#include <omp.h>

infMatterSPBasis::infMatterSPBasis(double densityIn, int tzMaxIn, int shellMaxIn, int NparticlesIn){
  this->density = densityIn;
  this->tzMax = tzMaxIn;
  this->shellMax = shellMaxIn;
  this->Nparticles = NparticlesIn;
  if(this->tzMax == 1){
    this->chanTzMax = -2;
    this->modChanTzMax = 0;
  } else if(this->tzMax == 2){
    this->chanTzMax = 2;
    this->modChanTzMax = 2;
  } else {
    this->chanTzMax = 8000;
    this->modChanTzMax = 8000;
    std::cout << "wrong tzMax" << std::endl;
  }

  this->EMax = 0;
  for(int i = 0; i < this->shellMax; this->EMax++){
    if(isASumOfThreeSquares(this->EMax)){
      i++;
    }
  }

  this->nMax = smallestSquareRootAtLeast(this->EMax);

  this->Nspstates = 0;
  int nMaxActual = 0;
  for( int nx = -this->nMax; nx <= this->nMax; nx++){
    for( int ny = -this->nMax; ny <= this->nMax; ny++){
      for( int nz = -this->nMax; nz <= this->nMax; nz++){
        if(nx*nx + ny*ny + nz*nz <= this->EMax){
          if(nx > nMaxActual){
            nMaxActual = nx;
          }
          for( int sz = -1; sz <= 1; sz = sz+2){
            for( int tz = -1; tz <= this->tzMax - 1; tz = tz+2){
              this->Nspstates++;
            } // end tz loop
          } // end sz loop
        }
      } // end nz loop
    } // end ny loop
  } // end nx loop

  this->nMax = nMaxActual;
}


void infMatterSPBasis::generateIndexMap(){

  this->indexMap = new int* [this->Nspstates];
  for(int i = 0; i < this->Nspstates; i++){
    this->indexMap[i] = new int[5];
  }

  int *shells = new int [this->EMax + 2];
  int E;

  for( E = 0; E <= this->EMax + 1; E++){
    shells[E] = 0;
  }

  //Determine how many single particle states are in each shell.
  for( int nx = -this->nMax; nx <= this->nMax; nx++){
    for( int ny = -this->nMax; ny <= this->nMax; ny++){
      for( int nz = -this->nMax; nz <= this->nMax; nz++){
        E = nx*nx + ny*ny + nz*nz;
        if( E <= this->EMax){
          for( int sz = -1; sz <= 1; sz = sz+2){
            for( int tz = -1; tz <= this->tzMax - 1; tz = tz+2){
              shells[E + 1]++;
            } // end tz loop
          } // end sz loop
        } // end if
      } // end nz loop
    } // end ny loop
  } // end nx loop

  //Determine offsets into index array for each energy level
  for( E = 1; E <= this->EMax; E++){
    shells[E] += shells[E - 1];
  }

  //Construct the index map
  for( int nx = -this->nMax; nx <= this->nMax; nx++){
    for( int ny = -this->nMax; ny <= this->nMax; ny++){
      for( int nz = -this->nMax; nz <= this->nMax; nz++){
        E = nx*nx + ny*ny + nz*nz;
        if( E <= this->EMax){
          for( int sz = -1; sz <= 1; sz = sz+2){
            for( int tz = -1; tz <= this->tzMax - 1; tz = tz+2){
              this->indexMap[shells[E]][0] = nx;
              this->indexMap[shells[E]][1] = ny;
              this->indexMap[shells[E]][2] = nz;
              this->indexMap[shells[E]][3] = sz;
              this->indexMap[shells[E]][4] = tz;
              shells[E]++;
            } // end tz loop
          } // end sz loop
        } // end if
      } // end nz loop
    } // end ny loop
  } // end nx loop
  delete[] shells;
} // end generateIndexMap


// requires generate index map first
void infMatterSPBasis::generateBasis(){
  this->spEnergy = new double[this->Nspstates];
  double hbarc = 197.3269788; // MeVfm
  //double m_neutronc2 = 939.565378; // MeV // accurate value
  double m_neutronc2 = 939.565; // MeV Gaute's value
  //double m_protonc2 = 938.272046; // MeV
  //double m_symmc2 = (m_neutronc2 + m_protonc2)*0.5;
  double massc2;
  if( tzMax == 1){
    massc2 = m_neutronc2;
  } else if( tzMax == 2){
    //massc2 = m_symmc2;
    massc2 = m_neutronc2;
  } else {
    std::cout << "tzMax = 1 or 2 ONLY" << std::endl;
  }

  double prefactor = hbarc*hbarc/(2.*massc2);
  double L = pow(Nparticles/density, 1./3.);
  std::cout << "mass*c^2: " << massc2 << "MeV" << std::endl;
  std::cout << "hbar^2/(2m) = " << prefactor << "MeV*fm^2" << std::endl;
  std::cout << "L = " << L << "fm V= " << L*L*L << "fm^3" << std::endl;
  std::cout << "k_fermi " << std::setprecision(16) << pow((6.*M_PI*M_PI*density/(2.*tzMax)),1./3.) << "fm^-1" << std::endl;

   // vector<state> psi; // vector of sp states
  for(int i = 0; i < this->Nspstates; i++){
    this->spEnergy[i] = (prefactor*4*M_PI*M_PI/(L*L) ) *
                        (this->indexMap[i][0]*this->indexMap[i][0] +
                         this->indexMap[i][1]*this->indexMap[i][1] +
                         this->indexMap[i][2]*this->indexMap[i][2]);
  }
} // end generateBasis


int infMatterSPBasis::checkSympqrs(int p, int q, int r, int s){
  int result = 0;
  //int ** indexMap = SPbasis.indexMap;
  // check spin and isospin
  if( (this->indexMap[p][3] + this->indexMap[q][3] ==
        this->indexMap[r][3] + this->indexMap[s][3]) &&
      (this->indexMap[p][4] + this->indexMap[q][4] ==
       this->indexMap[r][4] + this->indexMap[s][4]) ){
    // check momentum conservation
    if( (this->indexMap[p][0] + this->indexMap[q][0] ==
          this->indexMap[r][0] + this->indexMap[s][0]) &&
        (this->indexMap[p][1] + this->indexMap[q][1] ==
         this->indexMap[r][1] + this->indexMap[s][1]) &&
        (this->indexMap[p][2] + this->indexMap[q][2] ==
         this->indexMap[r][2] + this->indexMap[s][2] ) ) {
      result = 1;
    }
  }
  return result;
} // end checkSympqrs

int infMatterSPBasis::checkModSympqrs(int p, int q, int r, int s){
  int result = 0;
  // check spin and isospin
  if( (this->indexMap[p][3] - this->indexMap[q][3] ==
        this->indexMap[r][3] - this->indexMap[s][3]) &&
      (this->indexMap[p][4] - this->indexMap[q][4] ==
       this->indexMap[r][4] - this->indexMap[s][4]) ){
    // check momentum conservation
    if( (this->indexMap[p][0] - this->indexMap[q][0] ==
          this->indexMap[r][0] - this->indexMap[s][0]) &&
        (this->indexMap[p][1] - this->indexMap[q][1] ==
         this->indexMap[r][1] - this->indexMap[s][1]) &&
        (this->indexMap[p][2] - this->indexMap[q][2] ==
         this->indexMap[r][2] - this->indexMap[s][2] ) ) {
      result = 1;
    }
  }
  return result;
} // end checkSympqrs

int infMatterSPBasis::checkChanSym(int p, int q, int ichan){
  int result = 0;

  if( (this->indexMap[p][0] + this->indexMap[q][0] ==
        this->chanValue[ichan].chanNx) &&
      (this->indexMap[p][1] + this->indexMap[q][1] ==
       this->chanValue[ichan].chanNy) &&
      (this->indexMap[p][2] + this->indexMap[q][2] ==
       this->chanValue[ichan].chanNz) &&
      (this->indexMap[p][3] + this->indexMap[q][3] ==
       this->chanValue[ichan].chanSz) &&
      (this->indexMap[p][4] + this->indexMap[q][4] ==
       this->chanValue[ichan].chanTz) ){
    result = 1;
  }
  return result;
} // end checkChanSym

int infMatterSPBasis::checkChanSym2(int p, int q, int chanNx, int chanNy, int chanNz, int chanSz, int chanTz){
  int result = 0;

  if( (this->indexMap[p][0] + this->indexMap[q][0] == chanNx) &&
      (this->indexMap[p][1] + this->indexMap[q][1] == chanNy) &&
      (this->indexMap[p][2] + this->indexMap[q][2] == chanNz) &&
      (this->indexMap[p][3] + this->indexMap[q][3] == chanSz) &&
      (this->indexMap[p][4] + this->indexMap[q][4] == chanTz) ){

    result = 1;
  }
  return result;
} // end checkChanSym

int infMatterSPBasis::checkChanModSym2(int p, int q, int chanNx, int chanNy, int chanNz, int chanSz, int chanTz){
  int result = 0;

  if( (this->indexMap[p][0] - this->indexMap[q][0] == chanNx) &&
      (this->indexMap[p][1] - this->indexMap[q][1] == chanNy) &&
      (this->indexMap[p][2] - this->indexMap[q][2] == chanNz) &&
      (this->indexMap[p][3] - this->indexMap[q][3] == chanSz) &&
      (this->indexMap[p][4] - this->indexMap[q][4] == chanTz) ){

    result = 1;
  }
  return result;
} // end checkChanSym


int infMatterSPBasis::checkChanModSym(int p, int q, int ichan){
  int result = 0;

  if( (this->indexMap[p][0] - this->indexMap[q][0] ==
        this->chanModValue[ichan].chanNx) &&
      (this->indexMap[p][1] - this->indexMap[q][1] ==
       this->chanModValue[ichan].chanNy) &&
      (this->indexMap[p][2] - this->indexMap[q][2] ==
       this->chanModValue[ichan].chanNz) &&
      (this->indexMap[p][3] - this->indexMap[q][3] ==
       this->chanModValue[ichan].chanSz) &&
      (this->indexMap[p][4] - this->indexMap[q][4] ==
       this->chanModValue[ichan].chanTz) ){
    result = 1;
  }
  return result;
} // end checkChanModSym

void infMatterSPBasis::setUpTwoStateChannels_andCalcDims(){
}
// special for mbpt2 
void infMatterSPBasis::setUpTwoStateChannels(){
  int channelNmax = 2*this->nMax;
  int channelTzMax = this->chanTzMax;  
  int channelSzMax = 2;



  // first count how many unique channels there are
  // could potentially nuke these loops and just overallocate
  // upper bound is known
  //int NumChannels = 0;
  //int iMerge = 0;
  int NumChannels = (2*channelNmax+1)*(2*channelNmax+1)*(2*channelNmax+1)*3*(2+channelTzMax/2);
  printf("NumChannels: %d\n",NumChannels);


  //int nonZeroChan;
  this->Nchannels = NumChannels;
  this->chanValue = new channelBundle[NumChannels];
  this->threeBodyChanValue = new threeBodyChannelBundle[this->Nspstates];

  int channelCount = 0;
  for(int ichanNx = -channelNmax; ichanNx <= channelNmax; ichanNx++){
    for(int ichanNy = -channelNmax; ichanNy <= channelNmax; ichanNy++){
      for(int ichanNz = -channelNmax; ichanNz <= channelNmax; ichanNz++){
        for(int ichanSz = -channelSzMax; ichanSz <= channelSzMax; ichanSz = ichanSz +2){
          for(int ichanTz = -2; ichanTz <= channelTzMax; ichanTz = ichanTz +2){
            this->chanValue[channelCount].chanNx = ichanNx;
            this->chanValue[channelCount].chanNy = ichanNy;
            this->chanValue[channelCount].chanNz = ichanNz;
            this->chanValue[channelCount].chanSz = ichanSz;
            this->chanValue[channelCount].chanTz = ichanTz;
            this->chanValue[channelCount].ppDim = 0;
            this->chanValue[channelCount].hhDim = 0;
            this->chanValue[channelCount].hpDim = 0;
            this->chanValue[channelCount].phDim = 0;

            channelCount++;
          }
        }
      }
    }
  }



  int temp_chanNx, temp_chanNy, temp_chanNz, temp_chanSz, temp_chanTz;
  int iTB;
  // set up dimension
  for(int p = 0; p < this->Nspstates; p++){
    for(int q = 0; q < this->Nspstates; q++){
      temp_chanNx = this->indexMap[p][0] + this->indexMap[q][0];
      temp_chanNy = this->indexMap[p][1] + this->indexMap[q][1];
      temp_chanNz = this->indexMap[p][2] + this->indexMap[q][2];
      temp_chanSz = this->indexMap[p][3] + this->indexMap[q][3];
      temp_chanTz = this->indexMap[p][4] + this->indexMap[q][4];
      iTB = TBchanIndexFunction(temp_chanNx, temp_chanNy, temp_chanNz, temp_chanSz, temp_chanTz);

      if( p < this->Nparticles ) {
        if( q < this->Nparticles ) {
          this->chanValue[iTB].hhDim++;
        } else {
          this->chanValue[iTB].hpDim++;
        }
      } else {
        if( q < this->Nparticles ){
          this->chanValue[iTB].phDim++;
        } else {
          this->chanValue[iTB].ppDim++;
        }
      }
    }
  }



  // allocate maps
  for(int ichan=0; ichan < NumChannels; ichan++){
    // printf("ichan: %d, ipp: %d, iph: %d, ihp: %d , ihh: %d\n",ichan,this->chanValue[ichan].ppDim,this->chanValue[ichan].phDim,this->chanValue[ichan].hpDim,this->chanValue[ichan].hhDim);


    this->chanValue[ichan].ppMap = new int*[this->chanValue[ichan].ppDim];
    for(int i=0; i<this->chanValue[ichan].ppDim; i++){
      this->chanValue[ichan].ppMap[i] = new int[2];
    }

    this->chanValue[ichan].hhMap = new int*[this->chanValue[ichan].hhDim];
    for(int i=0; i<this->chanValue[ichan].hhDim; i++){
      this->chanValue[ichan].hhMap[i] = new int[2];
    }

    this->chanValue[ichan].hpMap = new int*[this->chanValue[ichan].hpDim];
    for(int i=0; i<this->chanValue[ichan].hpDim; i++){
      this->chanValue[ichan].hpMap[i] = new int[2];
    }

    this->chanValue[ichan].phMap = new int*[this->chanValue[ichan].phDim];
    for(int i=0; i<this->chanValue[ichan].hpDim; i++){
      this->chanValue[ichan].phMap[i] = new int[2];
    }
    this->chanValue[ichan].ppDim = 0;
    this->chanValue[ichan].hhDim = 0;
    this->chanValue[ichan].hpDim = 0;
    this->chanValue[ichan].phDim = 0;

  }



  // set up maps
  for(int p = 0; p < this->Nspstates; p++){
    for(int q = 0; q < this->Nspstates; q++){
      temp_chanNx = this->indexMap[p][0] + this->indexMap[q][0];
      temp_chanNy = this->indexMap[p][1] + this->indexMap[q][1];
      temp_chanNz = this->indexMap[p][2] + this->indexMap[q][2];
      temp_chanSz = this->indexMap[p][3] + this->indexMap[q][3];
      temp_chanTz = this->indexMap[p][4] + this->indexMap[q][4];
      iTB = TBchanIndexFunction(temp_chanNx, temp_chanNy, temp_chanNz, temp_chanSz, temp_chanTz);
      // printf("ihh: %d, ihp: %d, iph: %d , ipp: %d\n",ihh,ihp,iph,ipp);
      //printf("p: %d, q: %d\n",p,q);
      if( p < this->Nparticles ) {
        if( q < this->Nparticles ) {
          //printf("iTB: %d, ihhD: %d, ihpD: %d, iphD: %d , ippD: %d\n",iTB,this->chanValue[iTB].hhDim, this->chanValue[iTB].hpDim, this->chanValue[iTB].phDim ,this->chanValue[iTB].ppDim);
          this->chanValue[iTB].hhMap[this->chanValue[iTB].hhDim][0] = p;
          this->chanValue[iTB].hhMap[this->chanValue[iTB].hhDim][1] = q;
          this->chanValue[iTB].hhDim++;
        } else {
          this->chanValue[iTB].hpMap[this->chanValue[iTB].hpDim][0] = p;
          this->chanValue[iTB].hpMap[this->chanValue[iTB].hpDim][1] = q;
          this->chanValue[iTB].hpDim++;
        }
      } else {
        if( q < this->Nparticles ){
          this->chanValue[iTB].phMap[this->chanValue[iTB].phDim][0] = p;
          this->chanValue[iTB].phMap[this->chanValue[iTB].phDim][1] = q;
          this->chanValue[iTB].phDim++;
        } else {
          this->chanValue[iTB].ppMap[this->chanValue[iTB].ppDim][0] = p;
          this->chanValue[iTB].ppMap[this->chanValue[iTB].ppDim][1] = q;
          this->chanValue[iTB].ppDim++;
        }
      }

    }
  } 
} // end setUpTwoStateChannels

void infMatterSPBasis::printBasis(){

  std::cout << "SPBasis:" << std::endl;
  std::cout << "i nx ny nz sz tz E" << std::endl;
  for( int p = 0; p < this->Nspstates; p++ ) {
    std::cout << p << " " << this->indexMap[p][0] << " " << this->indexMap[p][1] << " " << this->indexMap[p][2] << " " << this->indexMap[p][3] << " " << this->indexMap[p][4] << " " << this->spEnergy[p] << std::endl;
  }
} // end printBasis

void infMatterSPBasis::deallocate(){

  for(int i = 0; i < this->Nspstates; i++){
    delete [] this->indexMap[i];
  }
  delete [] this->indexMap;
  delete [] this->spEnergy;
  delete [] this->chanValue;
  delete [] this->chanModValue;
} // end deallocate



double infMatterSPBasis::calc_TBME(int p, int q, int r, int s){
  double vout = 0.;
  int *qi = this->indexMap[p];
  int *qj = this->indexMap[q];
  int *qk = this->indexMap[r];
  int *ql = this->indexMap[s];

  int initialMomx, initialMomy, initialMomz;
  int finalMomx, finalMomy, finalMomz;
  // Momentum Conservation Checks.
  initialMomx = qi[0] + qj[0];
  initialMomy = qi[1] + qj[1];
  initialMomz = qi[2] + qj[2];
  finalMomx = qk[0] + ql[0];
  finalMomy = qk[1] + ql[1];
  finalMomz = qk[2] + ql[2];

  // maybe add a conservation of spin if here.
  if( initialMomx == finalMomx && initialMomy == finalMomy && initialMomz == finalMomz ){
    double L = pow(this->Nparticles/this->density,1./3.);
    double V_R, V_T, V_S;
    double V_0R, V_0T, V_0S;
    double kappa_R, kappa_T, kappa_S;
    double *ki = new double[3];
    double *kj = new double[3];
    double *kk = new double[3];
    double *kl = new double[3];
    double *relMomBra = new double[3];
    double *relMomKet = new double[3];
    double *relMomTransf = new double[3];
    double qSquared, spinEx, isoSpinEx;
    double IsIt, PsIt, PsPt, IsPt;
    V_0R = 200; //MeV
    V_0T = 178; //MeV
    V_0S = 91.85; //MeV
    kappa_R = 1.487; //fm^-2
    kappa_T = 0.639; //fm^-2
    kappa_S = 0.465; //fm^-2

    qSquared = 0.;
    for( int i = 0; i < 3; i++){
      ki[i] = 2*M_PI*qi[i]/L;
      kj[i] = 2*M_PI*qj[i]/L;
      kk[i] = 2*M_PI*qk[i]/L;
      kl[i] = 2*M_PI*ql[i]/L;
      relMomBra[i] = 0.5*(ki[i] - kj[i]);
      relMomKet[i] = 0.5*(kk[i] - kl[i]);
      relMomTransf[i] = relMomBra[i] - relMomKet[i];
      qSquared += relMomTransf[i]*relMomTransf[i];
    }

    V_R = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5)*exp(-qSquared/(4*kappa_R));
    V_T = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5)*exp(-qSquared/(4*kappa_T));
    V_S = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5)*exp(-qSquared/(4*kappa_S));

    spinEx = spinExchangeMtxEle(qi[3],qj[3],qk[3],ql[3]);
    isoSpinEx = spinExchangeMtxEle(qi[4],qj[4],qk[4],ql[4]);

    // 4 terms, IsIt, PsIt, PsPt, IsPt
    // identity spin, identity isospin
    IsIt = kron_del(qi[3],qk[3])*kron_del(qj[3],ql[3])*kron_del(qi[4],qk[4])*kron_del(qj[4],ql[4]);
    // Exchange spin, identity isospin
    PsIt = spinEx*kron_del(qi[4],qk[4])*kron_del(qj[4],ql[4]);
    // Exchange spin, Exchange isospin
    PsPt = spinEx*isoSpinEx;
    // identity spin, Exchange isospin
    IsPt = kron_del(qi[3],qk[3])*kron_del(qj[3],ql[3])*isoSpinEx;

    vout = 0.5*(V_R + 0.5*V_T + 0.5*V_S)*IsIt
      + 0.25*(V_T - V_S)*PsIt
      - 0.5*(V_R + 0.5*V_T + 0.5*V_S)*PsPt
      - 0.25*(V_T - V_S)*IsPt;

    delete[] ki;
    delete[] kj;
    delete[] kk;
    delete[] kl;
    delete[] relMomBra;
    delete[] relMomKet;
    delete[] relMomTransf;
  } // end Momentum Conversation Check

  return vout;
}


int infMatterSPBasis::spinExchangeMtxEle(int i, int j, int k, int l){
  if( i == l && j == k ){
    return 1;
  } else {
    return 0;
  }
} // end spinEx



int infMatterSPBasis::kron_del(int i, int j){
  if(i != j){
    return 0;
  }
  return 1;
} // end kron_del


int infMatterSPBasis::TBchanIndexFunction(int chanNx, int chanNy, int chanNz, int chanSz, int chanTz){
  int index;
  int TBnMax = 4*this->nMax + 1;
  int TBsMax = 3;
  int TBtMax = 2+this->chanTzMax/2;// 1 or 3 pnm vs snm

  chanNx += TBnMax/2;
  chanNy += TBnMax/2;
  chanNz += TBnMax/2;
  chanSz = (chanSz+2)/2;
  chanTz = (chanTz+2)/2;  // 0 or 0,1,2
  index = chanNx*TBnMax*TBnMax*TBsMax*TBtMax
    + chanNy*TBnMax*TBsMax*TBtMax
    + chanNz*TBsMax*TBtMax
    + chanSz*TBtMax
    + chanTz;
  return index;
} // end indexMap

int infMatterSPBasis::TBmodChanIndexFunction(int modChanNx, int modChanNy, int modChanNz, int modChanSz, int modChanTz){
  int index;
  int TBnMax = 4*this->nMax + 1;
  int TBsMax = 3;
  int TBtMax = 1 + this->modChanTzMax; //

  modChanNx += TBnMax/2;
  modChanNy += TBnMax/2;
  modChanNz += TBnMax/2;
  modChanSz = (modChanSz+2)/2;
  modChanTz = (modChanTz+this->modChanTzMax)/2; // 0 or 0,1,2
  index = modChanNx*TBnMax*TBnMax*TBsMax*TBtMax
    + modChanNy*TBnMax*TBsMax*TBtMax
    + modChanNz*TBsMax*TBtMax
    + modChanSz*TBtMax
    + modChanTz;
  return index;
}

