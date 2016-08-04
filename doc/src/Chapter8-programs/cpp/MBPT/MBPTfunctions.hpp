#ifndef CCFUNCTIONS_H
#define CCFUNCTIONS_H

#include <iostream>
#include <math.h>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <omp.h>

//LAPACK functions
extern "C" void dgemm_(char* ta,char* tb,int* m,int* n,int* k,double* al,double* a,int* la,double* b,int* lb,double* be,double* c,int* lc);
#define RM_dgemm(A, B, C, m, n, k, alpha, beta, transf_A, transf_B) dgemm_(transf_B, transf_A, n, m, k, alpha, B, n, A, k, beta, C, n)
#define RMT_dgemm(A, B, C, m, n, k, alpha, beta, transf_A, transf_B) dgemm_(transf_B, transf_A, n, m, k, alpha, B, k, A, k, beta, C, n)

#define hbarc_MeVfm 197.3269788 // MeV fm
#define m_neutronc2 939.5654133 // MeV
#define m_protonc2 938.2720813 // MeV

const std::string PATH = "inputs/";

struct Input_Parameters;
struct State;
struct Model_Space;
struct Channels;

void MBPT2_0(const Input_Parameters &Parameters, const Model_Space &Space);
void MBPT2_1(const Input_Parameters &Parameters, const Model_Space &Space);
void MBPT2_2(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan);
void MBPT2_3(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan);
void MBPT2_4(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan);
void MBPT3_0(const Input_Parameters &Parameters, const Model_Space &Space);
void MBPT3_1(const Input_Parameters &Parameters, const Model_Space &Space);
void MBPT3_2(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan);
void MBPT3_3(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan);
void MBPT3_4(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan);
void plus(State &S, const State &S1, const State &S2);
void Setup_Channels_MBPT(const Input_Parameters &Parameters, const Model_Space &Space, Channels &Chan);
void Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space);
int Chan_2bInd(const Model_Space &Space, const State &State);
int kron_del(const int &i, const int &j);
int spinExchangeMtxEle(const int &i, const int &j, const int &k, const int &l);
double vint_Minnesota_Momentum(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L);

struct Input_Parameters{
  int Nshells; //number of neutrons shells
  int Pshells; //number of protons shells
  int N; //number of neutrons
  int P; //number of protons
  double density;
  int Shells; // total number of shells
  int MBPT_Approx; // 2 for MBPT2, 3 for MBPT3_pp
  int MBPT_Function; // 0 for serial, 1 for parallel, 2 for block serial, 3 for block parallel, 4 for block M-M
};

struct State{
  int t; // x2
  int m; // x2
  int nx;
  int ny;
  int nz;
  double energy;
  std::string type;
};

struct Model_Space{
  int indp; //number of proton orbits
  int indn; //number of neutron orbits
  int indpar; //number of particle orbits
  int indhol; //number of hole orbits
  int indtot; //number of total orbits

  State *qnums;
  struct State qmins;
  struct State qsizes;

  int Nmax;
  int nmax;
  int *map_2b;
  int size_2b;
};

struct Channels{
  int size;

  int* nhh;
  int* npp;
  int** hhvec;
  int** ppvec;
};

#endif
