// ranlxs.C
// An implementation of Luscher's ranlxs code.

#include "uniform_deviates.h"
#include "error.h"
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define BASE 0x1000000
#define MASK 0xffffff

#define STEP(pi,pj) \
        d=(*pj).c1.c1-(*pi).c1.c1-carry.c1; \
        (*pi).c2.c1+=(d<0); \
        d+=BASE; \
        (*pi).c1.c1=d&MASK; \
        d=(*pj).c1.c2-(*pi).c1.c2-carry.c2; \
        (*pi).c2.c2+=(d<0); \
        d+=BASE; \
        (*pi).c1.c2=d&MASK; \
        d=(*pj).c1.c3-(*pi).c1.c3-carry.c3; \
        (*pi).c2.c3+=(d<0); \
        d+=BASE; \
        (*pi).c1.c3=d&MASK; \
        d=(*pj).c1.c4-(*pi).c1.c4-carry.c4; \
        (*pi).c2.c4+=(d<0); \
        d+=BASE; \
        (*pi).c1.c4=d&MASK; \
        d=(*pj).c2.c1-(*pi).c2.c1; \
        carry.c1=(d<0); \
        d+=BASE; \
        (*pi).c2.c1=d&MASK; \
        d=(*pj).c2.c2-(*pi).c2.c2; \
        carry.c2=(d<0); \
        d+=BASE; \
        (*pi).c2.c2=d&MASK; \
        d=(*pj).c2.c3-(*pi).c2.c3; \
        carry.c3=(d<0); \
        d+=BASE; \
        (*pi).c2.c3=d&MASK; \
        d=(*pj).c2.c4-(*pi).c2.c4; \
        carry.c4=(d<0); \
        d+=BASE; \
        (*pi).c2.c4=d&MASK

Ranlxs::Ranlxs(int level)
{
//  init=0;
  LuxLevel(level);
  Init(1); 
}

Ranlxs::~Ranlxs()
{
}

void Ranlxs::Init(long seed)
{

  char* fname = "Ranlxs::Init(long)";
  int i,k,l;
  int ibit,jbit,xbit[31];
  int ix,iy;

  if ((INT_MAX<2147483647)||(FLT_RADIX!=2)||(FLT_MANT_DIG<24)) {
    ERR.General(fname,"Arithmetic on this machine is not suitable for generator.");
  }

  DefineConstants();

  i=seed;

  for (k=0;k<31;k++) {
    xbit[k]=i%2;
    i/=2;
  }

  if ((seed<=0)||(i!=0)) { ERR.General(fname,"Bad choice of seed (should be between 1 and 2^31-1)."); }

  ibit=0;
  jbit=18;

  for (i=0;i<4;i++) {
    for (k=0;k<24;k++) {
      ix=0;
      for (l=0;l<24;l++) {
        iy=xbit[ibit];
        ix=2*ix+iy;
        xbit[ibit]=(xbit[ibit]+xbit[jbit])%2;
        ibit=(ibit+1)%31;
        jbit=(jbit+1)%31;
      }
      if ((k%4)==i) { ix=16777215-ix; }
      x.num[4*k+i]=ix;
    }
 }

 carry.c1=0;
 carry.c2=0;
 carry.c3=0;
 carry.c4=0;

 ir=0;
 jr=7;
 is=95;
 is_old=0;
 prm=pr%12;
// init=1;
}

double Ranlxs::Run()
{
//  if (init==0) {
//    LuxLevel(0);
//    Init(1);
//  }

  is=next[is];
  if (is==is_old)
    Update();
  return one_bit*(float)(x.num[is]);      

}

void Ranlxs::LuxLevel(int level)
{
  char* fname = "Ranlxs::LuxLevel(int)";

  if (level==0)
    pr=109;
  else if (level==1)
    pr=202;
  else if (level==2)
    pr=397;
  else
    ERR.General(fname,"Bad choice of luxury level (should be 0,1 or 2).");
}

void Ranlxs::Update(void)
{
  int k,kmax,d;
  dble_vec_t *pmin,*pmax,*pi,*pj;

  kmax=pr;
  pmin=&x.vec[0];
  pmax=pmin+12;
  pi=&x.vec[ir];
  pj=&x.vec[jr];

  for (k=0;k<kmax;k++) {
    STEP(pi,pj);
    pi+=1;
    pj+=1;
    if (pi==pmax) { pi=pmin; }
    if (pj==pmax) { pj=pmin; }
  }

  ir+=prm;
  jr+=prm;
  if (ir>=12) { ir-=12; }
  if (jr>=12) { jr-=12; }
  is=8*ir;
  is_old=is;
}

void Ranlxs::DefineConstants(void)
{
  int k;

  one_bit=(float)(ldexp(1.0,-24));

  for (k=0;k<96;k++) { next[k]=(k+1)%96; }
}

int Ranlxs::StateSize()
{
  return(105);
}

void Ranlxs::GetState(long* state)
{
  char* fname = "Ranlxs::GetState(long*)";
  int k;  

//  if (init==0)
//    ERR.General(fname,"Undefined state (generator not initialized).");

  state[0]=StateSize();

  for (k=0;k<96;k++)
    state[k+1]=x.num[k];

  state[97]=carry.c1;
  state[98]=carry.c2;
  state[99]=carry.c3;
  state[100]=carry.c4;

  state[101]=pr;
  state[102]=ir;
  state[103]=jr;
  state[104]=is;
}

void Ranlxs::SetState(long* state)
{
  char* fname = "Ranlxs::SetState(long*)";
  int k;

  if ((INT_MAX<2147483647)||(FLT_RADIX!=2)||(FLT_MANT_DIG<24))
    ERR.General(fname,"Arithmetic on this machine is not suitable for generator.");

  DefineConstants();

  if (state[0]!=StateSize())
    ERR.General(fname,"Unexpected input data.");

  for (k=0;k<96;k++) {
    if ((state[k+1]<0)||(state[k+1]>=167777216))
      ERR.General(fname,"Unexpected input data.");

    x.num[k]=state[k+1];
  }

  if (((state[97]!=0)&&(state[97]!=1))||
     ((state[98]!=0)&&(state[98]!=1))||
     ((state[99]!=0)&&(state[99]!=1))||
     ((state[100]!=0)&&(state[100]!=1)))
    ERR.General(fname,"Unexpected input data.");

  carry.c1=state[97];
  carry.c2=state[98];
  carry.c3=state[99];
  carry.c4=state[100];

  pr=state[101];
  ir=state[102];
  jr=state[103];
  is=state[104];
  is_old=8*ir;
  prm=pr%12;
//  init=1;

  if (((pr!=109)&&(pr!=202)&&(pr!=397))||
     (ir<0)||(ir>11)||(jr<0)||(jr>11)||(jr!=((ir+7)%12))||
     (is<0)||(is>95))
    ERR.General(fname,"Unexpected input data.");
}

#undef BASE
#undef MASK
