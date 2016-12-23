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
//  MtK-Utilities.cpp  --  utilities stuff.
//


#include <cstdlib>
#include <iostream>
using namespace std;

//#include "MtK-Utilities.hh"
#include "MtK-Global_data.h"

//
// sorting routines:
//

void Sort_u_double2dim(int ntot, int NDIM, int isort, double dbl[]) {
	int    i,ir,j,l;
  double *dblex = new double[NDIM];
	int n;
	
	if (ntot < 2) return;
	l = ntot/2;
	ir = ntot-1;
	while(1) {
		if (l > 0) {
			l--;
			for(n=0; n<NDIM; ++n) dblex[n] = dbl[l*NDIM+n];
		} else {
			for(n=0; n<NDIM; ++n) dblex[n]         = dbl[ir*NDIM+n];
			for(n=0; n<NDIM; ++n) dbl[ir*NDIM+n]   =         dbl[n];
			if ((--ir) == 0) {for(n=0; n<NDIM; ++n) dbl[n] = dblex[n];
				break; }
		}
		i = l;
		j = l+l+1;
		while (j <= ir) {
			if ((j < ir) && (dbl[j*NDIM+isort] < dbl[(j+1)*NDIM+isort])) j++;
			if (dblex[isort]  < dbl[j*NDIM+isort]) {
				for(n=0; n<NDIM; ++n) dbl[i*NDIM+n] = dbl[j*NDIM+n];
				i = j;
				j = j+j+1;
			} else break;
		}
		for(n=0; n<NDIM; ++n) dbl[i*NDIM+n] = dblex[n];
    }
	
	delete [] dblex;  dblex = NULL;
	return;
}


void Sort_d_double2dim(int ntot, int NDIM, int isort, double dbl[]) {
	int    i,ir,j,l;
  double *dblex = new double[NDIM];
	int n;
	
	if (ntot < 2) return;
	l = ntot/2;
	ir = ntot-1;
	while(1) {
		if (l > 0) {
			l--;
			for(n=0; n<NDIM; ++n) dblex[n] = dbl[l*NDIM+n];
		} else {
			for(n=0; n<NDIM; ++n) dblex[n]         = dbl[ir*NDIM+n];
			for(n=0; n<NDIM; ++n) dbl[ir*NDIM+n]   =         dbl[n];
			if ((--ir) == 0) {for(n=0; n<NDIM; ++n) dbl[n] = dblex[n];
				break; }
		}
		i = l;
		j = l+l+1;
		while (j <= ir) {
			if ((j < ir) && (dbl[j*NDIM+isort] > dbl[(j+1)*NDIM+isort])) j++;
			if (dblex[isort]  > dbl[j*NDIM+isort]) {
				for(n=0; n<NDIM; ++n) dbl[i*NDIM+n] = dbl[j*NDIM+n];
				i = j;
				j = j+j+1;
			} else break;
		}
		for(n=0; n<NDIM; ++n) dbl[i*NDIM+n] = dblex[n];
    }
	
	delete [] dblex;  dblex = NULL;
	return;
}

