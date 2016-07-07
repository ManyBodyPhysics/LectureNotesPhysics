#!/bin/bash

SRCDIR=../
BUILDDIR=../
FFTWDIR=$HOME/libraries/fftw-3.2.2

CXX=g++
#/home/rdm/sfw/mpich2/mpich2-1.0.7-install/bin/mpicxx

LDFLAGS="-L$FFTWDIR/lib"
CXXFLAGS="-O2 -Wall -DMPICH_IGNORE_CXX_SEEK"
ARFLAGS=

INCLUDE_FLAGS="-I$SRCDIR/include -I$FFTWDIR/include"

$CXX $CXXFLAGS main.C $BUILDDIR/nrCPS.a $LDFLAGS $INCLUDE_FLAGS -lfftw3
