#!/bin/bash

SRCDIR=$PWD
BUILDDIR=$PWD
SRCDIR2=$HOME/Desktop/UnitaryFermions/nrCPS_2body/src
BUILDDIR2=$HOME/Desktop/UnitaryFermions/nrCPS_2body/include/objs
FFTWDIR=$HOME/fftw-3.3.4

CXX=g++
AR=ar

LDFLAGS="-L$FFTWDIR/lib"
CXXFLAGS="-O2 -w -DMPICH_IGNORE_CXX_SEEK"
ARFLAGS="ruv"

INCLUDE_FLAGS="-I$SRCDIR/include -I$FFTWDIR/include"

if [ -d $BUILDDIR2 ]; then
  rm -r $BUILDDIR2
fi

mkdir $BUILDDIR2
cd $BUILDDIR2

$CXX $CXXFLAGS $LDFLAGS $INCLUDE_FLAGS -c $SRCDIR2/*.C
$AR $ARFLAGS $BUILDDIR/nrCPS.a $BUILDDIR2/*.o
ranlib  $BUILDDIR/nrCPS.a
