#!/bin/bash

SRCDIR=../
BUILDDIR=../
FFTWDIR=$HOME/libraries/fftw-3.2.2

CXX=clang++

LDFLAGS="-L/usr/local//lib"
CXXFLAGS="-O2 -Wall -DMPICH_IGNORE_CXX_SEEK"

INCLUDE_FLAGS="-I$SRCDIR/include -I/usr/local/include"

$CXX $CXXFLAGS main.C $BUILDDIR/nrCPS.a $LDFLAGS $INCLUDE_FLAGS -lfftw3
