#!/bin/bash

SRCDIR=$PWD
BUILDDIR=$PWD
SRCDIR2=$HOME/ManyBodyPhysics/LectureNotesPhysics/doc/src/Chapter5-programs/src
BUILDDIR2=$HOME/ManyBodyPhysics/LectureNotesPhysics/doc/src/Chapter5-programs/include/obj

CXX=clang++
AR=ar

LDFLAGS="-L/usr/local//lib"
CXXFLAGS="-O2 -w -DMPICH_IGNORE_CXX_SEEK"
ARFLAGS="ruv"

INCLUDE_FLAGS="-I$SRCDIR/include -I/usr/local/include"

if [ -d $BUILDDIR2 ]; then
  rm -r $BUILDDIR2
fi

mkdir $BUILDDIR2
cd $BUILDDIR2

$CXX $CXXFLAGS $LDFLAGS $INCLUDE_FLAGS -c $SRCDIR2/*.C
$AR $ARFLAGS $BUILDDIR/nrCPS.a $BUILDDIR2/*.o
ranlib  $BUILDDIR/nrCPS.a
