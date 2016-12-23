 

 “MattersK code for ADC(3) calculations of infinite nucleonic matter”


 The code ‘MattersK’ found in this chapter allows you to calculate the
one-body propagator (a.k.a. Green’s function) for infinite symmetric nuclear
matter and pure neutron matter. It uses a cartesian basis in momentum 
space, which discretized using boxes with standard periodic boundary
conditions (PBC).  Refer to the LNP chapter mentioned below here for 
more details on the functioning of this software.


 The present version on the code is published under the GNU general public
licence and it is available at:

https://github.com/ManyBodyPhysics/LectureNotesPhysics/blob/master/doc/src/Chapter11-programs/Inf_Matter

 If you are using this software (or any part of it) for publication, we ask
you that its origin is formally acknowledged by citing the following chapter:

 C. Barbieri and A. Carbone, "Self-consistent Green's function approaches",
chapter 11 of Lecture Notes in Physics "An advanced course in computational
nuclear physics: Bridging the scales from quarks to neutron stars",
Edited by Hjorth-Jensen, M. P. Lombardo, U. van Kolck, Editors
ISBN:
http://arxiv.org/abs/1611.03923


The source files included in this distribution are:

   MtK-SpProppagator.h / MtK-SpProppagator.cpp
                                    --  class to store the spectral (or Lehmann)
                                       representation; this is meant for *both* the
                                       single-particle propagator and the self-energy.

   MtK-Potential_minnesota.h / MtK-Potential_minnesota.cpp 
                                    --  to generate the matrix elements of the
                                       interaction in the momentum basis.

   MtK-Physics_constants.h          --  global definitions of relevant physics constants.

   MtK-Global_data.h                --  global variables to be used throughout the code.

   MtK-Main.cpp                     --  This is the Main function.

   MtK-Bases.h / MtK-Bases.cpp      --  classes for the s.-p. and 2p1h/2h1p configurations.

   MtK-Utilities.cpp                --  utilities stuff.

   MtK-Dyson.cpp                    --  Diagonalization of Dyson matrix; this is where
                                       the action is ;-).

   MtK-Lanczos.cpp                  --  Lanczos algorithm to reduce the 'C' and 'D' 
                                       sub-matrices in the Dyson eigenvalue problem.


See folder ‘Example’ for a brief tutorial in running the code.


Compile (for example) as:  g++  -O3  MtK-*cpp  -llapack  -lblas -o MtK.exe


 You will need to link to a BLAS and LAPACK library to compile this program.


(c) C. Barbieri and A. Carbone,   Surrey/Darmstadt,   December 2016.



