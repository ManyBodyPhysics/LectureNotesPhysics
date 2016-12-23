 


 "Green's function and coupled cluster calculations of spectral function
for a pairing model"


 The code in this folder is based on the lectures of the TALENT school, course #2,
held at GANIL on July 6th-24th, 2015 and it has been later improved to perform
the benchmark calculations presented in chapters 8 and 11 of the Lecture Notes
in Physics


 This version uses the SCGF and CC methods to solve the toy pairing model
with 4 levels and 4 spin-1/2 particles that was introduced during the first
week of the Talent lectures (see
 http://nucleartalent.github.io/Course2ManyBodyMethods/doc/web/course.html )
and in chapters 8 and 11 of the above book. The present version on the code
is publishes under the GNU general public licence and it is available at:
https://github.com/ManyBodyPhysics/LectureNotesPhysics/blob/master/doc/src/Chapter11-programs/Pair_Model

 If you are using this software (or any part of it) for publication, we ask
you that its origin is formally acknowledged by citing the following chapter:

 C. Barbieri and A. Carbone, "Self-consistent Green's function approaches",
chapter 11 of Lecture Notes in Physics "An advanced course in computational
nuclear physics: Bridging the scales from quarks to neutron stars",
Edited by Hjorth-Jensen, M. P. Lombardo, U. van Kolck, Editors
ISBN:
http://arxiv.org/abs/1611.03923




 Calculations can be done at different levels of SCGF theory by using the ADC(2)
approach at second order, by adding all order resummations of the two-particle
and two-hole ladders [the so-called 2p1h-TDA, which is an extension of ADC(2)],
in the full ADC(3) approximation and finally by adding CCD corrections to
obtain the ADC(3)-D scheme.

 It also iterates the energy independent self-energy (i.e. the 'correlated HF'
diagram with the dressed propagator to achieve self-consistency at the sc0
level.

 The coupled cluster equations can be solved by either including on only pp
and hh ladders or in the full CCD approximation scheme.

 Note that for this particular model does not allow for ph (ring) TDA-like
resummation in any of the methods.


 The ground state is then calculated in FCI to obtain the exact fundamental
energy and correlation energy, for comparison with the Koltun Sr result. It
is up to you, as an exercise, to do FCI for the N=3 and N=5 particle systems
and to extract the exact spectral function to compare with the results from
the ADC(n) method.

 The coupled cluster equations are also solved to resum the pp and hh ladders.
These diagrams are a subset of the full CCD method that could be compared to the
correlations added in the 2p1h-ADC(2).

 The particular pairing Hamiltonian used in this exercise allows us to use
the following two simplifications:

 - The spin-up and spin-down parts of the Dyson equation completely decouple
  and are equal because of symmetry. Hence, we need to calculate only one
  of the two.

 - The particle-hole contribution from the pairing interactions is zero in
  the ring summation. So, we don't need to worry about rings.




Compile as:    gfortran  Pairing_GF_CC.f90  -llapack -o Pairing.exe


 You will need to link to the LAPACK library to compile this program.


(c) C. Barbieri and A. Carbone,   Surrey/Darmstadt,   October 2016.



