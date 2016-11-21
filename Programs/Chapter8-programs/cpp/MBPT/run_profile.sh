#!/bin/sh

module load openblas

method=3
function=2
nmax=28
OMP_NUM_THREADS=1 ./mbpt.x 0.08 $nmax 2 2 ${method} ${function}
gprof ./mbpt.x > profile.out
