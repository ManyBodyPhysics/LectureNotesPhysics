#!/bin/sh

module load openblas

method=3
function=5
for nmax in 10 28; do
    for nn in 1 2 3; do
        OMP_NUM_THREADS=1 ./mbpt.x 0.08 $nmax 0 $nn ${method} ${function}
        echo
    done
    for nn in 1 2 3; do
        OMP_NUM_THREADS=1 ./mbpt.x 0.08 $nmax $nn $nn ${method} ${function}
        echo
    done
done
