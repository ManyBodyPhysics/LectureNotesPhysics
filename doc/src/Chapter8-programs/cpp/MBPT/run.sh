#!/bin/sh

module load openblas

for method in 2 3; do
    for function in 0 1 2 3 4 5 6; do
        OMP_NUM_THREADS=2 ./mbpt.x 0.08 28 0 2 ${method} ${function}
        echo
    done
done
