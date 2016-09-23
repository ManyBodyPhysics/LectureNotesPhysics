#!/bin/sh

module load openblas

method=2

for function in 0 1 2 3 4 5; do
    OMP_NUM_THREADS=1 ./MBPT 0.08 28 0 2 2 ${function}
    echo
done
