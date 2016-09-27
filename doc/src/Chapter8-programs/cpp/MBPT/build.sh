#!/bin/sh

module load openblas

make clean
OPTIMIZE=1 make
