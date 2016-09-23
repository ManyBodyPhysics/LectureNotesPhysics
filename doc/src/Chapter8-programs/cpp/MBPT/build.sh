#!/bin/sh

module load openblas

make clean
PROFILE=1 make
