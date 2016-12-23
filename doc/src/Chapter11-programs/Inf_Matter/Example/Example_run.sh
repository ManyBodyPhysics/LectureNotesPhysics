#!/bin/bash

#
#  This file is contains a simple example tutorial on how to run the
# 'MattersK' code and generate a 3D figure of the spectral function.
#
#  Read the explanation below to run the code in trearctively. Otherwise
# the whole calculation can be run automatically by executing this file.
# Type 'source Example_run.sh' on a terminal to do that.
#
#
#  Here we assume that the code hase been compiled with the name 'MtK.exe'
# and that ir is places in the parent directory (../MtK.exe).  Once the code
# is started, it asks for input parameters from terminal,  we will make a file
# called 'SNM_input' with all the required imputs.
#
#  First thing the user is asked to build the basis as follows:
#
#  "Enter in order: max for nsq = sum_i ni^2 ,  n_max, and chrg_max (0 for PNM and 1for SNM) ?"
#
#  There are two ways to constrain the basis: (1) to set the maximum value of N_{sq} (see right
# below Eq.(11.55)) which is a spherical truncation in momentum space; and (2) to set the maximum
# absolute values of n_x, n_y or n_z, this is a cubic box truncation in k-space.
# The code can do both and asks for both parameters. To use one simply put a slightly larger value
# for the other constraint.
#
#  The minimum charge os a nucleon in 0 (fpr neutron). So chrg_max=0 for pure neutron matter
# and =1 for symmetric nucleon matter (both protons and neutrons).
#
#  There we want to try a calculation of SNM in N_{sq}=10, thus:
#

rm -f  SNM_input
echo ' 10 ' >> SNM_input
echo ' 20 ' >> SNM_input
echo '  1 ' >> SNM_input

#
#  Next is the muber of occupied orbits:
#
#  "Number of nucleons in the Fermi Sea? "
#
#  We chose 76, whic is 3 values of |k| for the Fermi sea, 38 protons and 38 neutrons (see
# Table 8.1)
#

echo ' 76 ' >> SNM_input

#
#  Next is the density, we chose nominal saturation density:
#
#  "Density [in fm^-3]? "
#

echo ' 0.16 ' >> SNM_input

#
#   At this point the code replies:
#
#  "density = 0.16/fm^3 ,     K_F = 1.33302/fm ,   L_box = 7.80245 fm"
#
#  Then it builds the sp basis and 2p1h/2h1p ISCs...    and return the HF energy
# and free kinetic energy:
#
#  "E_kin/A = 21.2169 MeV   (kinetic energy of free Fermi gas)"
#  "E_HF/A  = -24.9396 MeV   (uncorrelated HF energy)"
#
#
#  The you are asked whether you want the files with the self-energy, spectral function, etcc...
#
#  "Would you like to generate data files for plotting the spectral function and self energy [y/Y or n/N] ? "
#
#  We willmake a plot f these, so say yes:
#

echo ' Y ' >> SNM_input

#
#  The last question is how many Lanczos iterations. About 300 is even more than redundant in
# most cases (but it is goo to  play with it and try many options)...
#
#
#  "How many Lanczos iterations [< 0 for no Lanczos] ? "
#

echo ' 300 ' >> SNM_input
echo  >> SNM_input

#
#  At this point the code diagonalizes Dyson for all independent parial waves and then output
# the results.
#
#   Let's run it:
#

echo
echo
echo ' Now running Mtk.exe (it will take a minute or two...)'
echo

time ../MtK.exe < SNM_input > SNM_output

#
#  The results for the total number of particle in the box, kientin energy an dKoltun SR
# are plotted last:
#

echo
echo
echo ' Final results for the SNM test are as follows (see file SNM_output for full output):'
echo ' ------------------------------------------------------------------------------------'
echo
echo

tail -n18  SNM_output

#
#  A large number of files have been created, one for each k channel in the basis. Those that
# allow to make the 3D plot are called 'SpectFunct**3D**.dat'. The following gnuplot file is
# an example og how to get the plot:
#

gnuplot Plot_3D_SpectFunct_SNM.gpl  --persist

epspdf spect_funct_3D_snm.eps







