#!/usr/bin/env python

#------------------------------------------------------------------------------
# plot_correlation_energy.py
#
# author:   H. Hergert 
# version:  1.0.0
# date:     Dec 6, 2016
# 
# tested with Python v2.7
# 
#------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import glob

from   pylab import *
from   matplotlib import rc

rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)

def Hamiltonian(delta,g):

    H = array(
        [[2*delta-g,    -0.5*g,     -0.5*g,     -0.5*g,    -0.5*g,          0.],
         [   -0.5*g, 4*delta-g,     -0.5*g,     -0.5*g,        0.,     -0.5*g ], 
         [   -0.5*g,    -0.5*g,  6*delta-g,         0.,    -0.5*g,     -0.5*g ], 
         [   -0.5*g,    -0.5*g,         0.,  6*delta-g,    -0.5*g,     -0.5*g ], 
         [   -0.5*g,        0.,     -0.5*g,     -0.5*g, 8*delta-g,     -0.5*g ], 
         [       0.,    -0.5*g,     -0.5*g,     -0.5*g,    -0.5*g, 10*delta-g ]]
      )

    return H

def myLabels(x, pos):
    return '$%s$'%x

# pairing strengths
glist  = [ -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,  0.0,  0.1,  0.2,  0.3,  
            0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0 ]

# uncorrelated energy
uncorr = [  3.0, 2.9, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 
            1.3, 1.2, 1.1, 1.0 ]

# exact
eigenvals  = [eigvalsh(Hamiltonian(1.0,g))[0] for g in glist]
exact      = [ a - b for a, b in zip(eigenvals,uncorr)]
  

# CCD correlation energy
ccd    = [ -0.21895200, -0.18230800, -0.14813200, -0.11667400, -0.08821240, -0.06305620, -0.04154780,
           -0.02406310, -0.01101140, -0.00283397,  0.00000000, -0.00300047, -0.01233930, -0.02852180,
           -0.05204090, -0.08336220, -0.12290800, -0.17104400, -0.22806500, -0.29418800, -0.36955000 ]

# MBPT4 correlation energy
mbpt4  = [ -0.44400000, -0.28809800, -0.19581600, -0.13686300, -0.09602420, -0.06570970, -0.04228900,
           -0.02421400, -0.01102850, -0.00283444,  0.00000000, -0.00300012, -0.01232950, -0.02845710,
           -0.05180360, -0.08273260, -0.12154800, -0.16849800, -0.22377500, -0.28752700, -0.35985600 ]

# IMSRG, MBPT2, MBPT3
mbpt2  = [ ]
mbpt3  = [ ]
white  = [ ]
wegner = [ ]
imtime = [ ]

# White generator flows - the MBPT numbers are contained in all flows, regardless of generator

for g in glist:
    filename = glob.glob("results/imsrg-white*g%+3.1f*.flow"%(g))[0]
    data = np.loadtxt(filename, skiprows=2)
    if g != 0.0:
      	mbpt2.append(data[0,2])		        	  # correlation energy: just MBPT2
      	mbpt3.append(data[0,2] + data[0,3])		# correlation energy: MBPT2+3 corrections
      	white.append(data[-1,1] - data[0,1])
    else:
      	mbpt2.append(0.0)
      	mbpt3.append(0.0)
      	white.append(0.0)

# Wegner generator flows
for g in glist:
    filename = glob.glob("results/imsrg-wegner*g%+3.1f*.flow"%(g))[0]
    data = np.loadtxt(filename, skiprows=2)
    if g != 0.0:
      	wegner.append(data[-1,1] - data[0,1])
    else:
        wegner.append(0.0)

# imaginary time generator flows
for g in glist:
    filename = glob.glob("results/imsrg-imtime*g%+3.1f*.flow"%(g))[0]
    data = np.loadtxt(filename, skiprows=2)
    if g != 0.0:
      	imtime.append(data[-1,1] - data[0,1])
    else:
        imtime.append(0.0)



#------------------------------------------------------------------------------
# Comparison of methods
#------------------------------------------------------------------------------
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig = figure(figsize=(8,6))

# pl, ax = plt.subplots()
ax = fig.gca()
ax.tick_params(axis='both',which='major',width=1.5,length=8)
ax.tick_params(axis='both',which='minor',width=1.5,length=5)
ax.tick_params(axis='both',width=2,length=10,labelsize=20)
ax.minorticks_on()
for s in ['left', 'right', 'top', 'bottom']:
  	ax.spines[s].set_linewidth(2)
ax.set_xlim([-1.05,1.05])  
ax.set_ylim([-0.5,0.06])  

ax.xaxis.set_major_formatter(FuncFormatter(myLabels))
ax.yaxis.set_major_formatter(FuncFormatter(myLabels))

plt.xlabel('$g\,\mathrm{[a.u.]}$', fontsize=20)
plt.ylabel('$E_\mathrm{corr}\, \mathrm{[a.u.]}$', fontsize=20)
pl_exact = plt.plot(glist, exact, color='black',linestyle='-',  linewidth = 2.0, label = 'exact')
pl_mbpt2 = plt.plot(glist, mbpt2, marker='^', markersize=8, color='gold', linestyle='-',  linewidth = 2.0, label = 'MBPT(2)')
pl_mbpt3 = plt.plot(glist, mbpt3, marker='v', markersize=8, color='orange', linestyle='-',  linewidth = 2.0, label = 'MBPT(3)')
pl_mbpt4 = plt.plot(glist, mbpt4, marker='D', markersize=8, color='red',    linestyle='--',  linewidth = 2.0, label = 'MBPT(4)')
pl_ccd   = plt.plot(glist,   ccd, marker='s', markersize=8, color='green',  dashes=[8,6], linewidth = 2.0, label = 'CCD')
pl_white = plt.plot(glist, white, marker='o', markersize=8, color='blue',   linestyle='-', linewidth = 2.0, label = 'IMSRG(2)')

plt.legend(bbox_to_anchor=(0.35, 0.05), loc=3, borderaxespad=0.5)
plt.savefig("correlation_energy.pdf", bbox_inches="tight", pad_inches=0.05)
plt.show()
plt.close()

#------------------------------------------------------------------------------
# Comparison of generators
#------------------------------------------------------------------------------
fig = figure(figsize=(8,6))

# pl, ax = plt.subplots()
ax = fig.gca()
ax.tick_params(axis='both',which='major',width=1.5,length=8)
ax.tick_params(axis='both',which='minor',width=1.5,length=5)
ax.tick_params(axis='both',width=2,length=10,labelsize=20)
ax.minorticks_on()
for s in ['left', 'right', 'top', 'bottom']:
    ax.spines[s].set_linewidth(2)
ax.set_xlim([-1.05,1.05])  
ax.set_ylim([-0.5,0.02])  

ax.xaxis.set_major_formatter(FuncFormatter(myLabels))
ax.yaxis.set_major_formatter(FuncFormatter(myLabels))

plt.xlabel('$g\,\mathrm{[a.u.]}$', fontsize=20)
plt.ylabel('$E_\mathrm{corr}\, \mathrm{[a.u.]}$', fontsize=20)
pl_exact  = plt.plot(glist, exact, color='black',linestyle='-',  linewidth = 2.0, label = 'exact')
pl_imtime = plt.plot(glist, imtime, marker='D', markersize=8, color='green', linestyle='solid',  linewidth = 2.0, label = 'imag. time')
pl_wegner = plt.plot(glist, wegner, marker='s', markersize=8, color='red', dashes=[8,6],  linewidth = 2.0, label = 'Wegner')
pl_white  = plt.plot(glist, white, marker='o', markersize=8, color='blue', dashes=[2,2,4],linewidth = 2.0, label = 'White')

plt.legend(bbox_to_anchor=(0.05, 0.05), loc=3, borderaxespad=0.5)
plt.savefig("correlation_energy_generators.pdf", bbox_inches="tight", pad_inches=0.05)
# plt.show()
