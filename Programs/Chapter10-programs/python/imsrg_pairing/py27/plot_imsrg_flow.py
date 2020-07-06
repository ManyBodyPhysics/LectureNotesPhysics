#!/usr/bin/env python

#------------------------------------------------------------------------------
# plot_imsrg_flow.py
#
# author:   H. Hergert 
# version:  1.0.1
# date:     Jul 6, 2020
# 
# tested with Python v2.7 and v3.7
# 
#------------------------------------------------------------------------------

from sys import argv

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import SymLogNorm, Normalize
from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable

import numpy as np
from numpy import array, dot, diag, reshape

#------------------------------------------------------------------------------
# plot helpers
#------------------------------------------------------------------------------
# format tick labels using LaTeX-like math fonts
def myLabels(x, pos):
  return '$%s$'%x

def myLogLabels(x, pos):
  return '$10^{%d}$'%(np.log10(x))

# save these settings for use in both following plots
def myPlotSettings(ax):
  ax.minorticks_on()
  ax.tick_params(axis='both',which='major',width=1.5,length=8)
  ax.tick_params(axis='both',which='minor',width=1.5,length=5)
  ax.tick_params(axis='both',width=2,length=10,labelsize=20)
  for s in ['left', 'right', 'top', 'bottom']:
    ax.spines[s].set_linewidth(2)
  return

#------------------------------------------------------------------------------
# plot flow
#------------------------------------------------------------------------------
def plot_energies(data, exact, filename):

  # diagonals vs. eigenvalues on absolute scale
  fig, ax = plt.subplots()

  plt.semilogx([1.0e-8,1.0e-4,1.0,100], [exact,exact,exact,exact], linewidth=2, 
    color='black', linestyle='dashed', dashes=(10,5))

  plt.semilogx(data[:,0], data[:,1], color='blue', marker='o', markersize=9, label='$E$')
  plt.semilogx(data[:,0], data[:,1]+data[:,2], color='red', marker='s', markersize=9, label='$+\Delta E^{(2)}$')
  plt.semilogx(data[:,0], data[:,1]+data[:,2]+data[:,3], color='green', marker='D', markersize=9,label='$+\Delta E^{(3)}$')

  myPlotSettings(ax)
  ax.xaxis.set_major_formatter(FuncFormatter(myLogLabels))
  ax.yaxis.set_major_formatter(FuncFormatter(myLabels))
  ax.set_xlim([0.00006,13])

  ymin,ymax=ax.get_ylim()
  ax.set_ylim(ymin-0.005,ymax+0.005)

  plt.xlabel('$s$', fontsize=20)
  plt.ylabel('$E\,\mathrm{[a.u.]}$', fontsize=20)

  # plt.legend(bbox_to_anchor=(0.35, 0.05), loc=3, borderaxespad=0.5)
  plt.legend(loc=1, borderaxespad=0.5)

  plt.savefig("%s.pdf"%(filename), bbox_inches="tight", pad_inches=0.05)
  plt.show()
  plt.close()

  return

def plot_norms_loglog(data, filename):

  # diagonals vs. eigenvalues on absolute scale
  fig, ax = plt.subplots()

  plt.loglog(data[:,0], data[:,6], basex=10, color='blue', marker='o', markersize=9, label='$||\eta||$')
  plt.loglog(data[:,0], data[:,8], basex=10, color='red',  marker='s', markersize=9, label='$||\Gamma_{od}||$')

  myPlotSettings(ax)
  ax.xaxis.set_major_formatter(FuncFormatter(myLogLabels))
  ax.yaxis.set_major_formatter(FuncFormatter(myLogLabels))
  plt.xlabel('$s$', fontsize=20)
  plt.ylabel('$||\eta||, ||\Gamma_{od}||\, [\mathrm{a.u.}]$', fontsize=20)

  plt.legend(bbox_to_anchor=(0.05, 0.05), loc=3, borderaxespad=0.5)
  
  plt.savefig("%s.norms.pdf"%(filename.rsplit(".",1)[0]), bbox_inches="tight", pad_inches=0.05)
  plt.show()
  plt.close()

  return

def plot_norms_semilog(data, filename):

  # diagonals vs. eigenvalues on absolute scale
  fig, ax = plt.subplots()

  plt.semilogy(data[:,0], data[:,6], basey=10, color='blue', marker='o', markersize=9, label='$||\eta||$')
  plt.semilogy(data[:,0], data[:,8], basey=10, color='red',  marker='s', markersize=9, label='$||\Gamma_{od}||$')

  myPlotSettings(ax)
  ax.xaxis.set_major_formatter(FuncFormatter(myLabels))
  ax.yaxis.set_major_formatter(FuncFormatter(myLogLabels))
  plt.xlabel('$s$', fontsize=20)
  plt.ylabel('$||\eta||, ||\Gamma_{od}||\, [\mathrm{a.u.}]$', fontsize=20)

  plt.legend(bbox_to_anchor=(0.05, 0.05), loc=3, borderaxespad=0.5)
  
  plt.savefig("%s.norms.semilog.pdf"%(filename.rsplit(".",1)[0]), bbox_inches="tight", pad_inches=0.05)
  plt.show()
  plt.close()

  return


#------------------------------------------------------------------------------
# main
#------------------------------------------------------------------------------
def main():
  filename = argv[1]
  exact    = argv[2]

  # read data from file
  data = np.loadtxt(filename, skiprows=2)  

  plot_energies(data, exact, filename)
  plot_norms_loglog(data,filename)
  plot_norms_semilog(data,filename)

  return

#------------------------------------------------------------------------------
# make executable
#------------------------------------------------------------------------------
if __name__ == "__main__": 
  main()
