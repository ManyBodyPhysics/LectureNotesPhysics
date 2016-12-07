#!/usr/bin/env python

#------------------------------------------------------------------------------
# srg_pairing.py
#
# author:   H. Hergert 
# version:  1.1.0
# date:     Nov 18, 2016
# 
# tested with Python v2.7
# 
# Solves the pairing model for four particles in a basis of four doubly 
# degenerate states by means of a Similarity Renormalization Group (SRG)
# flow.
#
#------------------------------------------------------------------------------


import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import SymLogNorm, Normalize
from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable

import numpy as np
from numpy import array, dot, diag, reshape
from scipy.linalg import eigvalsh
from scipy.integrate import odeint

#------------------------------------------------------------------------------
# plot helpers
#------------------------------------------------------------------------------
# format tick labels using LaTeX-like math fonts
def myLabels(x, pos):
  return '$%s$'%x

# save these settings for use in both following plots
def myPlotSettings(ax, formatter):
  ax.xaxis.set_major_formatter(formatter)
  ax.yaxis.set_major_formatter(formatter)
  ax.tick_params(axis='both',which='major',width=1.5,length=8)
  ax.tick_params(axis='both',which='minor',width=1.5,length=5)
  ax.tick_params(axis='both',width=2,length=10,labelsize=20)
  for s in ['left', 'right', 'top', 'bottom']:
    ax.spines[s].set_linewidth(2)
  ax.set_xlim([0.0007,13])  
  return

#------------------------------------------------------------------------------
# plot eigenvalues and diagonals
#------------------------------------------------------------------------------
def plot_diagonals(data, eigenvalues, flowparams, delta, g):
  dim       = len(data)
  formatter = FuncFormatter(myLabels)
  markers   = ['o' for i in range(dim)]
  cols      = ['blue', 'red', 'purple', 'green', 'orange', 'deepskyblue']

  # diagonals vs. eigenvalues on absolute scale
  fig, ax = plt.subplots()
  for i in range(dim):
    plt.semilogx(flowparams, [eigenvalues[i] for e in range(flowparams.shape[0])], color=cols[i], linestyle='solid')
    plt.semilogx(flowparams, data[i], color=cols[i], linestyle='dashed', marker=markers[i], markersize=10)

  myPlotSettings(ax, formatter)

  plt.savefig("srg_pairing_diag_delta%2.1f_g%2.1f.pdf"%(delta, g), bbox_inches="tight", pad_inches=0.05)
  plt.show()

  # difference between diagonals and eigenvalues
  fig, ax = plt.subplots()
  for i in range(dim):
    plot_diff = plt.semilogx(flowparams, data[i]-eigenvalues[i], color=cols[i], linestyle='solid', marker=markers[i], markersize=10)

  myPlotSettings(ax, formatter)

  plt.savefig("srg_pairing_diag-eval_delta%2.1f_g%2.1f.pdf"%(delta, g), bbox_inches="tight", pad_inches=0.05)
  plt.show()
  return

#------------------------------------------------------------------------------
# plot matrix snapshots
#------------------------------------------------------------------------------
def plot_snapshots(Hs, flowparams, delta, g):
  fig  = plt.figure(1, (10., 5.))
  grid = AxesGrid(fig, 111,                       # similar to subplot(111)
                   nrows_ncols=(2, Hs.shape[0]/2),  # creates grid of axes
                   axes_pad=0.25,                 # pad between axes in inch.
                   label_mode='L',                # put labels on left, bottom
                   cbar_mode='single',            # one color bar (default: right of last image in grid)
                   cbar_pad=0.20,                 # insert space between plots and color bar
                   cbar_size='10%'                # size of colorbar relative to last image
                   )

  # create individual snapshots - figures are still addressed by single index,
  # despite multi-row grid
  for s in range(Hs.shape[0]):
    img = grid[s].imshow(Hs[s], 
                          cmap=plt.get_cmap('RdBu_r'),                                  # choose color map
                          interpolation='nearest',       
                          norm=SymLogNorm(linthresh=1e-10,vmin=-0.5*g,vmax=10*delta),   # normalize 
                          vmin=-0.5*g,                                                  # min/max values for data
                          vmax=10*delta
                          )

    # tune plots: switch off tick marks, ensure that plots retain aspect ratio
    grid[s].set_title('$s=%s$'%flowparams[s])
    grid[s].tick_params(      
      bottom='off',      
      top='off',
      left='off',      
      right='off'
      )
    grid[s].set_xticks([0,1,2,3,4,5])
    grid[s].set_yticks([0,1,2,3,4,5])
    grid[s].set_xticklabels(['$0$','$1$','$2$','$3$','$4$','$5$'])
    grid[s].set_yticklabels(['$0$','$1$','$2$','$3$','$4$','$5$'])

  cbar = grid.cbar_axes[0]
  plt.colorbar(img, cax=cbar, 
    ticks=[ -1.0e-1, -1.0e-3, -1.0e-5, -1.0e-7, -1.09e-9 , 0., 
             1.0e-9, 1.0e-7, 1.0e-5, 1.0e-3, 0.1, 10.0]
    )
  cbar.axes.set_yticklabels(['$-10^{-1}$', '$-10^{-3}$', '$-10^{-5}$', '$-10^{-7}$', 
                             '$-10^{-9}$', '$0.0$', '$10^{-9}$', '$10^{-7}$', '$10^{-5}$', 
                             '$10^{-3}$', '$10^{-1}$', '$10$'])
  cbar.set_ylabel('$\mathrm{[a. u.]}$') 


  plt.savefig("srg_pairing_delta%2.1f_g%2.1f.pdf"%(delta, g), bbox_inches="tight", pad_inches=0.05)
  plt.show()

  return

#------------------------------------------------------------------------------
# SRG 
#------------------------------------------------------------------------------

# Hamiltonian for the pairing model
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

# commutator of matrices
def commutator(a,b):
  return dot(a,b) - dot(b,a)

# derivative / right-hand side of the flow equation
def derivative(y, t, dim):

  # reshape the solution vector into a dim x dim matrix
  H = reshape(y, (dim, dim))

  # extract diagonal Hamiltonian...
  Hd  = diag(diag(H))

  # ... and construct off-diagonal the Hamiltonian
  Hod = H-Hd

  # calculate the generator
  eta = commutator(Hd, Hod)

  # dH is the derivative in matrix form 
  dH  = commutator(eta, H)

  # convert dH into a linear array for the ODE solver
  dydt = reshape(dH, -1)
    
  return dydt



#------------------------------------------------------------------------------
# Main program
#------------------------------------------------------------------------------

def main():
  g     = 0.5
  delta = 1

  H0    = Hamiltonian(delta, g)
  dim   = H0.shape[0]

  # calculate exact eigenvalues
  eigenvalues = eigvalsh(H0)
  print eigenvalues
  exit()

  # turn initial Hamiltonian into a linear array
  y0  = reshape(H0, -1)                 

  # flow parameters for snapshot images
  flowparams = array([0.,0.001,0.01,0.05,0.1, 1., 5., 10.])

  # integrate flow equations - odeint returns an array of solutions,
  # which are 1d arrays themselves
  ys  = odeint(derivative, y0, flowparams, args=(dim,))

  # reshape individual solution vectors into dim x dim Hamiltonian
  # matrices
  Hs  = reshape(ys, (-1, dim,dim))

  # print Hs[-1]
  # print eigvalsh(Hs[-1])

  data = []
  for h in Hs:
    data.append(diag(h))
  data = zip(*data)

  plot_diagonals(data, eigenvalues, flowparams, delta, g)
  plot_snapshots(Hs, flowparams, delta, g)

#------------------------------------------------------------------------------
# make executable
#------------------------------------------------------------------------------
if __name__ == "__main__": 
  main()
