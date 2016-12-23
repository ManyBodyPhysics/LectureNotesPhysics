#!/usr/bin/env python

#------------------------------------------------------------------------------
# srg_nn.py
#
# author:   H. Hergert 
# version:  1.1.0
# date:     Nov 21, 2016
# 
# tested with Python v2.7
# 
# SRG evolution of a chiral NN interaction with cutoff Lambda in the deuteron 
# partial waves, using a Gauss-Legendre momentum mesh.
#
#------------------------------------------------------------------------------


import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import SymLogNorm, Normalize
from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable

import numpy as np
from numpy import array, dot, diag, reshape, sqrt
from math import sqrt, pi
from scipy.linalg import eigvalsh
from scipy.integrate import ode


#------------------------------------------------------------------------------
# constants
#------------------------------------------------------------------------------
hbarm = 41.4710570772


#------------------------------------------------------------------------------
# helpers
#------------------------------------------------------------------------------
def find_nearest(array, value):
    distance  = np.absolute(array-value)
    indices = np.where(distance == np.min(distance))
    return indices[0]


#------------------------------------------------------------------------------
# plot matrix snapshots
#------------------------------------------------------------------------------

def plot_snapshots(Hs, flowparams, momenta, qMax):
    fig  = plt.figure(1, (10., 50.))

    nplots = len(flowparams)
    ncols  = 1
    nrows  = nplots

    grid = AxesGrid(fig, 111,                       # similar to subplot(111)
                     nrows_ncols=(nrows, ncols),    # creates grid of axes
                     axes_pad=1.,                   # pad between axes in inch.
                     label_mode='all',              # put labels on left, bottom
                     cbar_mode='each',              # color bars 
                     cbar_pad=0.20,                 # insert space between plots and color bar
                     cbar_size='10%'                # size of colorbar relative to last image
                     )
    hmax = 0.0
    hmin = 0.0
    for h in Hs:
      hmax = max(hmax, np.ma.max(h))
      hmin = min(hmin, np.ma.min(h))

    # get indices of max. momenta
    cmax, ccmax  = find_nearest(momenta, qMax)
    edge         = len(momenta)/2

    # create individual snapshots - figures are still addressed by single index,
    # despite multi-row grid
    for s in range(Hs.shape[0]):
        h = np.vstack((np.hstack((Hs[s,0:cmax,0:cmax], Hs[s,0:cmax,edge:ccmax])), 
                     np.hstack((Hs[s,edge:ccmax,0:cmax], Hs[s,edge:ccmax,edge:ccmax]))
                    ))
        img = grid[s].imshow(h,
                              cmap=plt.get_cmap('RdBu_r'),                            # choose color map
                              interpolation='bicubic',
                              # filterrad=10,
                              norm=SymLogNorm(linthresh=0.0001, vmax=2., vmin=-2.0),   # normalize 
                              vmin=-2.0,                                               # min/max values for data
                              vmax=2.0
                              )

    # contours
    levels = np.arange(-2, 1, 0.12)
    grid[s].contour(h, levels, colors='black', ls="-", origin='lower',linewidths=1)

    # plot labels, tick marks etc.
    grid[s].set_title('$\\lambda=%s\,\mathrm{fm}^{-1}$'%flowparams[s])
    
    grid[s].set_xticks([0,20,40,60,80,100,120,140,160])
    grid[s].set_yticks([0,20,40,60,80,100,120,140,160])
    grid[s].set_xticklabels(['$0$','$1.0$','$2.0$','$3.0$','$4.0$','$1.0$','$2.0$','$3.0$','$4.0$'])
    grid[s].set_yticklabels(['$0$','$1.0$','$2.0$','$3.0$','$4.0$','$1.0$','$2.0$','$3.0$','$4.0$'])

    grid[s].tick_params(axis='both',which='major',width=1.5,length=5)
    grid[s].tick_params(axis='both',which='minor',width=1.5,length=5)

    grid[s].axvline(x=[80],color="black", ls="--")
    grid[s].axhline(y=[80],color="black", ls="--")
    grid[s].xaxis.set_label_text("$q\,[\mathrm{fm}^{-1}]$")
    grid[s].yaxis.set_label_text("$q'\,[\mathrm{fm}^{-1}]$")
    
    # color bar
    cbar = grid.cbar_axes[s]
    plt.colorbar(img, cax=cbar, 
      ticks=[ -2.0, -1.0, -1.0e-1, -1.0e-2, -1.0e-3, -1.0e-4, 0., 
         1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0]
      )
    cbar.axes.set_yticklabels(['$-2.0$', '$-1.0$', '$-10^{-1}$', '$-10^{-2}$', 
                               '$-10^{-3}$', '$-10^{-4}$','$0.0$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', 
                               '$10^{-1}$', '$1.0$'])
    cbar.set_ylabel("$V(q,q')\,\mathrm{[fm]}$") 

    # save figure
    plt.savefig("srg_n3lo500.pdf", bbox_inches="tight", pad_inches=0.05)
    plt.savefig("srg_n3lo500.png", bbox_inches="tight", pad_inches=0.05)
    #plt.show()

    return

#------------------------------------------------------------------------------
# matrix element I/O, mesh functions
#------------------------------------------------------------------------------
def uniform_weights(momenta):
    weights = np.ones_like(momenta)
    weights *= abs(momenta[1]-momenta[0])
    return weights

def read_mesh(filename):
    data = np.loadtxt(filename, comments="#")  
    dim  = data.shape[1]
    
    momenta = data[0,:dim]

    return momenta
  
def read_interaction(filename):
    data = np.loadtxt(filename, comments="#")  
    dim  = data.shape[1]
    V = data[1:,:dim]
    return V

#------------------------------------------------------------------------------
# commutator
#------------------------------------------------------------------------------
def commutator(a,b):
    return dot(a,b) - dot(b,a)

#------------------------------------------------------------------------------
# flow equation (right-hand side)
#------------------------------------------------------------------------------
def derivative(lam, y, T):
    dim = T.shape[0]

    # reshape the solution vector into a dim x dim matrix
    V = reshape(y, (dim, dim))

    # calculate the generator
    eta = commutator(T, V)

    # dV is the derivative in matrix form 
    dV  = -4.0/(lam**5) * commutator(eta, T+V)

    # convert dH into a linear array for the ODE solver
    dy = reshape(dV, -1)
      
    return dy


#------------------------------------------------------------------------------
# Main program
#------------------------------------------------------------------------------
def main():

    # duplicate the mesh points (and weights, see below) because we have a 
    # coupled-channel problem
    mom_tmp = read_mesh("n3lo500_3s1.meq")
    momenta = np.concatenate([mom_tmp,mom_tmp])
    weights = uniform_weights(momenta)
    dim     = len(momenta)

    # set up p^2 (kinetic energy in units where h^2/2\mu = 1)
    T = diag(momenta*momenta)

    # set up interaction matrix in coupled channels:
    #
    #    /   V_{3S1}            V_{3S1-3D1}  \
    #    \   V_{3S1-3D1}^\dag   V_{3D1}      /
    
    # read individual partial waves
    partial_waves=[]
    for filename in ["n3lo500_3s1.meq", "n3lo500_3d1.meq", "n3lo500_3sd1.meq"]:
      partial_waves.append(read_interaction(filename))
      # print partial_waves[-1].shape

    # assemble coupled channel matrix
    V = np.vstack((np.hstack((partial_waves[0], partial_waves[2])), 
                   np.hstack((np.transpose(partial_waves[2]), partial_waves[1]))
                  ))

    # switch to scattering units
    V = V/hbarm
  
    # set up conversion matrix for V: this is used to absorb momentum^2 and
    # weight factors into V, so that we can use the commutator routines for
    # eta and the derivative as is
    conversion_matrix = np.zeros_like(T)
    for i in range(dim):
        for j in range(dim):
            # Regularize the conversion matrix at zero momentum - set elements
            # to machine precision so we can invert the matrix for plots etc.
            # Note that momentum values are positive, by construction.
            qiqj = max(np.finfo(float).eps, momenta[i]*momenta[j])
            conversion_matrix[i,j] = qiqj*sqrt(weights[i]*weights[j])

    V *= conversion_matrix
    
    # turn initial interaction into a linear array
    y0  = reshape(V, -1)                 


    # flow parameters for snapshot images - the initial lambda should be
    # infinity, we use something reasonably large
    lam_initial = 20.0
    lam_final = 1.5

    # integrate using scipy.ode instead of scipy.odeint - this gives
    # us more control over the solver
    solver = ode(derivative,jac=None)

    # equations may get stiff, so we use VODE and Backward Differentiation
    solver.set_integrator('vode', method='bdf', order=5, nsteps=1000)
    solver.set_f_params(T)
    solver.set_initial_value(y0, lam_initial)


    print("%-8s   %-14s"%("s", "E_deuteron [MeV]"))
    print("-----------------------------------------------------------------------------------------------------------------")

    # calculate exact eigenvalues
    print("%8.5f %14.8f"%(solver.t, eigvalsh((T + V)*hbarm)[0]))

    flowparams=([lam_initial])
    Vs=([V])
    while solver.successful() and solver.t > lam_final:
        # adjust the step size in different regions of the flow parameter
        if solver.t >= 6.0:
            ys = solver.integrate(solver.t-1.0)
        elif solver.t < 6.0 and solver.t >= 2.5:
            ys = solver.integrate(solver.t-0.5)
        elif solver.t < 2.5 and solver.t >= lam_final:
            ys = solver.integrate(solver.t-0.1)

    # add evolved interactions to the list
    flowparams.append(solver.t)
    Vtmp = reshape(ys,(dim,dim))
    Vs.append(Vtmp)

    print("%8.5f %14.8f"%(solver.t, eigvalsh((T + Vtmp)*hbarm)[0]))

    # generate snapshots of the evolution
    plot_snapshots((Vs[-2:]/conversion_matrix), flowparams[-2:], momenta, 4.0)

    return

#------------------------------------------------------------------------------
# make executable
#------------------------------------------------------------------------------
if __name__ == "__main__": 
    main()