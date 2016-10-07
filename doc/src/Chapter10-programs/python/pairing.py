#!/usr/bin/env python

from numpy import array, dot, diag, reshape
from scipy.linalg import eigvalsh
from scipy.integrate import odeint

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


H0  = Hamiltonian(1.0,0.5)

print eigvalsh(H0)

dim = 6
y0  = reshape(H0, -1)

ys  = odeint(derivative, y0, array([0.,1.,2.,5.]), args=(dim,))
Hs  = reshape(ys, (-1, dim,dim))

print Hs[-1]
print eigvalsh(Hs[-1])
