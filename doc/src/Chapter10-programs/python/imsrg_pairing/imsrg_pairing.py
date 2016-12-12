#!/usr/bin/env python

#------------------------------------------------------------------------------
# imsrg_pairing.py
#
# author:   H. Hergert 
# version:  1.5.0
# date:     Dec 6, 2016
# 
# tested with Python v2.7
# 
# Solves the pairing model for four particles in a basis of four doubly 
# degenerate states by means of an In-Medium Similarity Renormalization 
# Group (IMSRG) flow.
#
#------------------------------------------------------------------------------

import numpy as np
from numpy import array, dot, diag, reshape, transpose
from scipy.linalg import eigvalsh
from scipy.integrate import odeint, ode

from sys import argv

#-----------------------------------------------------------------------------------
# basis and index functions
#-----------------------------------------------------------------------------------

def construct_basis_2B(holes, particles):
  basis = []
  for i in holes:
    for j in holes:
      basis.append((i, j))

  for i in holes:
    for a in particles:
      basis.append((i, a))

  for a in particles:
    for i in holes:
      basis.append((a, i))

  for a in particles:
    for b in particles:
      basis.append((a, b))

  return basis


def construct_basis_ph2B(holes, particles):
  basis = []
  for i in holes:
    for j in holes:
      basis.append((i, j))

  for i in holes:
    for a in particles:
      basis.append((i, a))

  for a in particles:
    for i in holes:
      basis.append((a, i))

  for a in particles:
    for b in particles:
      basis.append((a, b))

  return basis


#
# We use dictionaries for the reverse lookup of state indices
#
def construct_index_2B(bas2B):
  index = { }
  for i, state in enumerate(bas2B):
    index[state] = i

  return index



#-----------------------------------------------------------------------------------
# transform matrices to particle-hole representation
#-----------------------------------------------------------------------------------
def ph_transform_2B(Gamma, bas2B, idx2B, basph2B, idxph2B):
  dim = len(basph2B)
  Gamma_ph = np.zeros((dim, dim))

  for i1, (a,b) in enumerate(basph2B):
    for i2, (c, d) in enumerate(basph2B):
      Gamma_ph[i1, i2] -= Gamma[idx2B[(a,d)], idx2B[(c,b)]]

  return Gamma_ph

def inverse_ph_transform_2B(Gamma_ph, bas2B, idx2B, basph2B, idxph2B):
  dim = len(bas2B)
  Gamma = np.zeros((dim, dim))

  for i1, (a,b) in enumerate(bas2B):
    for i2, (c, d) in enumerate(bas2B):
      Gamma[i1, i2] -= Gamma_ph[idxph2B[(a,d)], idxph2B[(c,b)]]
  
  return Gamma

#-----------------------------------------------------------------------------------
# commutator of matrices
#-----------------------------------------------------------------------------------
def commutator(a,b):
  return dot(a,b) - dot(b,a)

#-----------------------------------------------------------------------------------
# norms of off-diagonal Hamiltonian pieces
#-----------------------------------------------------------------------------------
def calc_fod_norm(f, user_data):
  particles = user_data["particles"]
  holes     = user_data["holes"]
  
  norm = 0.0
  for a in particles:
    for i in holes:
      norm += f[a,i]**2 + f[i,a]**2

  return np.sqrt(norm)

def calc_Gammaod_norm(Gamma, user_data):
  particles = user_data["particles"]
  holes     = user_data["holes"]
  idx2B     = user_data["idx2B"]

  norm = 0.0
  for a in particles:    
    for b in particles:
      for i in holes:
        for j in holes:
          norm += Gamma[idx2B[(a,b)],idx2B[(i,j)]]**2 + Gamma[idx2B[(i,j)],idx2B[(a,b)]]**2

  return np.sqrt(norm)

#-----------------------------------------------------------------------------------
# occupation number matrices
#-----------------------------------------------------------------------------------
def construct_occupation_1B(bas1B, holes, particles):
  dim = len(bas1B)
  occ = np.zeros(dim)

  for i in holes:
    occ[i] = 1.

  return occ

# diagonal matrix: n_a - n_b
def construct_occupationA_2B(bas2B, occ1B):
  dim = len(bas2B)
  occ = np.zeros((dim,dim))

  for i1, (i,j) in enumerate(bas2B):
    occ[i1, i1] = occ1B[i] - occ1B[j]

  return occ


# diagonal matrix: 1 - n_a - n_b
def construct_occupationB_2B(bas2B, occ1B):
  dim = len(bas2B)
  occ = np.zeros((dim,dim))

  for i1, (i,j) in enumerate(bas2B):
    occ[i1, i1] = 1. - occ1B[i] - occ1B[j]

  return occ

# diagonal matrix: n_a * n_b
def construct_occupationC_2B(bas2B, occ1B):
  dim = len(bas2B)
  occ = np.zeros((dim,dim))

  for i1, (i,j) in enumerate(bas2B):
    occ[i1, i1] = occ1B[i] * occ1B[j]

  return occ

#-----------------------------------------------------------------------------------
# generators
#-----------------------------------------------------------------------------------
def eta_brillouin(f, Gamma, user_data):
  dim1B     = user_data["dim1B"]
  particles = user_data["particles"]
  holes     = user_data["holes"]
  idx2B     = user_data["idx2B"]

  # one-body part of the generator
  eta1B  = np.zeros_like(f)

  for a in particles:
    for i in holes:
      # (1-n_a)n_i - n_a(1-n_i) = n_i - n_a
      eta1B[a, i] =  f[a,i]
      eta1B[i, a] = -f[a,i]

  # two-body part of the generator
  eta2B = np.zeros_like(Gamma)

  for a in particles:
    for b in particles:
      for i in holes:
        for j in holes:
          val = Gamma[idx2B[(a,b)], idx2B[(i,j)]]

          eta2B[idx2B[(a,b)],idx2B[(i,j)]] = val
          eta2B[idx2B[(i,j)],idx2B[(a,b)]] = -val

  return eta1B, eta2B

def eta_imtime(f, Gamma, user_data):
  dim1B     = user_data["dim1B"]
  particles = user_data["particles"]
  holes     = user_data["holes"]
  idx2B     = user_data["idx2B"]

  # one-body part of the generator
  eta1B  = np.zeros_like(f)

  for a in particles:
    for i in holes:
      dE = f[a,a] - f[i,i] + Gamma[idx2B[(a,i)], idx2B[(a,i)]]
      val = np.sign(dE)*f[a,i]
      eta1B[a, i] =  val
      eta1B[i, a] = -val 

  # two-body part of the generator
  eta2B = np.zeros_like(Gamma)

  for a in particles:
    for b in particles:
      for i in holes:
        for j in holes:
          dE = ( 
            f[a,a] + f[b,b] - f[i,i] - f[j,j]  
            + Gamma[idx2B[(a,b)],idx2B[(a,b)]] 
            + Gamma[idx2B[(i,j)],idx2B[(i,j)]]
            - Gamma[idx2B[(a,i)],idx2B[(a,i)]] 
            - Gamma[idx2B[(a,j)],idx2B[(a,j)]] 
            - Gamma[idx2B[(b,i)],idx2B[(b,i)]] 
            - Gamma[idx2B[(b,j)],idx2B[(b,j)]] 
          )

          val = np.sign(dE)*Gamma[idx2B[(a,b)], idx2B[(i,j)]]

          eta2B[idx2B[(a,b)],idx2B[(i,j)]] = val
          eta2B[idx2B[(i,j)],idx2B[(a,b)]] = -val

  return eta1B, eta2B


def eta_white(f, Gamma, user_data):
  dim1B     = user_data["dim1B"]
  particles = user_data["particles"]
  holes     = user_data["holes"]
  idx2B     = user_data["idx2B"]

  # one-body part of the generator
  eta1B  = np.zeros_like(f)

  for a in particles:
    for i in holes:
      denom = f[a,a] - f[i,i] + Gamma[idx2B[(a,i)], idx2B[(a,i)]]
      val = f[a,i]/denom
      eta1B[a, i] =  val
      eta1B[i, a] = -val 

  # two-body part of the generator
  eta2B = np.zeros_like(Gamma)

  for a in particles:
    for b in particles:
      for i in holes:
        for j in holes:
          denom = ( 
            f[a,a] + f[b,b] - f[i,i] - f[j,j]  
            + Gamma[idx2B[(a,b)],idx2B[(a,b)]] 
            + Gamma[idx2B[(i,j)],idx2B[(i,j)]]
            - Gamma[idx2B[(a,i)],idx2B[(a,i)]] 
            - Gamma[idx2B[(a,j)],idx2B[(a,j)]] 
            - Gamma[idx2B[(b,i)],idx2B[(b,i)]] 
            - Gamma[idx2B[(b,j)],idx2B[(b,j)]] 
          )

          val = Gamma[idx2B[(a,b)], idx2B[(i,j)]] / denom

          eta2B[idx2B[(a,b)],idx2B[(i,j)]] = val
          eta2B[idx2B[(i,j)],idx2B[(a,b)]] = -val

  return eta1B, eta2B


def eta_white_mp(f, Gamma, user_data):
  dim1B     = user_data["dim1B"]
  particles = user_data["particles"]
  holes     = user_data["holes"]
  idx2B     = user_data["idx2B"]

  # one-body part of the generator
  eta1B  = np.zeros_like(f)

  for a in particles:
    for i in holes:
      denom = f[a,a] - f[i,i]
      val = f[a,i]/denom
      eta1B[a, i] =  val
      eta1B[i, a] = -val 

  # two-body part of the generator
  eta2B = np.zeros_like(Gamma)

  for a in particles:
    for b in particles:
      for i in holes:
        for j in holes:
          denom = ( 
            f[a,a] + f[b,b] - f[i,i] - f[j,j]  
          )

          val = Gamma[idx2B[(a,b)], idx2B[(i,j)]] / denom

          eta2B[idx2B[(a,b)],idx2B[(i,j)]] = val
          eta2B[idx2B[(i,j)],idx2B[(a,b)]] = -val

  return eta1B, eta2B

def eta_white_atan(f, Gamma, user_data):
  dim1B     = user_data["dim1B"]
  particles = user_data["particles"]
  holes     = user_data["holes"]
  idx2B     = user_data["idx2B"]

  # one-body part of the generator
  eta1B  = np.zeros_like(f)

  for a in particles:
    for i in holes:
      denom = f[a,a] - f[i,i] + Gamma[idx2B[(a,i)], idx2B[(a,i)]]
      val = 0.5 * np.arctan(2 * f[a,i]/denom)
      eta1B[a, i] =  val
      eta1B[i, a] = -val 

  # two-body part of the generator
  eta2B = np.zeros_like(Gamma)

  for a in particles:
    for b in particles:
      for i in holes:
        for j in holes:
          denom = ( 
            f[a,a] + f[b,b] - f[i,i] - f[j,j] 
            + Gamma[idx2B[(a,b)],idx2B[(a,b)]] 
            + Gamma[idx2B[(i,j)],idx2B[(i,j)]] 
            - Gamma[idx2B[(a,i)],idx2B[(a,i)]] 
            - Gamma[idx2B[(a,j)],idx2B[(a,j)]] 
            - Gamma[idx2B[(b,i)],idx2B[(b,i)]] 
            - Gamma[idx2B[(b,j)],idx2B[(b,j)]] 
          )

          val = 0.5 * np.arctan(2 * Gamma[idx2B[(a,b)], idx2B[(i,j)]] / denom)

          eta2B[idx2B[(a,b)],idx2B[(i,j)]] = val
          eta2B[idx2B[(i,j)],idx2B[(a,b)]] = -val

  return eta1B, eta2B


def eta_wegner(f, Gamma, user_data):

  dim1B     = user_data["dim1B"]
  holes     = user_data["holes"]
  particles = user_data["particles"]
  bas2B     = user_data["bas2B"]
  basph2B   = user_data["basph2B"]
  idx2B     = user_data["idx2B"]
  idxph2B   = user_data["idxph2B"]
  occB_2B   = user_data["occB_2B"]
  occC_2B   = user_data["occC_2B"]
  occphA_2B = user_data["occphA_2B"]


  # split Hamiltonian in diagonal and off-diagonal parts
  fd      = np.zeros_like(f)
  fod     = np.zeros_like(f)
  Gammad  = np.zeros_like(Gamma)
  Gammaod = np.zeros_like(Gamma)

  for a in particles:
    for i in holes:
      fod[a, i] = f[a,i]
      fod[i, a] = f[i,a]
  fd = f - fod

  for a in particles:
    for b in particles:
      for i in holes:
        for j in holes:
          Gammaod[idx2B[(a,b)], idx2B[(i,j)]] = Gamma[idx2B[(a,b)], idx2B[(i,j)]]
          Gammaod[idx2B[(i,j)], idx2B[(a,b)]] = Gamma[idx2B[(i,j)], idx2B[(a,b)]]
  Gammad = Gamma - Gammaod


  #############################        
  # one-body flow equation  
  eta1B  = np.zeros_like(f)

  # 1B - 1B
  eta1B += commutator(fd, fod)

  # 1B - 2B
  for p in range(dim1B):
    for q in range(dim1B):
      for i in holes:
        for a in particles:
          eta1B[p,q] += (
            fd[i,a]  * Gammaod[idx2B[(a, p)], idx2B[(i, q)]] 
            - fd[a,i]  * Gammaod[idx2B[(i, p)], idx2B[(a, q)]] 
            - fod[i,a] * Gammad[idx2B[(a, p)], idx2B[(i, q)]] 
            + fod[a,i] * Gammad[idx2B[(i, p)], idx2B[(a, q)]]
          )

  # 2B - 2B
  # n_a n_b nn_c + nn_a nn_b n_c = n_a n_b + (1 - n_a - n_b) * n_c
  GammaGamma = dot(Gammad, dot(occB_2B, Gammaod))
  for p in range(dim1B):
    for q in range(dim1B):
      for i in holes:
        eta1B[p,q] += (
          0.5*GammaGamma[idx2B[(i,p)], idx2B[(i,q)]] 
          - transpose(GammaGamma)[idx2B[(i,p)], idx2B[(i,q)]]
        )

  GammaGamma = dot(Gammad, dot(occC_2B, Gammaod))
  for p in range(dim1B):
    for q in range(dim1B):
      for r in range(dim1B):
        eta1B[p,q] += (
          0.5*GammaGamma[idx2B[(r,p)], idx2B[(r,q)]] 
          + transpose(GammaGamma)[idx2B[(r,p)], idx2B[(r,q)]] 
        )


  #############################        
  # two-body flow equation  
  eta2B = np.zeros_like(Gamma)

  # 1B - 2B
  for p in range(dim1B):
    for q in range(dim1B):
      for r in range(dim1B):
        for s in range(dim1B):
          for t in range(dim1B):
            eta2B[idx2B[(p,q)],idx2B[(r,s)]] += (
              fd[p,t] * Gammaod[idx2B[(t,q)],idx2B[(r,s)]] 
              + fd[q,t] * Gammaod[idx2B[(p,t)],idx2B[(r,s)]] 
              - fd[t,r] * Gammaod[idx2B[(p,q)],idx2B[(t,s)]] 
              - fd[t,s] * Gammaod[idx2B[(p,q)],idx2B[(r,t)]]
              - fod[p,t] * Gammad[idx2B[(t,q)],idx2B[(r,s)]] 
              - fod[q,t] * Gammad[idx2B[(p,t)],idx2B[(r,s)]] 
              + fod[t,r] * Gammad[idx2B[(p,q)],idx2B[(t,s)]] 
              + fod[t,s] * Gammad[idx2B[(p,q)],idx2B[(r,t)]]
            )

  
  # 2B - 2B - particle and hole ladders
  # Gammad.occB.Gammaod
  GammaGamma = dot(Gammad, dot(occB_2B, Gammaod))

  eta2B += 0.5 * (GammaGamma - transpose(GammaGamma))

  # 2B - 2B - particle-hole chain
  
  # transform matrices to particle-hole representation and calculate 
  # Gammad_ph.occA_ph.Gammaod_ph
  Gammad_ph = ph_transform_2B(Gammad, bas2B, idx2B, basph2B, idxph2B)
  Gammaod_ph = ph_transform_2B(Gammaod, bas2B, idx2B, basph2B, idxph2B)

  GammaGamma_ph = dot(Gammad_ph, dot(occphA_2B, Gammaod_ph))

  # transform back to standard representation
  GammaGamma    = inverse_ph_transform_2B(GammaGamma_ph, bas2B, idx2B, basph2B, idxph2B)

  # commutator / antisymmetrization
  work = np.zeros_like(GammaGamma)
  for i1, (i,j) in enumerate(bas2B):
    for i2, (k,l) in enumerate(bas2B):
      work[i1, i2] -= (
        GammaGamma[i1, i2] 
        - GammaGamma[idx2B[(j,i)], i2] 
        - GammaGamma[i1, idx2B[(l,k)]] 
        + GammaGamma[idx2B[(j,i)], idx2B[(l,k)]]
      )
  GammaGamma = work

  eta2B += GammaGamma


  return eta1B, eta2B


#-----------------------------------------------------------------------------------
# derivatives 
#-----------------------------------------------------------------------------------
def flow_imsrg2(eta1B, eta2B, f, Gamma, user_data):

  dim1B     = user_data["dim1B"]
  holes     = user_data["holes"]
  particles = user_data["particles"]
  bas2B     = user_data["bas2B"]
  idx2B     = user_data["idx2B"]
  basph2B   = user_data["basph2B"]
  idxph2B   = user_data["idxph2B"]
  occB_2B   = user_data["occB_2B"]
  occC_2B   = user_data["occC_2B"]
  occphA_2B = user_data["occphA_2B"]

  #############################        
  # zero-body flow equation
  dE = 0.0

  for i in holes:
    for a in particles:
      dE += eta1B[i,a] * f[a,i] - eta1B[a,i] * f[i,a]

  for i in holes:
    for j in holes:
      for a in particles:
        for b in particles:
          dE += 0.5 * eta2B[idx2B[(i,j)], idx2B[(a,b)]] * Gamma[idx2B[(a,b)], idx2B[(i,j)]]


  #############################        
  # one-body flow equation  
  df  = np.zeros_like(f)

  # 1B - 1B
  df += commutator(eta1B, f)

  # 1B - 2B
  for p in range(dim1B):
    for q in range(dim1B):
      for i in holes:
        for a in particles:
          df[p,q] += (
            eta1B[i,a] * Gamma[idx2B[(a, p)], idx2B[(i, q)]] 
            - eta1B[a,i] * Gamma[idx2B[(i, p)], idx2B[(a, q)]] 
            - f[i,a] * eta2B[idx2B[(a, p)], idx2B[(i, q)]] 
            + f[a,i] * eta2B[idx2B[(i, p)], idx2B[(a, q)]]
          )

  # 2B - 2B
  # n_a n_b nn_c + nn_a nn_b n_c = n_a n_b + (1 - n_a - n_b) * n_c
  etaGamma = dot(eta2B, dot(occB_2B, Gamma))
  for p in range(dim1B):
    for q in range(dim1B):
      for i in holes:
        df[p,q] += 0.5*(
          etaGamma[idx2B[(i,p)], idx2B[(i,q)]] 
          + transpose(etaGamma)[idx2B[(i,p)], idx2B[(i,q)]]
        )

  etaGamma = dot(eta2B, dot(occC_2B, Gamma))
  for p in range(dim1B):
    for q in range(dim1B):
      for r in range(dim1B):
        df[p,q] += 0.5*(
          etaGamma[idx2B[(r,p)], idx2B[(r,q)]] 
          + transpose(etaGamma)[idx2B[(r,p)], idx2B[(r,q)]] 
        )


  #############################        
  # two-body flow equation  
  dGamma = np.zeros_like(Gamma)

  # 1B - 2B
  for p in range(dim1B):
    for q in range(dim1B):
      for r in range(dim1B):
        for s in range(dim1B):
          for t in range(dim1B):
            dGamma[idx2B[(p,q)],idx2B[(r,s)]] += (
              eta1B[p,t] * Gamma[idx2B[(t,q)],idx2B[(r,s)]] 
              + eta1B[q,t] * Gamma[idx2B[(p,t)],idx2B[(r,s)]] 
              - eta1B[t,r] * Gamma[idx2B[(p,q)],idx2B[(t,s)]] 
              - eta1B[t,s] * Gamma[idx2B[(p,q)],idx2B[(r,t)]]
              - f[p,t] * eta2B[idx2B[(t,q)],idx2B[(r,s)]] 
              - f[q,t] * eta2B[idx2B[(p,t)],idx2B[(r,s)]] 
              + f[t,r] * eta2B[idx2B[(p,q)],idx2B[(t,s)]] 
              + f[t,s] * eta2B[idx2B[(p,q)],idx2B[(r,t)]]
            )

  
  # 2B - 2B - particle and hole ladders
  # eta2B.occB.Gamma
  etaGamma = dot(eta2B, dot(occB_2B, Gamma))

  dGamma += 0.5 * (etaGamma + transpose(etaGamma))

  # 2B - 2B - particle-hole chain
  
  # transform matrices to particle-hole representation and calculate 
  # eta2B_ph.occA_ph.Gamma_ph
  eta2B_ph = ph_transform_2B(eta2B, bas2B, idx2B, basph2B, idxph2B)
  Gamma_ph = ph_transform_2B(Gamma, bas2B, idx2B, basph2B, idxph2B)

  etaGamma_ph = dot(eta2B_ph, dot(occphA_2B, Gamma_ph))

  # transform back to standard representation
  etaGamma    = inverse_ph_transform_2B(etaGamma_ph, bas2B, idx2B, basph2B, idxph2B)

  # commutator / antisymmetrization
  work = np.zeros_like(etaGamma)
  for i1, (i,j) in enumerate(bas2B):
    for i2, (k,l) in enumerate(bas2B):
      work[i1, i2] -= (
        etaGamma[i1, i2] 
        - etaGamma[idx2B[(j,i)], i2] 
        - etaGamma[i1, idx2B[(l,k)]] 
        + etaGamma[idx2B[(j,i)], idx2B[(l,k)]]
      )
  etaGamma = work

  dGamma += etaGamma


  return dE, df, dGamma


#-----------------------------------------------------------------------------------
# derivative wrapper
#-----------------------------------------------------------------------------------
def get_operator_from_y(y, dim1B, dim2B):
  
  # reshape the solution vector into 0B, 1B, 2B pieces
  ptr = 0
  zero_body = y[ptr]

  ptr += 1
  one_body = reshape(y[ptr:ptr+dim1B*dim1B], (dim1B, dim1B))

  ptr += dim1B*dim1B
  two_body = reshape(y[ptr:ptr+dim2B*dim2B], (dim2B, dim2B))

  return zero_body,one_body,two_body


def derivative_wrapper(t, y, user_data):

  dim1B = user_data["dim1B"]
  dim2B = dim1B*dim1B


  holes     = user_data["holes"]
  particles = user_data["particles"]
  bas1B     = user_data["bas1B"]
  bas2B     = user_data["bas2B"]
  basph2B   = user_data["basph2B"]
  idx2B     = user_data["idx2B"]
  idxph2B   = user_data["idxph2B"]
  occA_2B   = user_data["occA_2B"]
  occB_2B   = user_data["occB_2B"]
  occC_2B   = user_data["occC_2B"]
  occphA_2B = user_data["occphA_2B"]
  calc_eta  = user_data["calc_eta"]
  calc_rhs  = user_data["calc_rhs"]

  # extract operator pieces from solution vector
  E, f, Gamma = get_operator_from_y(y, dim1B, dim2B)


  # calculate the generator
  eta1B, eta2B = calc_eta(f, Gamma, user_data)

  # calculate the right-hand side
  dE, df, dGamma = calc_rhs(eta1B, eta2B, f, Gamma, user_data)

  # convert derivatives into linear array
  dy   = np.append([dE], np.append(reshape(df, -1), reshape(dGamma, -1)))

  # share data
  user_data["dE"] = dE
  user_data["eta_norm"] = np.linalg.norm(eta1B,ord='fro')+np.linalg.norm(eta2B,ord='fro')
  
  return dy

#-----------------------------------------------------------------------------------
# pairing Hamiltonian
#-----------------------------------------------------------------------------------
def pairing_hamiltonian(delta, g, user_data):
  bas1B = user_data["bas1B"]
  bas2B = user_data["bas2B"]
  idx2B = user_data["idx2B"]

  dim = len(bas1B)
  H1B = np.zeros((dim,dim))

  for i in bas1B:
    H1B[i,i] = delta*np.floor_divide(i, 2)

  dim = len(bas2B)
  H2B = np.zeros((dim, dim))

  # spin up states have even indices, spin down the next odd index
  for (i, j) in bas2B:
    if (i % 2 == 0 and j == i+1):
      for (k, l) in bas2B:
        if (k % 2 == 0 and l == k+1):
          H2B[idx2B[(i,j)],idx2B[(k,l)]] = -0.5*g
          H2B[idx2B[(j,i)],idx2B[(k,l)]] = 0.5*g
          H2B[idx2B[(i,j)],idx2B[(l,k)]] = 0.5*g
          H2B[idx2B[(j,i)],idx2B[(l,k)]] = -0.5*g
  
  return H1B, H2B

#-----------------------------------------------------------------------------------
# normal-ordered pairing Hamiltonian
#-----------------------------------------------------------------------------------
def normal_order(H1B, H2B, user_data):
  bas1B     = user_data["bas1B"]
  bas2B     = user_data["bas2B"]
  idx2B     = user_data["idx2B"]
  particles = user_data["particles"]
  holes     = user_data["holes"]

  # 0B part
  E = 0.0
  for i in holes:
    E += H1B[i,i]

  for i in holes:
    for j in holes:
      E += 0.5*H2B[idx2B[(i,j)],idx2B[(i,j)]]  

  # 1B part
  f = H1B
  for i in bas1B:
    for j in bas1B:
      for h in holes:
        f[i,j] += H2B[idx2B[(i,h)],idx2B[(j,h)]]  

  # 2B part
  Gamma = H2B

  return E, f, Gamma

#-----------------------------------------------------------------------------------
# Perturbation theory
#-----------------------------------------------------------------------------------
def calc_mbpt2(f, Gamma, user_data):
  DE2 = 0.0

  particles = user_data["particles"]
  holes     = user_data["holes"]
  idx2B     = user_data["idx2B"]

  for i in holes:
    for j in holes:
      for a in particles:
        for b in particles:
          denom = f[i,i] + f[j,j] - f[a,a] - f[b,b]
          me    = Gamma[idx2B[(a,b)],idx2B[(i,j)]]
          DE2  += 0.25*me*me/denom

  return DE2

def calc_mbpt3(f, Gamma, user_data):
  particles = user_data["particles"]
  holes     = user_data["holes"]
  idx2B     = user_data["idx2B"]

  # DE3 = 0.0

  DE3pp = 0.0
  DE3hh = 0.0
  DE3ph = 0.0

  for a in particles:
    for b in particles:
      for c in particles:
        for d in particles:
          for i in holes:
            for j in holes:
              denom = (f[i,i] + f[j,j] - f[a,a] - f[b,b])*(f[i,i] + f[j,j] - f[c,c] - f[d,d])
              me    = Gamma[idx2B[(i,j)],idx2B[(a,b)]]*Gamma[idx2B[(a,b)],idx2B[(c,d)]]*Gamma[idx2B[(c,d)],idx2B[(i,j)]]
              DE3pp += 0.125*me/denom

  for i in holes:
    for j in holes:
      for k in holes:
        for l in holes:
          for a in particles:
            for b in particles:
              denom = (f[i,i] + f[j,j] - f[a,a] - f[b,b])*(f[k,k] + f[l,l] - f[a,a] - f[b,b])
              me    = Gamma[idx2B[(a,b)],idx2B[(k,l)]]*Gamma[idx2B[(k,l)],idx2B[(i,j)]]*Gamma[idx2B[(i,j)],idx2B[(a,b)]]
              DE3hh += 0.125*me/denom

  for i in holes:
    for j in holes:
      for k in holes:
        for a in particles:
          for b in particles:
            for c in particles:
              denom = (f[i,i] + f[j,j] - f[a,a] - f[b,b])*(f[k,k] + f[j,j] - f[a,a] - f[c,c])
              me    = Gamma[idx2B[(i,j)],idx2B[(a,b)]]*Gamma[idx2B[(k,b)],idx2B[(i,c)]]*Gamma[idx2B[(a,c)],idx2B[(k,j)]]
              DE3ph -= me/denom
  return DE3pp+DE3hh+DE3ph

#------------------------------------------------------------------------------
# Main program
#------------------------------------------------------------------------------

def main():
  # g          = 0.5
  g          = float(argv[1])
  delta      = 1

  particles  = 4

  # setup shared data
  dim1B     = 8

  # this defines the reference state
  # 1st state
  # holes     = [0,1,2,3]
  # particles = [4,5,6,7]

  # 2nd state
  holes     = [0,1,4,5]
  particles = [2,3,6,7]

  # 3rd state
  # holes     = [0,1,6,7]
  # particles = [2,3,4,5]

  # basis definitions
  bas1B     = range(dim1B)
  bas2B     = construct_basis_2B(holes, particles)
  basph2B   = construct_basis_ph2B(holes, particles)

  idx2B     = construct_index_2B(bas2B)
  idxph2B   = construct_index_2B(basph2B)

  # occupation number matrices
  occ1B     = construct_occupation_1B(bas1B, holes, particles)
  occA_2B   = construct_occupationA_2B(bas2B, occ1B)
  occB_2B   = construct_occupationB_2B(bas2B, occ1B)
  occC_2B   = construct_occupationC_2B(bas2B, occ1B)

  occphA_2B = construct_occupationA_2B(basph2B, occ1B)

  # store shared data in a dictionary, so we can avoid passing the basis
  # lookups etc. as separate parameters all the time
  user_data  = {
    "dim1B":      dim1B, 
    "holes":      holes,
    "particles":  particles,
    "bas1B":      bas1B,
    "bas2B":      bas2B,
    "basph2B":    basph2B,
    "idx2B":      idx2B,
    "idxph2B":    idxph2B,
    "occ1B":      occ1B,
    "occA_2B":    occA_2B,
    "occB_2B":    occB_2B,
    "occC_2B":    occC_2B,
    "occphA_2B":  occphA_2B,

    "eta_norm":   0.0,                # variables for sharing data between ODE solver
    "dE":         0.0,                # and main routine


    "calc_eta":   eta_imtime,          # specify the generator (function object)    
    # "calc_eta":   eta_white,          # specify the generator (function object)
    "calc_rhs":   flow_imsrg2         # specify the right-hand side and truncation
  }

  # set up initial Hamiltonian
  H1B, H2B = pairing_hamiltonian(delta, g, user_data)

  E, f, Gamma = normal_order(H1B, H2B, user_data) 

  # reshape Hamiltonian into a linear array (initial ODE vector)
  y0   = np.append([E], np.append(reshape(f, -1), reshape(Gamma, -1)))

  # integrate flow equations 
  solver = ode(derivative_wrapper,jac=None)
  solver.set_integrator('vode', method='bdf', order=5, nsteps=1000)
  solver.set_f_params(user_data)
  solver.set_initial_value(y0, 0.)

  sfinal = 50
  ds = 0.1

  print "%-8s   %-14s   %-14s   %-14s   %-14s   %-14s   %-14s   %-14s   %-14s"%(
    "s", "E" , "DE(2)", "DE(3)", "E+DE", "dE/ds", 
    "||eta||", "||fod||", "||Gammaod||")
  # print "-----------------------------------------------------------------------------------------------------------------"
  print "-" * 148
  while solver.successful() and solver.t < sfinal:
    ys = solver.integrate(sfinal, step=True)
    
    dim2B = dim1B*dim1B
    E, f, Gamma = get_operator_from_y(ys, dim1B, dim2B)

    DE2 = calc_mbpt2(f, Gamma, user_data)
    DE3 = calc_mbpt3(f, Gamma, user_data)

    norm_fod     = calc_fod_norm(f, user_data)
    norm_Gammaod = calc_Gammaod_norm(Gamma, user_data)

    print "%8.5f %14.8f   %14.8f   %14.8f   %14.8f   %14.8f   %14.8f   %14.8f   %14.8f"%(
      solver.t, E , DE2, DE3, E+DE2+DE3, user_data["dE"], user_data["eta_norm"], norm_fod, norm_Gammaod)
    if abs(DE2/E) < 10e-8: break


#    solver.integrate(solver.t + ds)

#------------------------------------------------------------------------------
# make executable
#------------------------------------------------------------------------------
if __name__ == "__main__": 
  main()
