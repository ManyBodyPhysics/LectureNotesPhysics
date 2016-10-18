#!/usr/bin/env python

#------------------------------------------------------------------------------
# imsrg_pairing.py
#
# author:   H. Hergert 
# version:  1.0.1
# date:     Oct 12, 2016
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
from scipy.integrate import odeint


#-----------------------------------------------------------------------------------
# commutator of matrices
#-----------------------------------------------------------------------------------
def commutator(a,b):
  return dot(a,b) - dot(b,a)

#-----------------------------------------------------------------------------------
# basis and index functions
#-----------------------------------------------------------------------------------

def construct_basis_2B(holes, particles):
  basis = []
  for i in holes:
    for j in holes:
      basis.append((i, j))

  for i in holes:
    for j in particles:
      basis.append((i, j))

  for i in particles:
    for j in holes:
      basis.append((i, j))

  for i in particles:
    for j in particles:
      basis.append((i, j))

  return basis


def construct_basis_3B(holes, particles):
  basis = []
  for i in holes:
    for j in holes:
      for k in holes:
        basis.append((i, j, k))

  for i in holes:
    for j in holes:
      for k in particles:
        basis.append((i, j, k))

  for i in holes:
    for j in particles:
      for k in holes:
        basis.append((i, j, k))

  for i in particles:
    for j in holes:
      for k in holes:
        basis.append((i, j, k))

  for i in holes:
    for j in particles:
      for k in particles:
        basis.append((i, j, k))

  for i in particles:
    for j in holes:
      for k in particles:
        basis.append((i, j, k))

  for i in particles:
    for j in particles:
      for k in holes:
        basis.append((i, j, k))

  for i in particles:
    for j in particles:
      for k in particles:
        basis.append((i, j, k))

  return basis


def construct_basis_ph2B(holes, particles):
  basis = []
  for i in holes:
    for j in holes:
      basis.append((i, j))

  for i in holes:
    for j in particles:
      basis.append((i, j))

  for i in particles:
    for j in holes:
      basis.append((i, j))

  for i in particles:
    for j in particles:
      basis.append((i, j))

  return basis


#
# We use dictionaries for the reverse lookup of state indices
#
def construct_index_2B(bas2B):
  index = { }
  for i, state in enumerate(bas2B):
    index[state] = i

  return index

def construct_index_3B(bas3B):
  index = { }
  for i, state in enumerate(bas3B):
    index[state] = i

  return index


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
# transform matrices to particle-hole representation
#-----------------------------------------------------------------------------------
def ph_transform_2B(Gamma, bas2B, idx2B, basph2B, idxph2B):
  dim = len(basph2B)
  Gamma_ph = np.zeros((dim, dim))

  for i1, (a,b) in enumerate(basph2B):
    for i2, (c, d) in enumerate(basph2B):
      Gamma_ph[i1, i2] += -1.0 * Gamma[idx2B[(a,d)], idx2B[(c,b)]]

  return Gamma_ph

def inverse_ph_transform_2B(Gamma_ph, bas2B, idx2B, basph2B, idxph2B):
  dim = len(bas2B)
  Gamma = np.zeros((dim, dim))

  for i1, (a,b) in enumerate(bas2B):
    for i2, (c, d) in enumerate(bas2B):
      Gamma[i1, i2]                     -= Gamma_ph[idxph2B[(a,d)], idxph2B[(c,b)]]
  
  return Gamma



#-----------------------------------------------------------------------------------
# generators
#-----------------------------------------------------------------------------------
def eta_brillouin(f, Gamma, user_data):
  dim1B     = user_data["dim"]
  particles = user_data["particles"]
  holes     = user_data["holes"]
  idx2B     = user_data["idx2B"]

  # one-body part of the generator
  eta1B  = np.zeros_like(f)

  for i in particles:
    for j in holes:
      # (1-n_i)n_j - n_i(1-n_j) = n_j - n_i
      eta1B[i, j] =  f[i,j]
      eta1B[j, i] = -f[i,j]

  # two-body part of the generator
  eta2B = np.zeros_like(Gamma)

  for i in particles:
    for j in particles:
      for k in holes:
        for l in holes:
          val = Gamma[idx2B[(i,j)], idx2B[(k,l)]]

          eta2B[idx2B[(i,j)],idx2B[(k,l)]] = val
          eta2B[idx2B[(k,l)],idx2B[(i,j)]] = -val

  return eta1B, eta2B


def eta_white(f, Gamma, user_data):
  dim1B     = user_data["dim"]
  particles = user_data["particles"]
  holes     = user_data["holes"]
  idx2B     = user_data["idx2B"]

  # one-body part of the generator
  eta1B  = np.zeros_like(f)

  for i in particles:
    for j in holes:
      denom = f[i,i] - f[j,j] + Gamma[idx2B[(i,j)], idx2B[(i,j)]]
      val = f[i,j]/denom
      eta1B[i, j] =  val
      eta1B[j, i] = -val 

  # two-body part of the generator
  eta2B = np.zeros_like(Gamma)

  for i in particles:
    for j in particles:
      for k in holes:
        for l in holes:
          denom = ( 
            f[i,i] + f[j,j] - f[k,k] - f[l,l] 
            + Gamma[idx2B[(i,j)],idx2B[(i,j)]] 
            + Gamma[idx2B[(k,l)],idx2B[(k,l)]]
            - Gamma[idx2B[(i,k)],idx2B[(i,k)]] 
            - Gamma[idx2B[(i,l)],idx2B[(i,l)]] 
            - Gamma[idx2B[(j,k)],idx2B[(j,k)]] 
            - Gamma[idx2B[(j,l)],idx2B[(j,l)]] 
          )

          val = Gamma[idx2B[(i,j)], idx2B[(k,l)]] / denom

          eta2B[idx2B[(i,j)],idx2B[(k,l)]] = val
          eta2B[idx2B[(k,l)],idx2B[(i,j)]] = -val

  return eta1B, eta2B


def eta_white_atan(f, Gamma, user_data):
  dim1B     = user_data["dim"]
  particles = user_data["particles"]
  holes     = user_data["holes"]
  idx2B     = user_data["idx2B"]

  # one-body part of the generator
  eta1B  = np.zeros_like(f)

  for i in particles:
    for j in holes:
      denom = f[i,i] - f[j,j] + Gamma[idx2B[(i,j)], idx2B[(i,j)]]
      val = 0.5 * np.arctan(2 * f[i,j]/denom)
      eta1B[i, j] =  val
      eta1B[j, i] = -val 

  # two-body part of the generator
  eta2B = np.zeros_like(Gamma)

  for i in particles:
    for j in particles:
      for k in holes:
        for l in holes:
          denom = ( 
            f[i,i] + f[j,j] - f[k,k] - f[l,l] 
            + Gamma[idx2B[(i,j)],idx2B[(i,j)]] 
            + Gamma[idx2B[(k,l)],idx2B[(k,l)]] 
            - Gamma[idx2B[(i,k)],idx2B[(i,k)]] 
            - Gamma[idx2B[(i,l)],idx2B[(i,l)]] 
            - Gamma[idx2B[(j,k)],idx2B[(j,k)]] 
            - Gamma[idx2B[(j,l)],idx2B[(j,l)]] 
          )

          val = 0.5 * np.arctan(2 * Gamma[idx2B[(i,j)], idx2B[(k,l)]] / denom)

          eta2B[idx2B[(i,j)],idx2B[(k,l)]] = val
          eta2B[idx2B[(k,l)],idx2B[(i,j)]] = -val

  return eta1B, eta2B

# def eta1B_brillouin(f, Gamma, idx2B, occ1B, occB_2B, occC_2B, holes, particles):
#   dim = f.shape[0]
#   eta1B  = np.zeros_like(f)

#   for i in particles:
#     for j in holes:
#       # (1-n_i)n_j - n_i(1-n_j) = n_j - n_i
#       eta1B[i, j] =  f[i,j]
#       eta1B[j, i] = -f[i,j]

#   return eta1B

# def eta2B_brillouin(f, Gamma, bas2B, idx2B, holes, particles):
#   # dim = Gamma.shape[0]
#   eta2B = np.zeros_like(Gamma)

#   for i in particles:
#     for j in particles:
#       for k in holes:
#         for l in holes:
#           val = Gamma[idx2B[(i,j)], idx2B[(k,l)]]

#           eta2B[idx2B[(i,j)],idx2B[(k,l)]] = val
#           eta2B[idx2B[(k,l)],idx2B[(i,j)]] = -val

#   return eta2B


# def eta1B_white_atan(f, Gamma, idx2B, occ1B, occB_2B, occC_2B, holes, particles):
#   dim = f.shape[0]
#   eta1B  = np.zeros_like(f)

#   for i in particles:
#     for j in holes:
#       denom = f[i,i] - f[j,j] + Gamma[idx2B[(i,j)], idx2B[(i,j)]]
#       val = 0.5 * np.arctan(2 * f[i,j]/denom)
#       eta1B[i, j] =  val
#       eta1B[j, i] = -val 

#   return eta1B

# def eta2B_white_atan(f, Gamma, bas2B, idx2B, holes, particles):
#   # dim = Gamma.shape[0]
#   eta2B = np.zeros_like(Gamma)

#   for i in particles:
#     for j in particles:
#       for k in holes:
#         for l in holes:
#           denom = ( 
#             f[i,i] + f[j,j] - f[k,k] - f[l,l] 
#             + Gamma[idx2B[(i,j)],idx2B[(i,j)]] 
#             + Gamma[idx2B[(k,l)],idx2B[(k,l)]] 
#             - Gamma[idx2B[(i,k)],idx2B[(i,k)]] 
#             - Gamma[idx2B[(i,l)],idx2B[(i,l)]] 
#             - Gamma[idx2B[(j,k)],idx2B[(j,k)]] 
#             - Gamma[idx2B[(j,l)],idx2B[(j,l)]] 
#           )

#           val = 0.5 * np.arctan(2 * Gamma[idx2B[(i,j)], idx2B[(k,l)]] / denom)

#           eta2B[idx2B[(i,j)],idx2B[(k,l)]] = val
#           eta2B[idx2B[(k,l)],idx2B[(i,j)]] = -val

#   return eta2B


# def eta1B_white(f, Gamma, idx2B, occ1B, occB_2B, occC_2B, holes, particles):
#   dim = f.shape[0]
#   eta1B  = np.zeros_like(f)


#   for i in particles:
#     for j in holes:
#       denom = f[i,i] - f[j,j] + Gamma[idx2B[(i,j)], idx2B[(i,j)]]
#       val = f[i,j]/denom
#       eta1B[i, j] =  val
#       eta1B[j, i] = -val 

#   return eta1B

# def eta2B_white(f, Gamma, bas2B, idx2B, holes, particles):
#   # dim = Gamma.shape[0]
#   eta2B = np.zeros_like(Gamma)

#   for i in particles:
#     for j in particles:
#       for k in holes:
#         for l in holes:
#           denom = ( 
#             f[i,i] + f[j,j] - f[k,k] - f[l,l] 
#             + Gamma[idx2B[(i,j)],idx2B[(i,j)]] 
#             + Gamma[idx2B[(k,l)],idx2B[(k,l)]]
#             - Gamma[idx2B[(i,k)],idx2B[(i,k)]] 
#             - Gamma[idx2B[(i,l)],idx2B[(i,l)]] 
#             - Gamma[idx2B[(j,k)],idx2B[(j,k)]] 
#             - Gamma[idx2B[(j,l)],idx2B[(j,l)]] 
#           )

#           val = Gamma[idx2B[(i,j)], idx2B[(k,l)]] / denom

#           eta2B[idx2B[(i,j)],idx2B[(k,l)]] = val
#           eta2B[idx2B[(k,l)],idx2B[(i,j)]] = -val

#   return eta2B

#-----------------------------------------------------------------------------------
# derivatives 
#-----------------------------------------------------------------------------------
def imsrg2_dE(eta1B, eta2B, f, Gamma, idx2B, holes, particles):
  dE = 0.0

  for a in holes:
    for b in particles:
      dE += eta1B[a,b] * f[b, a] - eta1B[b,a] * f[a, b]

  for a in holes:
    for b in holes:
      for c in particles:
        for d in particles:
          dE += 0.5 * eta2B[idx2B[(a,b)], idx2B[(c,d)]] * Gamma[idx2B[(c,d)], idx2B[(a,b)]]

  return dE

def imsrg2_df(eta1B, eta2B, f, Gamma, idx2B, occB_2B, occC_2B, holes, particles):
  dim = f.shape[0]
  df  = np.zeros_like(f)

  # 1B - 1B
  df += commutator(eta1B, f)

  # 1B - 2B
  for i in range(dim):
    for j in range(dim):
      for a in holes:
        for b in particles:
          df[i,j] += (
            eta1B[a,b] * Gamma[idx2B[(b, i)], idx2B[(a, j)]] 
            - eta1B[b,a] * Gamma[idx2B[(a, i)], idx2B[(b, j)]] 
            - f[a,b] * eta2B[idx2B[(b, i)], idx2B[(a, j)]] 
            + f[b,a] * eta2B[idx2B[(a, i)], idx2B[(b, j)]]
          )

  # 2B - 2B
  # n_a n_b nn_c + nn_a nn_b n_c = n_a n_b + (1 - n_a - n_b) * n_c
  etaGamma = dot(eta2B, dot(occB_2B, Gamma))
  for i in range(dim):
    for j in range(dim):
      for c in holes:
        df[i,j] += (
          0.5*etaGamma[idx2B[(c,i)], idx2B[(c,j)]] 
          + transpose(etaGamma)[idx2B[(c,i)], idx2B[(c,j)]]
        )

  etaGamma = dot(eta2B, dot(occC_2B, Gamma))
  for i in range(dim):
    for j in range(dim):
      for c in range(dim):
        df[i,j] += (
          0.5*etaGamma[idx2B[(c,i)], idx2B[(c,j)]] 
          + transpose(etaGamma)[idx2B[(c,i)], idx2B[(c,j)]] 
        )
  
  return df



def imsrg2_dGamma(eta1B, eta2B, f, Gamma, bas2B, idx2B, basph2B, idxph2B, occB_2B, occphA_2B, occ1B):
  dim    = f.shape[0]
  dGamma = np.zeros_like(Gamma)

  # 1B - 2B
  for i in range(dim):
    for j in range(dim):
      for k in range(dim):
        for l in range(dim):
          for a in range(dim):
            dGamma[idx2B[(i,j)],idx2B[(k,l)]] += (
              eta1B[i,a] * Gamma[idx2B[(a,j)],idx2B[(k,l)]] 
              + eta1B[j,a] * Gamma[idx2B[(i,a)],idx2B[(k,l)]] 
              - eta1B[a,k] * Gamma[idx2B[(i,j)],idx2B[(a,l)]] 
              - eta1B[a,l] * Gamma[idx2B[(i,j)],idx2B[(k,a)]]
              - f[i,a] * eta2B[idx2B[(a,j)],idx2B[(k,l)]] 
              - f[j,a] * eta2B[idx2B[(i,a)],idx2B[(k,l)]] 
              + f[a,k] * eta2B[idx2B[(i,j)],idx2B[(a,l)]] 
              + f[a,l] * eta2B[idx2B[(i,j)],idx2B[(k,a)]]
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


  # etaGamma2 = np.zeros_like(dGamma)
  # for i1, (i, j) in enumerate(bas2B):
  #   for i2, (k, l) in enumerate(bas2B):
  #     # term1 = 0.
  #     # term2 = 0.
  #     # term3 = 0.
  #     # term4 = 0.
  #     for a in range(dim):
  #       for b in range(dim):
  #         # etaGamma2[i1,i2] += occA_2B[idx2B[(a,b)],idx2B[(a,b)]] * (
  #         # eta2B[idx2B[(i,a)], idx2B[(k,b)]] * Gamma[idx2B[(b,j)], idx2B[(a,l)]]
  #         # -  Gamma[idx2B[(i,a)], idx2B[(k,b)]] * eta2B[idx2B[(b,j)], idx2B[(a,l)]]
  #         # -  eta2B[idx2B[(j,a)], idx2B[(k,b)]] * Gamma[idx2B[(b,i)], idx2B[(a,l)]]
  #         # +  Gamma[idx2B[(j,a)], idx2B[(k,b)]] * eta2B[idx2B[(b,i)], idx2B[(a,l)]])
  #         # term1 += (occ1B[a] - occ1B[b]) * eta2B[idx2B[(i,a)], idx2B[(k,b)]] * Gamma[idx2B[(b,j)], idx2B[(a,l)]]
  #         # term2 -= (occ1B[a] - occ1B[b]) * Gamma[idx2B[(i,a)], idx2B[(k,b)]] * eta2B[idx2B[(b,j)], idx2B[(a,l)]]
  #         # term3 -= (occ1B[a] - occ1B[b]) * eta2B[idx2B[(j,a)], idx2B[(k,b)]] * Gamma[idx2B[(b,i)], idx2B[(a,l)]]
  #         # term4 += (occ1B[a] - occ1B[b]) * Gamma[idx2B[(j,a)], idx2B[(k,b)]] * Gamma[idx2B[(b,i)], idx2B[(a,l)]]

  #         etaGamma2[i1,i2] += (occ1B[a] - occ1B[b]) * (
  #         eta2B[idx2B[(i,a)], idx2B[(k,b)]] * Gamma[idx2B[(b,j)], idx2B[(a,l)]]
  #         -  Gamma[idx2B[(i,a)], idx2B[(k,b)]] * eta2B[idx2B[(b,j)], idx2B[(a,l)]]
  #         -  eta2B[idx2B[(j,a)], idx2B[(k,b)]] * Gamma[idx2B[(b,i)], idx2B[(a,l)]]
  #         +  Gamma[idx2B[(j,a)], idx2B[(k,b)]] * eta2B[idx2B[(b,i)], idx2B[(a,l)]])

  #         # print i1, i2, (i,j), (k,l), a, b, etaGamma2[i1,i2]
  #         # etaGamma2[i1,i2] += (
  #         # eta2B[idx2B[(i,a)], idx2B[(k,b)]] * Gamma[idx2B[(b,j)], idx2B[(a,l)]]
  #         # -  Gamma[idx2B[(i,a)], idx2B[(k,b)]] * eta2B[idx2B[(b,j)], idx2B[(a,l)]]
  #         # -  eta2B[idx2B[(j,a)], idx2B[(k,b)]] * Gamma[idx2B[(b,i)], idx2B[(a,l)]]
  #         # +  Gamma[idx2B[(j,a)], idx2B[(k,b)]] * eta2B[idx2B[(b,i)], idx2B[(a,l)]])
  
  # etaGamma3_ph = np.zeros_like(dGamma)
  # for i1, (i, j) in enumerate(basph2B):
  #   for i2, (k, l) in enumerate(basph2B):
  #     for a in range(dim):
  #       for b in range(dim):
  #         # etaGamma3_ph[i1,i2] += (occ1B[a] - occ1B[b]) * (
  #         # eta2B_ph[idxph2B[(i,j)], idxph2B[(a,b)]] * Gamma_ph[idxph2B[(a,b)], idxph2B[(k,l)]]
  #         # -  Gamma_ph[idxph2B[(i,j)], idxph2B[(a,b)]] * eta2B_ph[idxph2B[(a,b)], idxph2B[(k,l)]]
  #         # -  eta2B_ph[idxph2B[(l,j)], idxph2B[(a,b)]] * Gamma_ph[idxph2B[(a,b)], idxph2B[(k,i)]]
  #         # +  Gamma_ph[idxph2B[(l,j)], idxph2B[(a,b)]] * eta2B_ph[idxph2B[(a,b)], idxph2B[(k,i)]])
  #         etaGamma3_ph[i1,i2] += (occ1B[a] - occ1B[b]) * (
  #         eta2B_ph[idxph2B[(i,j)], idxph2B[(a,b)]] * Gamma_ph[idxph2B[(a,b)], idxph2B[(k,l)]]
  #         -  Gamma_ph[idxph2B[(i,j)], idxph2B[(a,b)]] * eta2B_ph[idxph2B[(a,b)], idxph2B[(k,l)]])

  #         # print i1, i2, (i,j), (k,l), a, b, etaGamma3_ph[i1,i2]

  #         # etaGamma2[i1,i2] += (
  #         # eta2B[idx2B[(i,a)], idx2B[(k,b)]] * Gamma[idx2B[(b,j)], idx2B[(a,l)]]
  #         # -  Gamma[idx2B[(i,a)], idx2B[(k,b)]] * eta2B[idx2B[(b,j)], idx2B[(a,l)]]
  #         # -  eta2B[idx2B[(j,a)], idx2B[(k,b)]] * Gamma[idx2B[(b,i)], idx2B[(a,l)]]
  #         # +  Gamma[idx2B[(j,a)], idx2B[(k,b)]] * eta2B[idx2B[(b,i)], idx2B[(a,l)]])

  # etaGamma3 = inverse_ph_transform_2B(etaGamma3_ph, bas2B, idx2B, basph2B, idxph2B)
  # print np.linalg.norm(etaGamma3-transpose(etaGamma3),ord='fro')
  # work = np.zeros_like(etaGamma)
  # for i1, (i,j) in enumerate(bas2B):
  #   for i2, (k,l) in enumerate(bas2B):
  #     work[i1, i2] -= 0.5*(etaGamma3[i1, i2] - etaGamma3[idx2B[(j,i)], i2] 
  #     - etaGamma3[i1, idx2B[(l,k)]] - etaGamma3[idx2B[(j,i)], idx2B[(l,k)]])
  # etaGamma3 = work  

  work = np.zeros_like(etaGamma)
  for i1, (i,j) in enumerate(bas2B):
    for i2, (k,l) in enumerate(bas2B):
      work[i1, i2] -= (
        etaGamma[i1, i2] 
        - etaGamma[idx2B[(j,i)], i2] 
        - etaGamma[i1, idx2B[(l,k)]] 
        + etaGamma[idx2B[(j,i)], idx2B[(l,k)]]
      )
      # print "%10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f"%(
      # etaGamma[i1, i2], etaGamma[idx2B[(j,i)], i2], etaGamma[i1, idx2B[(l,k)]], etaGamma[idx2B[(j,i)], idx2B[(l,k)]], 
      # work[i1, i2], work[idx2B[(j,i)], i2], work[i1, idx2B[(l,k)]], work[idx2B[(j,i)], idx2B[(l,k)]],
      # etaGamma2[i1, i2], etaGamma2[idx2B[(j,i)], i2], etaGamma2[i1, idx2B[(l,k)]], etaGamma2[idx2B[(j,i)], idx2B[(l,k)]]
      # )
  etaGamma = work

  # for i1, (i,j) in enumerate(bas2B):
  #   for i2, (k,l) in enumerate(bas2B):
  #     if etaGamma[i1,i2] != 0.0 or etaGamma2[i1,i2] != 0.0 or etaGamma3[i1,i2] != 0.0:
  #       print i1, i2, (i,j), (k,l)
  #       print "%10.8f %10.8f %10.8f %10.8f | %10.8f %10.8f %10.8f %10.8f | %10.8f %10.8f %10.8f %10.8f"%(
  #       etaGamma[i1, i2], etaGamma[idx2B[(j,i)], i2], etaGamma[i1, idx2B[(l,k)]], etaGamma[idx2B[(j,i)], idx2B[(l,k)]],
  #       etaGamma2[i1, i2], etaGamma2[idx2B[(j,i)], i2], etaGamma2[i1, idx2B[(l,k)]], etaGamma2[idx2B[(j,i)], idx2B[(l,k)]],
  #       etaGamma3[i1, i2], etaGamma3[idx2B[(j,i)], i2], etaGamma3[i1, idx2B[(l,k)]], etaGamma3[idx2B[(j,i)], idx2B[(l,k)]]
  #       )
  # for i1, (i,j) in enumerate(bas2B):
  #   for i2, (k,l) in enumerate(bas2B):
  #     if etaGamma[i1,i2] != etaGamma2[i1,i2]:
  #       print i1, i2, (i,j), (k,l)
  #       print "%10.8f %10.8f %10.8f %10.8f | %10.8f %10.8f %10.8f %10.8f"%(
  #       etaGamma[i1, i2], etaGamma[idx2B[(j,i)], i2], etaGamma[i1, idx2B[(l,k)]], etaGamma[idx2B[(j,i)], idx2B[(l,k)]],
  #       etaGamma2[i1, i2], etaGamma2[idx2B[(j,i)], i2], etaGamma2[i1, idx2B[(l,k)]], etaGamma2[idx2B[(j,i)], idx2B[(l,k)]]
  #       )

  # print np.linalg.norm(etaGamma, ord='fro')
  # print np.linalg.norm(etaGamma2, ord='fro')
  # print np.linalg.norm(etaGamma3, ord='fro')
  # for i1, (i,j) in enumerate(bas2B):
  #   for i2, (k,l) in enumerate(bas2B):
  #     if etaGamma[i1,i2] != etaGamma2[i1,i2] or etaGamma2[i1,i2] != etaGamma3[i1,i2] or etaGamma[i1,i2] != etaGamma3[i1,i2] :
  #       print i1, i2, (i,j), (k,l), 
  #       print "%10.8f %10.8f %10.8f "%(
  #       etaGamma[i1,i2], etaGamma2[i1,i2], etaGamma3[i1,i2]      
  #       )

  # for i1, (i,j) in enumerate(bas2B):
  #   for i2, (k,l) in enumerate(bas2B):
  #     if work[i1,i2]  != etaGamma2[i1,i2]:
  #       print i1, i2, (i,j), (k,l), work[i1,i2], etaGamma2[i1,i2]
  #       print "%10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f %10.8f"%(      
  #       work[i1, i2], work[idx2B[(j,i)], i2], work[i1, idx2B[(l,k)]], work[idx2B[(j,i)], idx2B[(l,k)]],
  #       etaGamma2[i1, i2], etaGamma2[idx2B[(j,i)], i2], etaGamma2[i1, idx2B[(l,k)]], etaGamma2[idx2B[(j,i)], idx2B[(l,k)]]
  #       )

  # print np.linalg.norm(work - transpose(work), ord='fro')
  # print np.linalg.norm(etaGamma2 - transpose(etaGamma2), ord='fro')

  # antisymmetrize indices and add to derivative
  # for i1, (i,j) in enumerate(bas2B):
  #   for i2, (k,l) in enumerate(bas2B):
  #     # dGamma[i1, i2] += etaGamma2[i1, i2]
  #     dGamma[i1, i2] += etaGamma[i1, i2]



  # add to derivative
  # dGamma += etaGamma

  # for i1, (i,j) in enumerate(bas2B):
  #   for i2, (k,l) in enumerate(bas2B):
  #     if dGamma[i1,i2] != 0.0:
  #       print i1, i2, (i,j), (k,l)
  #       print "%10.8f %10.8f %10.8f %10.8f"%(      
  #       dGamma[i1, i2], dGamma[idx2B[(j,i)], i2], dGamma[i1, idx2B[(l,k)]], dGamma[idx2B[(j,i)], idx2B[(l,k)]]
  #       )

  return dGamma


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


def derivative_wrapper(y, t, user_data):
  dim1B = user_data["dim"]
  dim2B = dim*dim


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

  # extract operator pieces from solution vector
  E, f, Gamma = get_operator_from_y(y, dim, dim2B)


  # calculate the generator
  # eta1B = eta1B_white(f, Gamma, idx2B, occ1B, occB_2B, occC_2B, holes, particles)
  # eta2B = eta2B_white(f, Gamma, bas2B, idx2B, holes, particles)
  # eta1B = eta1B_brillouin(f, Gamma, idx2B, occ1B, occB_2B, occC_2B, holes, particles)
  # eta2B = eta2B_brillouin(f, Gamma, bas2B, idx2B, holes, particles)
  eta1B, eta2B = calc_eta(f, Gamma, user_data)

  # calculate the right-hand side
  dE     = imsrg2_dE(eta1B, eta2B, f, Gamma, idx2B, holes, particles)
  df     = imsrg2_df(eta1B, eta2B, f, Gamma, idx2B, occB_2B, occC_2B, holes, particles)
  dGamma = imsrg2_dGamma(eta1B, eta2B, f, Gamma, bas2B, idx2B, basph2B, idxph2B, occB_2B, occphA_2B, occ1B)


  # convert derivatives into linear array
  dydt   = np.append([dE], np.append(reshape(df, -1), reshape(dGamma, -1)))

  # print flow parameter etc. 
  print "%7.5f   %10.8f   %10.8f"%(t, E , np.linalg.norm(eta1B,ord='fro')+np.linalg.norm(eta2B,ord='fro'))
  
  return dydt

#-----------------------------------------------------------------------------------
# initialize normal-ordered pairing Hamiltonian
#-----------------------------------------------------------------------------------
#
# 0B
#
def E0(delta, g):
  return 2*delta - g

#
# 1B
#
def f0(delta, g, bas1B, holes, particles):
  dim = len(bas1B)
  f = np.zeros((dim, dim))

  for i in bas1B:
    if i in holes:
      f[i, i] = delta*np.floor_divide(i, 2) - 0.5*g
    elif i in particles:
      f[i, i] = delta*np.floor_divide(i, 2)

  return f

#
# 2B
#
def Gamma0(delta, g, bas2B, idx2B):
  dim = len(bas2B)
  Gamma = np.zeros((dim, dim))

  # spin up states have even indices, spin down the next odd index
  for (i, j) in bas2B:
    if (i % 2 == 0 and j == i+1):
      for (k, l) in bas2B:
        if (k % 2 == 0 and l == k+1):
          Gamma[idx2B[(i,j)],idx2B[(k,l)]] = -0.5*g
          Gamma[idx2B[(j,i)],idx2B[(k,l)]] = 0.5*g
          Gamma[idx2B[(i,j)],idx2B[(l,k)]] = 0.5*g
          Gamma[idx2B[(j,i)],idx2B[(l,k)]] = -0.5*g
  
  return Gamma

#
# 3B
#
def W0(delta, g, bas3B, idx3B):
  dim = len(bas3B)
  W   = np.zeros((dim, dim))
  return W


#------------------------------------------------------------------------------
# Main program
#------------------------------------------------------------------------------


g          = 0.5
delta      = 1
particles  = 4

# setup shared data
dim        = 8

holes     = range(particles)
particles = range(particles,dim)

bas1B     = range(dim)
bas2B     = construct_basis_2B(holes, particles)
bas3B     = construct_basis_3B(holes, particles)
basph2B   = construct_basis_ph2B(holes, particles)

idx2B     = construct_index_2B(bas2B)
idxph2B   = construct_index_2B(basph2B)

occ1B     = construct_occupation_1B(bas1B, holes, particles)
occA_2B   = construct_occupationA_2B(bas2B, occ1B)
occB_2B   = construct_occupationB_2B(bas2B, occ1B)
occC_2B   = construct_occupationC_2B(bas2B, occ1B)

occphA_2B = construct_occupationA_2B(basph2B, occ1B)

user_data  = {
  "dim":        dim, 
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
  "calc_eta":   eta_white_atan
}

# set up initial Hamiltonian
E     = E0(delta, g)
f     = f0(delta, g, bas1B, holes, particles)
Gamma = Gamma0(delta, g, bas2B, idx2B)


# reshape Hamiltonian into a linear array (initial ODE vector)
y0   = np.append([E], np.append(reshape(f, -1), reshape(Gamma, -1)))

# integrate flow equations 
ys  = odeint(derivative_wrapper, y0, [0, 7], args=(user_data,))


