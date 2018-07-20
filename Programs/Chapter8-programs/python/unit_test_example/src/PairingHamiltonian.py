from sympy import *
#import PairingBasisHardCoded as basis
import PairingBasisGenerated as pbg
basis = pbg.PairingBasisGen(8,4)

g = Symbol('g')

def h0(p,q):
    if p == q:
        p1, s1 = basis.states[p]
        return (p1 - 1)
    else:
        return 0

def f(p,q):
    if p == q:
        return 0
    s = h0(p,q)
    for i in basis.below_fermi:
        s += assym(p,i,q,i)
        return s


def assym(p,q,r,s):
    p1, s1 = basis.states[p]
    p2, s2 = basis.states[q]
    p3, s3 = basis.states[r]
    p4, s4 = basis.states[s]

    if p1 != p2 or p3 != p4:
        return 0
    if s1 == s2 or s3 == s4:
        return 0
    if s1 == s3 and s2 == s4:
        return -g/2.
    if s1 == s4 and s2 == s3:
        return g/2.

