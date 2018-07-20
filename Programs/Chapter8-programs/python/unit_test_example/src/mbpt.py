#!/usr/bin/python
from sympy import *
from pylab import *
import matplotlib.pyplot as plt

import PairingHamiltonian as ph
import Energies
from WhoKnowsWhat import *

import PairingBasisGenerated as pb

#print dog

correlation_energy_mbpt2 = Energies.corr_mbpt2()
correlation_energy_mbpt3 = Energies.corr_mbpt3()

ga = linspace(-1,1,20)
e1 = []
corr2 = []
corr3 = []

for g_val in ga:
    H1 = matrix([[2-g_val , -g_val/2.,  -g_val/2., -g_val/2., -g_val/2.,     0],
                         [-g_val/2.,   4-g_val,  -g_val/2., -g_val/2.,    0., -g_val/2.],
                         [-g_val/2., -g_val/2.,    6-g_val,     0, -g_val/2., -g_val/2.],
                 [-g_val/2., -g_val/2.,      0,   6-g_val, -g_val/2., -g_val/2.],
                 [-g_val/2.,     0,  -g_val/2., -g_val/2.,   8-g_val, -g_val/2.],
                 [0    , -g_val/2.,  -g_val/2., -g_val/2., -g_val/2.,  10-g_val]])

    u1, v1 = linalg.eig(H1)
    e1.append(min(u1))
    corr2.append((correlation_energy_mbpt2).subs(ph.g,g_val))
    corr3.append((correlation_energy_mbpt3).subs(ph.g,g_val))

exact = e1 - (2-ga)

plt.axis([-1,1,-0.5,0.05])
plt.xlabel(r'Interaction strength, $g$', fontsize=16)
plt.ylabel(r'Correlation energy', fontsize=16)
exact = plt.plot(ga, exact,'b-*',linewidth = 2.0, label = 'Exact')
mbpt2 = plt.plot(ga, corr2,'r:.', linewidth = 2.0, label = 'MBPT2')
mbpt3 = plt.plot(ga, corr3, 'm:v',linewidth = 2.0, label = 'MBPT3')
plt.legend()
plt.savefig('perturbationtheory.pdf', format='pdf')
plt.show()



    
