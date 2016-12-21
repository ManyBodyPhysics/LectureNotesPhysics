#!/usr/bin/python
from sympy import *
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
# Load in data file
data = np.loadtxt("data.txt")
ga = data[:,0]
exact = data[:,2]
corrCCD = data[:,8]
corr4 = data[:,7]


plt.axis([-1,1,-0.5,0.05])
plt.xlabel(r'Interaction strength, $g$', fontsize=16)
plt.ylabel(r'Correlation energy', fontsize=16)
exact = plt.plot(ga, exact,'b-*',linewidth = 2.0, label = 'Exact')
mbpt4 = plt.plot(ga, corr4,'r:.', linewidth = 2.0, label = 'MBPT4')
ccd = plt.plot(ga, corrCCD, 'm:v',linewidth = 2.0, label = 'CCD')
plt.legend()
plt.savefig('CCDMBPT4theory.pdf', format='pdf')
plt.show()



