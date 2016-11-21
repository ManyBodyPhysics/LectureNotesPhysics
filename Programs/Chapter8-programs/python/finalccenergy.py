#!/usr/bin/python
from sympy import *
from pylab import *
import matplotlib.pyplot as plt


density = [0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2]
ECCD = [6.468200479265218, 7.932148263480525,9.105590614303889,  10.0743311672404, 10.88453820307055,11.56523993876917,12.13635897512306,12.61239790477462,13.00438290346254]
EREF = [6.987522,8.580410,9.884711, 10.980483,11.909966, 12.699944, 13.369356,13.932540,14.400847]
plt.axis([0.039,0.205,6.0,15.5])
plt.xlabel(r'Density $\rho$, fm$^{-3}$', fontsize=16)
plt.ylabel(r'Energy per particle', fontsize=16)
ccd = plt.plot(density, ECCD,'r:.', linewidth = 2.0, label = 'CCD')
refenergy = plt.plot(density, EREF, 'm:v',linewidth = 2.0, label = 'Reference energy')
plt.legend()
plt.savefig('finalccenergy.pdf', format='pdf')
plt.show()


