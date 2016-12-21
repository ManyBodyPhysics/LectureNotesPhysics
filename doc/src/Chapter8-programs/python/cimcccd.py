#!/usr/bin/python
from sympy import *
from pylab import *
import matplotlib.pyplot as plt


density = [0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2]
ECIMC = [6.4735287261, 7.9714751994, 9.1357756052, 10.0941195187,10.899839927, 11.5775331731, 12.1433034824, 12.6242634614, 13.0118757745 ]
ECCD = [6.468200479265218, 7.932148263480525,9.105590614303889,  10.0743311672404, 10.88453820307055,11.56523993876917,12.13635897512306,12.61239790477462,13.00438290346254]
EGF = [6.5173, 7.97249, 9.1355, 10.0959, 10.8996, 11.5752, 12.1421, 12.6146, 13.0034]
EREF = [6.987522,8.580410,9.884711, 10.980483,11.909966, 12.699944, 13.369356,13.932540,14.400847]

plt.axis([0.039,0.205,6.0,15.5])
plt.xlabel(r'Density, fm$^{-3}$', fontsize=16)
plt.ylabel(r'Energy per particle', fontsize=16)
montecarlo = plt.plot(density, ECIMC,'b-*',linewidth = 2.0, label = 'CIMC')
ccd = plt.plot(density, ECCD,'r:.', linewidth = 2.0, label = 'CCD')
green = plt.plot(density, EGF,'r:.', linewidth = 2.0, label = 'GF')
refenergy = plt.plot(density, EREF, 'm:v',linewidth = 2.0, label = 'Reference energy')
plt.legend()
plt.savefig('cimcccd.pdf', format='pdf')
plt.show()



