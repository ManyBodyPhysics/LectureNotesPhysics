#!/usr/bin/env python
from sympy import *
from pylab import *
import matplotlib.pyplot as plt

fig = figure(figsize=(9,6))
ax  = fig.gca()

density = [0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2]
ECIMC = [6.4735287261, 7.9714751994, 9.1357756052, 10.0941195187,10.899839927, 11.5775331731, 12.1433034824, 12.6242634614, 13.0118757745 ]
ECCD = [6.468200479265218, 7.932148263480525,9.105590614303889,  10.0743311672404, 10.88453820307055,11.56523993876917,12.13635897512306,12.61239790477462,13.00438290346254]
EREF = [6.987522,8.580410,9.884711, 10.980483,11.909966, 12.699944, 13.369356,13.932540,14.400847]
EADC = [6.5173, 7.97249, 9.1355, 10.0959, 10.8996, 11.5752, 12.1421, 12.6146, 13.0034]

# plt.axis([0.039,0.205,6.0,15.5])
plt.xlabel(r'$\rho\, [\mathrm{fm}^{-3}]$', fontsize=18)
plt.ylabel(r'$E/A\,[\mathrm{MeV}]$', fontsize=18)
qmc   = plt.plot(density, ECIMC,'b-*',linewidth = 2.0, label = 'CIMC')
ccd   = plt.plot(density, ECCD,'r:.', linewidth = 2.0, label = 'CCD')
adc   = plt.plot(density, EADC,'r:.', linewidth = 2.0, label = 'ADC(3)')
# imsrg = plt.plot(density, EIMSRG-EREF,'r:.', linewidth = 2.0, label = 'IMSRG(2)')
# refenergy = plt.plot(density, EREF-EREF, 'm:v',linewidth = 2.0, label = 'Reference energy')


ax.set_xlim([0.03,0.21])  
ax.set_ylim([5.9,13.1])  

plt.legend(loc=4, borderaxespad=0.5)
# plt.savefig('imsrg_pnm.pdf', format='pdf',bbox_inches='tight')
plt.show()
plt.close()
exit()

ax.set_xlim([0.03,0.21])  
ax.set_ylim([-1.5,-0.4])  

plt.xlabel(r'$\rho\, [\mathrm{fm}^{-3}]$', fontsize=18)
plt.ylabel(r'$E_\mathrm{corr}/A\,[\mathrm{MeV}]$', fontsize=18)
qmc   = plt.plot(density, [a-b for a,b in zip(ECIMC, EREF)],'b-*',linewidth = 2.0, label = 'CIMC')
ccd   = plt.plot(density, [a-b for a,b in zip(ECCD,  EREF)],'r:.', linewidth = 2.0, label = 'CCD')
adc   = plt.plot(density, [a-b for a,b in zip(EADC,  EREF)],'r:.', linewidth = 2.0, label = 'ADC(3)')

plt.legend(loc=1, borderaxespad=0.5)
plt.savefig('imsrg_pnm_ecorr.pdf', format='pdf',bbox_inches='tight')
# plt.show()

