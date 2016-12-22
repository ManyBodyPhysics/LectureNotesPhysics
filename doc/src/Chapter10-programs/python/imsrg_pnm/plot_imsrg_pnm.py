#!/usr/bin/env python
import matplotlib.pyplot as plt

from sympy import *
from pylab import *
from matplotlib import rc

rc('font',**{'size':14, 'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)

def myLabels(x, pos):
  return '$%s$'%x


fig = figure(figsize=(9,6))
ax  = fig.gca()

density = [0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2]
ECIMC  = [6.4735287261, 7.9714751994, 9.1357756052, 10.0941195187,10.899839927, 
			11.5775331731, 12.1433034824, 12.6242634614, 13.0118757745 ]
ECCD   = [6.468200479265218, 7.932148263480525,9.105590614303889,  10.0743311672404, 
			10.88453820307055,11.56523993876917,12.13635897512306,12.61239790477462,13.00438290346254]
EREF   = [6.987522,8.580410,9.884711, 10.980483,11.909966, 12.699944, 13.369356,13.932540,14.400847]
EADC   = [6.5173, 7.97249, 9.1355, 10.0959, 10.8996, 11.5752, 12.1421, 12.6146, 13.0034]
EIMSRG = [6.5038401494, 7.9557647416, 9.1161816114, 10.074745431, 10.877348180, 11.552285089, 12.118922813,
			12.591409488, 12.980552843]

# plt.axis([0.039,0.205,6.0,15.5])
plt.xlabel(r'$\rho\, [\mathrm{fm}^{-3}]$', fontsize=18)
plt.ylabel(r'$E/A\,[\mathrm{MeV}]$', fontsize=18)
qmc   = plt.plot(density, ECIMC,	markersize=8, color='orange', marker='v', linestyle='-',  linewidth = 2.0, label = 'CIMC')
adc   = plt.plot(density, EADC,		markersize=7, color='green', marker='D',  linestyle='--', linewidth = 2.0, label = 'ADC(3)')
ccd   = plt.plot(density, ECCD,		markersize=8, color='red',  marker='s',   dashes=[8,6],   linewidth = 2.0, label = 'CCD')
imsrg = plt.plot(density, EIMSRG,	markersize=8, color='blue', marker='o',   linestyle=':',  linewidth = 2.0, label = 'IMSRG(2)')
refenergy = plt.plot(density, EREF, color='black' ,linewidth = 2.0, label = 'Reference energy')

ax.xaxis.set_major_formatter(FuncFormatter(myLabels))
ax.yaxis.set_major_formatter(FuncFormatter(myLabels))
ax.tick_params(axis='both',width=2,length=10,labelsize=18)

ax.tick_params(axis='both',which='major',width=1.5,length=8)
ax.tick_params(axis='both',which='minor',width=1.5,length=5)
ax.minorticks_on()
for s in ['left', 'right', 'top', 'bottom']:
	ax.spines[s].set_linewidth(2)
ax.set_xlim([0.03,0.21])  
ax.set_ylim([5.9,14.6])  

plt.legend(frameon=false, loc=4, borderaxespad=0.5)
plt.savefig('imsrg_pnm.pdf', format='pdf',bbox_inches='tight')
plt.show()
# plt.close()

fig = figure(figsize=(9,6))
ax  = fig.gca()
ax.xaxis.set_major_formatter(FuncFormatter(myLabels))
ax.yaxis.set_major_formatter(FuncFormatter(myLabels))
ax.tick_params(axis='both',width=2,length=10,labelsize=18)

ax.tick_params(axis='both',which='major',width=1.5,length=8)
ax.tick_params(axis='both',which='minor',width=1.5,length=5)
ax.minorticks_on()
for s in ['left', 'right', 'top', 'bottom']:
	ax.spines[s].set_linewidth(2)

ax.set_xlim([0.03,0.21])  
ax.set_ylim([-1.5,-0.4])  

plt.xlabel(r'$\rho\, [\mathrm{fm}^{-3}]$', fontsize=18)
plt.ylabel(r'$E_\mathrm{corr}/A\,[\mathrm{MeV}]$', fontsize=18)
qmc   = plt.plot(density, [a-b for a,b in zip(ECIMC, EREF)],   markersize=8, color='orange', marker='v',linestyle='-',  linewidth = 2.0, label = 'CIMC')
adc   = plt.plot(density, [a-b for a,b in zip(EADC,  EREF)],   markersize=7, color='green', marker='D', linestyle='--', linewidth = 2.0, label = 'ADC(3)')
ccd   = plt.plot(density, [a-b for a,b in zip(ECCD,  EREF)],   markersize=8, color='red',  marker='s',  dashes=[8,6],   linewidth = 2.0, label = 'CCD')
imsrg = plt.plot(density, [a-b for a,b in zip(EIMSRG,  EREF)], markersize=8, color='blue', marker='o',  linestyle=':',  linewidth = 2.0, label = 'IMSRG(2)')

plt.legend(frameon=false, loc=1, borderaxespad=0.5)
plt.savefig('imsrg_pnm_ecorr.pdf', format='pdf',bbox_inches='tight')
plt.show()

