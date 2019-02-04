from data import *
import numpy as np 
import astropy.constants as con
import matplotlib.pyplot as plt 
import matplotlib
matplotlib.style.use('classic')
matplotlib.rc('xtick.major', size=15, width=3)
matplotlib.rc('xtick.minor', size=7.5, width=3)
matplotlib.rc('ytick.major', size=15, width=3)
matplotlib.rc('ytick.minor', size=7.5, width=3)
matplotlib.rc('lines',linewidth=4)
matplotlib.rc('axes', linewidth=4)

fig=plt.figure(figsize = (20,10))
ax = fig.add_axes([0.12,0.12,0.83,0.83])

lamb_norm=1e4

data=np.genfromtxt(datapath+"K13_SED.dat",names=["lognu","logall"],)
lamb=con.c.value/(10**data['lognu'])*1e10
ax.plot(lamb,data['logall'],c='royalblue',label=r'$\rm Krawczyk+$ $\rm 2013$')
fIR = data['logall'][np.abs(lamb-lamb_norm)<=100]

data=np.genfromtxt(datapath+"R06_SED.dat",names=["lognu","logall","sigall","blue"],)
lamb=con.c.value/(10**data['lognu'])*1e10
rescale = fIR/data['blue'][np.abs(lamb-lamb_norm)<=100]
ax.plot(lamb,rescale*data['blue'],c='crimson',label=r'$\rm Richards+$ $\rm 2006$ ($\rm blue$ $\rm quasars$)')
rescale = fIR/data['logall'][np.abs(lamb-lamb_norm)<=100]
ax.plot(lamb,rescale*data['logall'],c='chocolate',label=r'$\rm Richards+$ $\rm 2006$ ($\rm all$ $\rm quasars$)')



prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\lambda[\AA]$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\nu L_{\nu}[{\rm erg}/{\rm s}])}$',fontsize=40,labelpad=5)

ax.set_xlim(1e6,3e3)
ax.set_xscale('log')
ax.set_ylim(43.8,46)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
#plt.savefig("../figs/sed.pdf",fmt='pdf')
plt.show()

