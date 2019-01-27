from data import *
import numpy as np 
import astropy.constants as con
import matplotlib.pyplot as plt 
import matplotlib
from sherpa.astro import xspec
matplotlib.style.use('classic')
matplotlib.rc('xtick.major', size=15, width=3)
matplotlib.rc('xtick.minor', size=7.5, width=3)
matplotlib.rc('ytick.major', size=15, width=3)
matplotlib.rc('ytick.minor', size=7.5, width=3)
matplotlib.rc('lines',linewidth=4)
matplotlib.rc('axes', linewidth=4)

fig=plt.figure(figsize = (20,10))
ax = fig.add_axes([0.12,0.12,0.83,0.83])

model=xspec.XSpexrav(PhoIndex=1.8,foldE=500,rel_refl=1,redshift=0,cosIncl=0.5)
print model

data=np.genfromtxt(datapath+"XRAY_SED.dat",names=["lognu","logall"],)
lamb=con.c.value/(10**data['lognu'])*1e10
L_HX =  45
ax.plot(lamb,(data['logall']-(-1.47406))+L_HX,c='seagreen',label=r'$\rm Xray$')

prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\lambda[\AA]$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\nu L_{\nu}[{\rm erg}/{\rm s}])}$',fontsize=40,labelpad=5)

#ax.set_xlim(1e6,1.)
ax.set_xscale('log')
#ax.set_ylim(43.8,46)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
#plt.show()

