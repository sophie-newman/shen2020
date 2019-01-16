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

data=np.genfromtxt(datapath+"K13_SED.dat",names=["lognu","logall"],)
lamb=con.c.value/(10**data['lognu'])*1e10
ax.plot(lamb,data['logall'],c='royalblue',label=r'$\rm Krawczyk+$ $\rm 2013$')

data=np.genfromtxt(datapath+"R06_SED.dat",names=["lognu","logall","sigall","blue"],)
lamb=con.c.value/(10**data['lognu'])*1e10
ax.plot(lamb,data['blue'],c='crimson',label=r'$\rm Richards+$ $\rm 2006$')
lamb=con.c.value/(10**data['lognu'])*1e10
ax.plot(lamb,data['logall'],c='chocolate')

data=np.genfromtxt(datapath+"XRAY_SED.dat",names=["lognu","logall"],)
lamb=con.c.value/(10**data['lognu'])*1e10
L_HX=45.1
ax.plot(lamb,(data['logall']-(-1.47406))+L_HX,c='seagreen',label=r'$\rm Xray$')

x=np.logspace(2.7,4.,1000)
ax.plot(x, 45.766+np.log10((912./x)**(-0.44+1)) ,'--', lw=2, c='cyan',label=r'$\rm Vanden$ $\rm Berk+$ $\rm 2001$')
x=np.logspace(np.log10(912),4.,1000)
ax.plot(x, 45.766+np.log10((912./x)**(-0.61+1)) ,'--', lw=2, c='gray',label=r'$\rm Lusso+$ $\rm 2015$')
x=np.logspace(np.log10(600),np.log10(912),1000)
ax.plot(x, 45.766+np.log10((912./x)**(-1.70+1)) ,'--', lw=2, c='gray')

ax.axvline(2500.,ymin=0.95)
ax.axvline(1216.,ymin=0.95)
ax.axvline(912.,ymin=0.95)
ax.axvline(con.c.value/(0.2*1000.*con.e.value/con.h.value)*1e10, ymin=0.95)
ax.axvline(con.c.value/(0.5*1000.*con.e.value/con.h.value)*1e10, ymin=0.95)

ax.axvspan(912,2500,color='seagreen',edgecolor='seagreen',alpha=0.1)
ax.text(0.45, 0.6, r'$\rm FUV$',color='gray',fontsize=25,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,rotation='vertical')
ax.axvspan(912,con.c.value/(0.2*1000.*con.e.value/con.h.value)*1e10,
	color='cyan',edgecolor='cyan',alpha=0.1)
ax.text(0.52, 0.6, r'$\rm EUV$',color='gray',fontsize=25,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,rotation='vertical')
ax.axvspan(con.c.value/(0.5*1000.*con.e.value/con.h.value)*1e10,con.c.value/(2*1000.*con.e.value/con.h.value)*1e10,
	color='navy',edgecolor='navy',alpha=0.1)
ax.text(0.78, 0.6, r'$\rm Soft$ $\rm Xray$',color='gray',fontsize=25,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,rotation='vertical')
ax.axvspan(con.c.value/(2*1000.*con.e.value/con.h.value)*1e10,con.c.value/(10*1000.*con.e.value/con.h.value)*1e10,
	color='darkorchid',edgecolor='darkorchid',alpha=0.1)
ax.text(0.88, 0.6, r'$\rm Hard$ $\rm Xray$',color='gray',fontsize=25,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,rotation='vertical')

prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\lambda[\AA]$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\nu L_{\nu}[{\rm erg}/{\rm s}])}$',fontsize=40,labelpad=5)

ax.set_xlim(1e6,1.)
ax.set_xscale('log')
ax.set_ylim(43.8,46)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/sed.pdf",fmt='pdf')
#plt.show()

