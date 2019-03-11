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

data=np.genfromtxt(datapath+"MySED.dat",names=["lamb","logall"])
ax.plot(data['lamb'],data['logall'],lw=5,c='royalblue',label=r'$\rm This$ $\rm work$')
myf2500 = data['logall'][(data['lamb']-2500)<0][0]

data=np.genfromtxt(datapath+"K13_SED.dat",names=["lognu","logall"],)
lamb=con.c.value/(10**data['lognu'])*1e10
scale = 0#myf2500 - data['logall'][(lamb-2500)<0][0]
ax.plot(lamb[lamb>912],data['logall'][lamb>912]+scale,'--',dashes=(15,9),c='cyan',label=r'$\rm Krawczyk+$ $\rm 2013$')

data=np.genfromtxt(datapath+"R06_SED.dat",names=["lognu","logall","sigall","blue"],)
lamb=con.c.value/(10**data['lognu'])*1e10
scale = 0#myf2500 - data['blue'][(lamb-2500)<0][0]
ax.plot(lamb[lamb>500],data['blue'][lamb>500]+scale,'--',dashes=(15,9),c='crimson',label=r'$\rm Richards+$ $\rm 2006$ ($\rm blue$ $\rm quasars$)'
	+'\n'+r'$\rm adopted$ $\rm in$ $\rm Hopkins+$ $\rm 2007$')
lamb=con.c.value/(10**data['lognu'])*1e10
scale = 0#myf2500 - data['logall'][(lamb-2500)<0][0]
ax.plot(lamb[lamb>500],data['logall'][lamb>500]+scale,'--',dashes=(15,9),c='chocolate',label=r'$\rm Richards+$ $\rm 2006$ ($\rm all$ $\rm quasars$)')
'''
data=np.genfromtxt(datapath+"OpticalSED.dat",names=["lognu","logall"])
lamb=con.c.value/(10**data['lognu'])*1e10
scale = myf2500 - data['logall'][(lamb-2500)<0][0]
ax.plot(lamb,data['logall']+scale,'--',c='darkorchid',label=r'$\rm Vanden$ $\rm Berk+$ $\rm 2001$')
'''

f2500 = np.power(10., myf2500)/(con.c.value/2500./1e-10)
alphaox = -0.137*np.log10(f2500)+2.638
f2kev=10**(alphaox/0.384)*f2500
data=np.genfromtxt(datapath+"XRAY_SED.dat",names=["lognu","logall"],)
lamb=con.c.value/(10**data['lognu'])*1e10
ratio = data['logall'][np.abs(10**data['lognu']*con.h.value/con.e.value/1000-2)<=0.05]-(-1.47406)
L_HX =  np.log10( f2kev*(2.*1000.*con.e.value/con.h.value) / 10**ratio )
ax.plot(lamb[lamb<50],(data['logall'][lamb<50]-(-1.47406))+L_HX,'--',dashes=(15,9),c='yellow',label=r'$\rm Xray$ $\rm Hopkins+$ $\rm 2007$')

x=np.logspace(2.7,4.,1000)
ax.plot(x, 45.766+np.log10((912./x)**(-0.44+1)) ,'-', lw=2, c='seagreen',label=r'$\rm Vanden$ $\rm Berk+$ $\rm 2001$')
x=np.logspace(np.log10(912),4.,1000)
ax.plot(x, 45.766+np.log10((912./x)**(-0.61+1)) ,'-', lw=2, c='k',label=r'$\rm Lusso+$ $\rm 2015$')
x=np.logspace(np.log10(600),np.log10(912),1000)
ax.plot(x, 45.766+np.log10((912./x)**(-1.70+1)) ,'-', lw=2, c='k')

#ax.axvline(2500.,ymin=0.95)
#ax.axvline(1216.,ymin=0.95)
#ax.axvline(912.,ymin=0.95)
ax.axvline(500.,ymax=0.86,ymin=0.76,color='gray')
ax.axvline(50., ymax=0.15,ymin=0.05,color='gray')
#ax.axvline(con.c.value/(0.2*1000.*con.e.value/con.h.value)*1e10, ymin=0.95)
#ax.axvline(con.c.value/(0.5*1000.*con.e.value/con.h.value)*1e10, ymin=0.95)

ax.axvspan(10000,1e6,color='tan',edgecolor='tan',alpha=0.1)
ax.text(0.25, 0.6, r'$\rm Infrared$',color='gray',fontsize=25,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,rotation='vertical')
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

