#from data import *
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
'''
f100kev=1
data=np.genfromtxt("Xspec.dat",names=["E","f"],)
lamb=con.c.value/(data['E']*1000*con.e.value/con.h.value)*1e10
sed = np.log10(data['f']*data['E'])
id = (lamb - con.c.value/(100.*1000.*con.e.value/con.h.value)*1e10)<0
scale =  np.log10( f100kev*(100.*1000.*con.e.value/con.h.value)) - sed[id][0]
ax.plot(lamb, sed + scale,'-',c='cyan',label=r'$\rm repro$',alpha=0.4)

f100kev=1
data=np.genfromtxt("Xspec_ueda14.dat",names=["E","Ef"],)
lamb=con.c.value/(data['E']*1000*con.e.value/con.h.value)*1e10
sed = np.log10(data['Ef'])
id = (lamb - con.c.value/(100.*1000.*con.e.value/con.h.value)*1e10)<0
scale =  np.log10( f100kev*(100.*1000.*con.e.value/con.h.value)) - sed[id][0]
ax.plot(lamb, sed + scale,'-',c='crimson',label=r'$\rm repro$',alpha=0.4)
'''
f10kev=1
data=np.genfromtxt("Xspec.dat",names=["E","f"],)
lamb=con.c.value/(data['E']*1000*con.e.value/con.h.value)*1e10
sed = np.log10(data['f']*data['E'])
id = (lamb - con.c.value/(10.*1000.*con.e.value/con.h.value)*1e10)<0
scale =  np.log10( f10kev*(10.*1000.*con.e.value/con.h.value)) - sed[id][0]
ax.plot(lamb, sed + scale,'-',c='cyan',label=r'$\rm repro$',alpha=0.4)

f10kev=1
data=np.genfromtxt("Xspec_aird15.dat",names=["E","Ef"],)
lamb=con.c.value/(data['E']*1000*con.e.value/con.h.value)*1e10
sed = np.log10(data['Ef'])
id = (lamb - con.c.value/(10.*1000.*con.e.value/con.h.value)*1e10)<0
scale =  np.log10( f10kev*(10.*1000.*con.e.value/con.h.value)) - sed[id][0]
ax.plot(lamb, sed + scale,'-',c='crimson',label=r'$\rm a15$',alpha=0.4)

'''
f100kev=1
data=np.genfromtxt("Xspec_1.dat",names=["E","f"],)
lamb=con.c.value/(data['E']*1000*con.e.value/con.h.value)*1e10
sed = np.log10(data['f']*data['E'])
id = (lamb - con.c.value/(100.*1000.*con.e.value/con.h.value)*1e10)<0
scale =  np.log10( f100kev*(100.*1000.*con.e.value/con.h.value)) - sed[id][0]
ax.plot(lamb, sed + scale,'-',c='cyan',label=r'$\rm repro$',alpha=0.4)


f100kev=1
data=np.genfromtxt(datapath+"XRAY_SED.dat",names=["lognu","logall"],)
lamb=con.c.value/(10**data['lognu'])*1e10
sed = data['logall']
id = (lamb - con.c.value/(100.*1000.*con.e.value/con.h.value)*1e10)<0
scale =  np.log10( f100kev*(100.*1000.*con.e.value/con.h.value) ) - sed[id][0]
ax.plot(lamb, sed + scale,'--',dashes=(15,9),c='k',label=r'$\rm Xray$ $\rm Hopkins+$ $\rm 2007$')
'''

ax.axvline(con.c.value/(0.5*1000.*con.e.value/con.h.value)*1e10, ymin=0.95)
ax.axvline(con.c.value/(2*1000.*con.e.value/con.h.value)*1e10, ymin=0.95)
ax.axvline(con.c.value/(10*1000.*con.e.value/con.h.value)*1e10, ymin=0.95)


prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\lambda[\AA]$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\nu L_{\nu}[{\rm erg}/{\rm s}])}$',fontsize=40,labelpad=5)

ax.set_xlim(1e3,0.001)
ax.set_xscale('log')
ax.set_ylim(17,20.5)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
#plt.savefig("../figs/sed.pdf",fmt='pdf')
plt.show()

