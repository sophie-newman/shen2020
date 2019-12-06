from data import *
import numpy as np 
from lf_shape import *
from new_load_kk18_lf_shape import *
import scipy.interpolate as inter
from scipy.integrate import quad
from scipy.integrate import romberg
from convolve import *
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import least_squares
import astropy.constants as con
# fit the luminosity function based on datasets at a given redshift
from ctypes import *
import ctypes
import sys

# parameters of the H07 model
parameters_init = np.array([0.41698725, 2.17443860, -4.82506430, 13.03575300, 0.63150872, -11.76356000, -14.24983300, -0.62298947, 1.45993930, -0.79280099])

# best-fit in our global fit
fit_evolve=np.genfromtxt("../../Fit_parameters/codes/zevolution_fit_global.dat",names=True)
paraid, pglobal, pglobal_err = fit_evolve['paraid'], fit_evolve['value'], (fit_evolve['uperr']+fit_evolve['loerr'])/2.
zlist=np.linspace(0.1,7,100)

#load the shared object file
c_extenstion = CDLL(homepath+'codes/c_lib/convolve.so')
convolve_c = c_extenstion.convolve
convolve_c.restype = ctypes.POINTER(ctypes.c_double * N_bol_grid)

def get_model_lf_global(parameters,redshift):
	zref = 2.
	p=parameters[paraid==0]
	gamma1 = polynomial(redshift,p,2)
	p=parameters[paraid==1]
	gamma2 = doublepower(redshift,(p[0],zref, p[1], p[2]))
	p=parameters[paraid==2]
	logphi = polynomial(redshift,p,1)
	p=parameters[paraid==3]	
	Lbreak = doublepower(redshift,(p[0],zref, p[1], p[2]))
	parameters_at_z = np.array([gamma1,gamma2,logphi,Lbreak])
	return get_model_lf(parameters_at_z)

def get_model_lf(parameters):
        return L_bol_grid, LF(L_bol_grid,parameters)

def cumulative_lum(L_bol,Phi_bol,L_limit_low,L_limit_up):
	logphi=inter.interp1d(L_bol,Phi_bol)
	def func(x):
		return np.power(10.,logphi(x)) * 10**(x+L_solar)
	#return quad(func,L_limit_low,L_limit_up)[0]
	return romberg(func,L_limit_low,L_limit_up,divmax=20)

eps = 0.1 #our fiducial choice of the radiative efficiency
def rhoBH(zstart=7):
	zlist = np.linspace(zstart,0,101)
	rholist = 0.0*zlist
	rholist[0] = 1e-30
	for i,z in enumerate(zlist):
		if i!=0:
			K = (1-eps)/(cosmo.H(z).value*(1000./1e6/con.pc.value)*(1+z))/eps/con.c.value**2 #prefactor
			L_bol, Phi_bol = get_model_lf_global(pglobal, z)
			rholist[i] = rholist[i-1] + K * cumulative_lum(L_bol, Phi_bol, 43-L_solar, 48-L_solar)*1e-7 * np.abs(zlist[i]-zlist[i-1])  
	rholist = rholist/(con.M_sun.value)
	return zlist, rholist

import matplotlib.pyplot as plt 
import matplotlib

matplotlib.style.use('classic')
matplotlib.rc('xtick.major', size=15, width=3)
matplotlib.rc('xtick.minor', size=7.5, width=3)
matplotlib.rc('ytick.major', size=15, width=3)
matplotlib.rc('ytick.minor', size=7.5, width=3)
matplotlib.rc('lines',linewidth=4)
matplotlib.rc('axes', linewidth=4)

fig=plt.figure(figsize = (15,10))
ax = fig.add_axes([0.13,0.12,0.79,0.83])

x,y = rhoBH()
ax.plot(x,np.log10(y),c='royalblue',label=r'$\rm This$ $\rm work$'+'\n'+r'$\rm integrated$ $\rm from$ $z_{\rm i}=7$')
ax.fill_between(x,y1=np.log10(y/(0.9/0.1)*(0.8/0.2)),y2=np.log10(y/(0.9/0.1)*(0.95/0.05)),color='royalblue',alpha=0.2)

x,y = rhoBH(10)
ax.plot(x,np.log10(y),c='crimson',label=r'$\rm integrated$ $\rm from$ $z_{\rm i}=10$')
ax.fill_between(x,y1=np.log10(y/(0.9/0.1)*(0.8/0.2)),y2=np.log10(y/(0.9/0.1)*(0.95/0.05)),color='crimson',alpha=0.2)

x,y = rhoBH(4)
ax.plot(x,np.log10(y),c='seagreen',label=r'$\rm integrated$ $\rm from$ $z_{\rm i}=4$')
ax.fill_between(x,y1=np.log10(y/(0.9/0.1)*(0.8/0.2)),y2=np.log10(y/(0.9/0.1)*(0.95/0.05)),color='seagreen',alpha=0.2)

ax.errorbar(0.,np.log10(4.81*1e5),yerr=([np.log10(4.81)-np.log10(4.81-0.99)],[np.log10(4.81+1.24)-np.log10(4.81)]),marker='o',c='gray',mec='gray',ms=15,capsize=9,capthick=4,label=r'$\rm Hopkins+$ $\rm 2007$',linestyle='none')

ax.errorbar(0-0.25,np.log10(4.2*1e5),yerr=([np.log10(4.2)-np.log10(4.2-1.1)],[np.log10(4.2+1.1)-np.log10(4.2)]),marker='o',c='navy',mec='navy',ms=15,capsize=9,capthick=4,label=r'$\rm Shankar+$ $\rm 2001$',linestyle='none')

ax.errorbar(0-0.50,np.log10(4.6*1e5),yerr=([np.log10(4.6)-np.log10(4.6-1.4)],[np.log10(4.6+1.9)-np.log10(4.6)]),marker='o',c='k',mec='k',ms=15,capsize=9,capthick=4,label=r'$\rm Marconi+$ $\rm 2004$',linestyle='none')

ax.errorbar(0-0.75,np.log10(4.41*1e5),yerr=([np.log10(4.41)-np.log10(4.41-1.67)],[np.log10(4.41+1.67)-np.log10(4.41)]),marker='o',c='darkorchid',mec='darkorchid',ms=15,capsize=9,capthick=4,label=r'$\rm Graham+$ $\rm 2007$',linestyle='none')

ax.errorbar(0-1.00,np.log10(3.8*1e5),yerr=([np.log10(3.8)-np.log10(3.8-0.6)],[np.log10(3.8+0.7)-np.log10(3.8)]),marker='o',c='chocolate',mec='chocolate',ms=15,capsize=9,capthick=4,label=r'$\rm Yu+$ $\rm 2008$',linestyle='none')

ax.plot(0-1.25, np.log10(4.5*1e5*0.1/0.075), marker='o',c='pink',mec='pink',ms=15,label=r'$\rm Shankar+$ $\rm 2009$',linestyle='none')

#######################################################################
ax.axvspan(-2,0,color='cyan',alpha=0.3)

prop = matplotlib.font_manager.FontProperties(size=21.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=1,ncol=1,frameon=False)
ax.set_xlabel(r'$\rm z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\rho_{\rm BH}\,[{\rm M}_{\odot}\,{\rm Mpc}^{-3}])}$',fontsize=40,labelpad=5)

#ax.text(0.25, 0.87, r'$\rm <-18$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=30,color='gray')

ax.set_xlim(-1.4,7.2)
ax.set_ylim(2.2,6.1)

from matplotlib.ticker import FixedLocator, FixedFormatter
#ax.minorticks_on()
x_formatter = FixedFormatter(["0", "1", "2", "3", "4", "5", "6", "7"])
x_locator   = FixedLocator([0, 1, 2, 3, 4, 5, 6, 7])
xm_locator  = FixedLocator(np.linspace(0,7,int(7/0.2+1)))
ax.xaxis.set_major_formatter(x_formatter)
ax.xaxis.set_major_locator(x_locator)
ax.xaxis.set_minor_locator(xm_locator)

ym_locator = FixedLocator(np.linspace(2.0,7.0,int((7.0-2.0)/0.1+1)))
ax.yaxis.set_minor_locator(ym_locator)

ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
plt.savefig("../figs/bh_mass_density.pdf",fmt='pdf')
#plt.show()

