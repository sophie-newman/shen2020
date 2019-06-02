from data import *
import numpy as np 
from lf_shape import *
import scipy.interpolate as inter
from scipy.integrate import quad
from convolve import *
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import least_squares
import astropy.constants as con
# fit the luminosity function based on datasets at a given redshift
from ctypes import *
import ctypes
import sys

parameters_init = np.array([0.41698725, 2.17443860, -4.82506430, 13.03575300, 0.63150872, -11.76356000, -14.24983300, -0.62298947, 1.45993930, -0.79280099])

fit_evolve=np.genfromtxt("../../Fit_parameters/codes/zevolution_fit_global.dat",names=True)
paraid, pglobal, pglobal_err = fit_evolve['paraid'], fit_evolve['value'], (fit_evolve['uperr']+fit_evolve['loerr'])/2.

#load the shared object file
c_extenstion = CDLL(homepath+'codes/c_lib/convolve_new.so')
convolve_c = c_extenstion.convolve
convolve_c.restype = ctypes.POINTER(ctypes.c_double * N_bol_grid)

def get_model_lf(parameters,nu,redshift):
        L_band = bolometric_correction(L_bol_grid,nu)
        nu_c = c_double(nu)
	redshift_c = c_double(redshift)
        input_c= np.power(10.,LF(L_bol_grid,parameters)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        res = convolve_c(input_c,nu_c,redshift_c)
        res = [i for i in res.contents]
        PHI_band = np.array(res,dtype=np.float64)
	return L_band, np.log10(PHI_band)

def get_model_lf_global(nu,redshift):
	parameters=pglobal
	p=parameters[paraid==0]
	gamma1 = polynomial(redshift,p,2)
	p=parameters[paraid==1]
	gamma2 = doublepower(redshift,p)
	p=parameters[paraid==2]
	logphi = polynomial(redshift,p,1)
	p=parameters[paraid==3]	
	Lbreak = doublepower(redshift,p) + 0.5
	parameters_at_z = np.array([gamma1,gamma2,logphi,Lbreak])
	return get_model_lf(parameters_at_z,nu,redshift)

def cumulative_emissivity(L_nu,Phi_nu,L_limit_low,L_limit_up,nu):
	logphi=inter.interp1d(L_nu, Phi_nu)
	def emis(x):
		return np.power(10.,logphi(x))*np.power(10.,x)/nu

	return quad(emis, L_limit_low, L_limit_up)[0]*np.power(10.,L_solar)

def to_be_integrate(z, nuobs):
	nuem = nuobs*(1+z)
        L_nu, PHI_nu = get_model_lf_global(nuem, z)
        emissivity = cumulative_emissivity(L_nu, PHI_nu, L_nu[0], L_nu[-1], nuem) 
	return emissivity/4./np.pi/(cosmo.luminosity_distance(z).value*1e6*con.pc.value*1e2)**2  * cosmo.differential_comoving_volume(z).value

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

E_list = np.logspace(-0.5,3,100)
nu_list = E_list*1000.*con.e.value/con.h.value

zbins = np.linspace(0,7,50)
zcenters = (zbins[1:]+zbins[:-1])/2.
deltaz= zbins[5]-zbins[4]

Intensity = 0.0*E_list
unit_convertion = 1e-7/(1000.*con.e.value)

for i in range(len(E_list)):
	Intensity[i]=0
	for j in range(len(zcenters)):
		Intensity[i] += to_be_integrate(zcenters[j],nu_list[i])*deltaz

	#Intensity[i]=quad( to_be_integrate, 0, 7, args=(nu_list[i]) )[0]/(con.pc.value*1e6*1e2)**2
	print i

#data=np.genfromtxt("ajello2008.dat",names=True)
#scale = np.max(data['CXB'])/np.max(nu_list * Intensity * unit_convertion)
ax.plot(E_list, nu_list * Intensity * unit_convertion, c='royalblue', label=r'$\rm This$ $\rm work$')
print nu_list * Intensity * unit_convertion

data=np.genfromtxt("ajello2008.dat",names=True)
ax.errorbar(data['E'],data['CXB'],xerr=data['dE'],yerr=data['dCXB'],c='crimson',mec='crimson',capsize=0,capthick=0,linestyle='none',marker='.',lw=2,ms=1,label=r'$\rm Ajello+$ $\rm 2008$')

data=np.genfromtxt("churazov2007.dat",names=True)
ax.errorbar(data['E'],data['CXB'],xerr=(data['E']-data['Elo'],data['Eup']-data['E']),yerr=(data['CXB']-data['lo'],data['up']-data['CXB']),c='seagreen',mec='seagreen',capsize=0,capthick=0,linestyle='none',marker='.',lw=2,ms=1,label=r'$\rm Churazov+$ $\rm 2007$')

prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\rm E$ [$\rm keV$]',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\nu I_{\rm XRB,\nu}$ [$\rm keV\,s^{-1}\,cm^{-2}\,sr^{-1}$]',fontsize=40,labelpad=5)

#ax.text(0.25, 0.64, r'$\rm <-21$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=30,color='gray')

ax.set_xlim(np.min(E_list),np.max(E_list))
#ax.set_ylim(-2.5,2)
ax.set_yscale('log')
ax.set_xscale('log')
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/CXB_new.pdf",fmt='pdf')
#plt.show()
