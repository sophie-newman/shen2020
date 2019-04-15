from data import *
import numpy as np 
from lf_shape import *
import scipy.interpolate as inter
from convolve import *
# fit the luminosity function based on datasets at a given redshift
from lf_fitter_data import *
from ctypes import *
import ctypes
import sys

redshift=float(sys.argv[1])

parameters_init = np.array([0.41698725, 2.17443860, -4.82506430, 13.03575300, 0.63150872, -11.76356000, -14.24983300, -0.62298947, 1.45993930, -0.79280099])

fit_res=np.genfromtxt("../../fitresult/fit_at_z.dat",names=True)
id=fit_res["z"]==redshift
parameters_fix_local=np.array([ fit_res["gamma1"][id],fit_res["gamma2"][id],fit_res["phi_s"][id],fit_res["L_s"][id]])

source = np.genfromtxt("../../Fit_parameters/codes/zevolution_fit_global.dat",names=True)
p=source['value'][ source['paraid']==0 ]
gamma1 = polynomial(redshift,p,2)
p=source['value'][ source['paraid']==1 ]
gamma2 = doublepower(redshift,p)
p=source['value'][ source['paraid']==2 ]
logphi = polynomial(redshift,p,1)
p=source['value'][ source['paraid']==3 ]
Lbreak = doublepower(redshift,p)
parameters_global_2 = np.array([gamma1,gamma2,logphi,Lbreak])

#load the shared object file
c_extenstion = CDLL(homepath+'codes/c_lib/convolve.so')
convolve_c = c_extenstion.convolve
convolve_c.restype = ctypes.POINTER(ctypes.c_double * N_bol_grid)

def lamb_vs_Lbol(Lbol):
	median = 0.469 * Lbol - 22.46
	sigma = 0.7
	return median, sigma

def Lbol_to_mass(Lbol):
	C = np.log10(1.26e38)
	median, sigma = lamb_vs_Lbol(Lbol)
	lamb = np.random.normal(median, sigma, size=nsample)
	return - lamb - C + Lbol

mock_volume = 1.
nsample = 10000
def generate_ensemble(Lbol, LogPhi):
	Masses = np.array([])
	Weights = np.array([])
	for i in range(len(Lbol)):
		Masses = np.append(Masses, Lbol_to_mass(Lbol[i]))
		Weights = np.append(Weights, np.ones(nsample)*10**LogPhi[i]*mock_volume)
	return Masses, Weights
	
def get_mass_function(Masses, Weights):
	bins = np.linspace(6., 11., 21)
	lenbin = bins[5]-bins[4]
	result,_ = np.histogram(Masses, weights=Weights, bins=bins)
	result = result/lenbin/mock_volume
	result = np.log10(result)
	return (bins[1:]+bins[:-1])/2., result

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

catalog = {"mass":0, "weight":0}

catalog["mass"], catalog["weight"] = generate_ensemble(L_bol_grid+L_solar, LF(L_bol_grid,parameters_global_2) )
x,y = get_mass_function(catalog["mass"],catalog["weight"])

ax.plot(x,y,c='royalblue')

prop = matplotlib.font_manager.FontProperties(size=22.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\log{(L_{\rm bol}[{\rm erg}\,{\rm s}^{-1}])}$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\phi[{\rm dex}^{-1}{\rm Mpc}^{-3}])}$',fontsize=40,labelpad=5)
ax.text(0.88, 0.92, r'${\rm z\sim'+str(redshift)+'}$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=40)

ax.set_xlim(5.5, 11.5)
#ax.set_ylim(-11.2,-2.3)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
#plt.savefig("../figs/bol_"+str(redshift)+".pdf",fmt='pdf')
plt.show()

