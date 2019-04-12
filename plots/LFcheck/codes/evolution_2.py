from data import *
import numpy as np 
from lf_shape import *
import scipy.interpolate as inter
from convolve import *
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import least_squares
# fit the luminosity function based on datasets at a given redshift
from lf_fitter_data import *
from ctypes import *
import ctypes
import sys

parameters_init = np.array([0.41698725, 2.17443860, -4.82506430, 13.03575300, 0.63150872, -11.76356000, -14.24983300, -0.62298947, 1.45993930, -0.79280099])

#load the shared object file
c_extenstion = CDLL(homepath+'codes/c_lib/convolve.so')
convolve_c = c_extenstion.convolve
convolve_c.restype = ctypes.POINTER(ctypes.c_double * N_bol_grid)

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

def plot_for_z(redshift,color):
        fit_res=np.genfromtxt("../../fitresult/fit_at_z_nofix.dat",names=True)
        id=fit_res["z"]==redshift
        parameters_free_local=np.array([ fit_res["gamma1"][id],fit_res["gamma2"][id],fit_res["phi_s"][id],fit_res["L_s"][id]])

        fit_res=np.genfromtxt("../../fitresult/fit_at_z.dat",names=True)
        id=fit_res["z"]==redshift
        parameters_fix_local=np.array([ fit_res["gamma1"][id],fit_res["gamma2"][id],fit_res["phi_s"][id],fit_res["L_s"][id]])

        fit_evolve=np.genfromtxt("../../Fit_parameters/codes/zevolution_fit.dat",names=['gamma1','gamma2','phis','Lbreak'])
        parameters_global_1 = pars_at_z(fit_evolve,redshift)

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

	x = L_bol_grid + L_solar 
	#y = LF(L_bol_grid,parameters_free_local)
	#ax.plot(x,y,'-',c=color,alpha=0.8,label=r'$\rm z=$'+str(redshift))
	y = LF(L_bol_grid,parameters_global_2)
	ax.plot(x,y,'-',c=color,label=r'$\rm z=$'+str(redshift))

	#x = L_bol_grid + L_solar 
	#y = LF_at_z_H07(L_bol_grid,parameters_init,redshift,"Fiducial")
	#ax.plot(x,y,'-',c=color,label=r'$\rm Hopkins+$ $\rm 2007$')

color_control = np.linspace(0,1,6)
i=0
for redshift in [0.2,0.4,0.8,1.2,1.8,2.4]:
	plot_for_z(redshift,(color_control[i],0,1-color_control[i]))
	i+=1

prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\log{(L_{\rm bol}[{\rm erg}\,{\rm s}^{-1}])}$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\phi[{\rm dex}^{-1}{\rm Mpc}^{-3}])}$',fontsize=40,labelpad=5)
ax.text(0.8, 0.9, r'$\rm Late$ $\rm phase$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=40)

ax.set_xlim(42.5,51.3)
ax.set_ylim(-11.2,-2.3)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/evolve_late.pdf",fmt='pdf')
#plt.show()

