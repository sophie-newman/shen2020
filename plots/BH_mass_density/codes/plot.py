from data import *
import numpy as np 
from lf_shape import *
from new_load_kk18_lf_shape import *
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

data=np.genfromtxt("../../fitresult/fit_at_z.dat",names=True)
pgamma1_fix, pgamma1_err_fix  = data["gamma1"], data["err1"]
pgamma2_fix, pgamma2_err_fix  = data["gamma2"], data["err2"]
plogphis_fix,plogphis_err_fix = data["phi_s"],  data["err3"]
pLbreak_fix, pLbreak_err_fix  = data["L_s"],    data["err3"]
pz_fix=data['z']
zpoints_fix=np.array(pz_fix)

data=np.genfromtxt("../../fitresult/fit_at_z_nofix.dat",names=True)
pgamma1_free, pgamma1_err_free  = data["gamma1"], data["err1"]
pgamma2_free, pgamma2_err_free  = data["gamma2"], data["err2"]
plogphis_free,plogphis_err_free = data["phi_s"],  data["err3"]
pLbreak_free, pLbreak_err_free  = data["L_s"],    data["err3"]
pz_free=data['z']
zpoints_free=np.array(pz_free)

fit_evolve=np.genfromtxt("../../Fit_parameters/codes/zevolution_fit_global.dat",names=True)
paraid, pglobal, pglobal_err = fit_evolve['paraid'], fit_evolve['value'], (fit_evolve['uperr']+fit_evolve['loerr'])/2.
zlist=np.linspace(0.1,7,100)

#load the shared object file
c_extenstion = CDLL(homepath+'codes/c_lib/convolve.so')
convolve_c = c_extenstion.convolve
convolve_c.restype = ctypes.POINTER(ctypes.c_double * N_bol_grid)

def get_model_lf_global(parameters,nu,redshift,magnitude=False):
	p=parameters[paraid==0]
	gamma1 = polynomial(redshift,p,2)
	p=parameters[paraid==1]
	gamma2 = doublepower(redshift,p)
	p=parameters[paraid==2]
	logphi = polynomial(redshift,p,1)
	p=parameters[paraid==3]	
	Lbreak = doublepower(redshift,p)
	parameters_at_z = np.array([gamma1,gamma2,logphi,Lbreak])
	return get_model_lf(parameters_at_z,nu,magnitude=magnitude)

def get_model_lf(parameters,nu,magnitude=False):
        L_band = bolometric_correction(L_bol_grid,nu)
        nu_c = c_double(nu)
        input_c= np.power(10.,LF(L_bol_grid,parameters)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        res = convolve_c(input_c,nu_c)
        res = [i for i in res.contents]
        PHI_band = np.array(res,dtype=np.float64)
        if magnitude==False:
                return L_band, np.log10(PHI_band)
        else:
                M_1450 = (M_sun_Bband_AB -2.5*L_band) + 0.706
                PHI_1450 = np.log10(PHI_band) - np.log10(2.5)
                return M_1450, PHI_1450

def cumulative_emissivity(L_band,Phi_band,L_limit_low,L_limit_up):
	Mband, Mlow, Mup = L_band, L_limit_low, L_limit_up
	logphi=inter.interp1d(Mband,Phi_band)
	def emis(x):
		fnu = np.power(10.,-0.4*x)*3631*1e-23*4*np.pi*(10*con.pc.value*100)**2
		fnu912 = fnu * (912./1450.)**(0.61)
		return np.power(10.,logphi(x))*fnu
	return quad(emis,Mlow,Mup)[0]

def Gamma_err(parameters,errs,L_limit_low,L_limit_up,redshift,global_fit=False):
	partials = 0.0 * parameters
	delta = 1e-6
	if global_fit==False:
		def fobjective(parameters):
			M_1450, PHI_1450 = get_model_lf(parameters, -1, magnitude=True)
        		return Gamma(cumulative_emissivity(M_1450, PHI_1450, L_limit_low, L_limit_up),redshift)
		for i in range(len(parameters)):
			parameters_add = parameters.copy()
			parameters_add[i] += delta
			partials[i] = ( fobjective(parameters_add) - fobjective(parameters))/delta
			partials[i] = np.abs(partials[i])
	else:
		def fobjective(parameters):
                        M_1450, PHI_1450 = get_model_lf_global(parameters, -1, redshift, magnitude=True)
                        return Gamma(cumulative_emissivity(M_1450, PHI_1450, L_limit_low, L_limit_up),redshift)
                for i in range(len(parameters)):
                        parameters_add = parameters.copy()
                        parameters_add[i] += delta
                        partials[i] = ( fobjective(parameters_add) - fobjective(parameters))/delta
                        partials[i] = np.abs(partials[i])

	return np.sqrt(np.sum((partials * errs)**2))

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

lowlimit=-35

result=np.zeros((len(zlist),2))
for i in range(len(zlist)):
	L_band = bolometric_correction(L_bol_grid,-1)
        nu_c = c_double(-1)
        input_c= np.power(10.,LF_at_z_H07(L_bol_grid,parameters_init,zlist[i],"Fiducial")).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        res = convolve_c(input_c,nu_c)
        res = [j for j in res.contents]
        PHI_band = np.array(res,dtype=np.float64)
	M_1450 = (M_sun_Bband_AB -2.5*L_band) + 0.706
	PHI_1450 = np.log10(PHI_band) - np.log10(2.5)

	result[i,0]= Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, -18),zlist[i])
	result[i,1]= Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, -21),zlist[i])
ax.plot(zlist,np.log10(result[:,0]),'--',dashes=(25,15),c='crimson',label=r'$\rm Hopkins+$ $\rm 2007$')

result=np.zeros((len(zlist),2))
uncertainty=np.zeros((len(zlist),2))
for i in range(len(zlist)):
	M_1450, PHI_1450 = get_model_lf_global(pglobal, -1, zlist[i], magnitude=True)
	result[i,0]= Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, -18),zlist[i])
	result[i,1]= Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, -21),zlist[i])
	uncertainty[i,0]= Gamma_err(pglobal, pglobal_err, lowlimit, -18, zlist[i], global_fit=True)
uperr = np.log10(result + uncertainty)-np.log10(result)
loerr = np.log10(result)-np.log10(result - uncertainty)
loerr[np.invert(np.isfinite(loerr))] = 100
ax.plot(zlist,np.log10(result[:,0]),'-',c='darkorchid',label=r'$\rm Global$ $\rm fit$')
ax.fill_between(zlist, y1=np.log10(result[:,0])+uperr[:,0] ,y2=np.log10(result[:,0])-loerr[:,0], color='darkorchid', alpha=0.4)

###### fit at a given redshift
result=np.zeros((len(zpoints_free),2))
uncertainty=np.zeros((len(zpoints_free),2))
for i in range(len(zpoints_free)):
        id = pz_free==zpoints_free[i]
        M_1450, PHI_1450 = get_model_lf([pgamma1_free[id],pgamma2_free[id],plogphis_free[id],pLbreak_free[id]], -1, magnitude=True)
        result[i,0]= Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, -18),zpoints_free[i])
        uncertainty[i,0]= Gamma_err(np.array([pgamma1_free[id],pgamma2_free[id],plogphis_free[id],pLbreak_free[id]]), [pgamma1_err_free[id],pgamma2_err_free[id],plogphis_err_free[id],pLbreak_err_free[id]], lowlimit, -18, zpoints_free[i])
        result[i,1]= Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, -21),zpoints_free[i])
uperr = np.log10(result + uncertainty)-np.log10(result)
loerr = np.log10(result)-np.log10(result - uncertainty)
loerr[np.invert(np.isfinite(loerr))] = 100
ax.errorbar(zpoints_free,np.log10(result[:,0]),yerr=(loerr[:,0],uperr[:,0]),linestyle='none',marker='o',c='gray',mec='gray',ms=15,capsize=9,capthick=4,alpha=0.5)


result=np.zeros((len(zpoints_fix),2))
uncertainty=np.zeros((len(zpoints_fix),2))
for i in range(len(zpoints_fix)):
        id= pz_fix==zpoints_fix[i]
        M_1450, PHI_1450 = get_model_lf([pgamma1_fix[id],pgamma2_fix[id],plogphis_fix[id],pLbreak_fix[id]], -1, magnitude=True)
        result[i,0]= Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, -18),zpoints_fix[i])
	uncertainty[i,0]= Gamma_err(np.array([pgamma1_fix[id],pgamma2_fix[id],plogphis_fix[id],pLbreak_fix[id]]), [pgamma1_err_fix[id],pgamma2_err_fix[id],plogphis_err_fix[id],pLbreak_err_fix[id]], lowlimit, -18, zpoints_fix[i])
	result[i,1]= Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, -21),zpoints_fix[i])
uperr = np.log10(result + uncertainty)-np.log10(result)
loerr = np.log10(result)-np.log10(result - uncertainty)
loerr[np.invert(np.isfinite(loerr))] = 100
ax.errorbar(zpoints_fix,np.log10(result[:,0]),yerr=(loerr[:,0],uperr[:,0]),linestyle='none',marker='o',c='royalblue',mec='royalblue',ms=15,capsize=9,capthick=4)

#######################################################################
prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\rm z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\rho_{\rm BH}$ [${\rm M}_{\odot}\,{\rm Mpc}^{-3}$]',fontsize=40,labelpad=5)

#ax.text(0.25, 0.87, r'$\rm <-18$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=30,color='gray')

ax.set_xlim(0.1,7)
ax.set_ylim(-2.5,0.3)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/bh_mass_density.pdf",fmt='pdf')
#plt.show()

