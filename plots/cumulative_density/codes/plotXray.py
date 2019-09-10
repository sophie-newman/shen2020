from data import *
import numpy as np 
from lf_shape import *
import scipy.interpolate as inter
from scipy.integrate import quad
from convolve import *
from convolve_h07 import *
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import least_squares
import astropy.constants as con
# fit the luminosity function based on datasets at a given redshift
from ctypes import *
import ctypes
import sys

parameters_init = np.array([0.41698725, 2.17443860, -4.82506430, 13.03575300, 0.63150872, -11.76356000, -14.24983300, -0.62298947, 1.45993930, -0.79280099])

data=np.genfromtxt("../../../codes/lf_fit/output/fit_at_z_fix.dat",names=True)
pgamma1, pgamma1_err  = data["gamma1"], data["err1"]
pgamma2, pgamma2_err  = data["gamma2"], data["err2"]
plogphis,plogphis_err = data["phi_s"],  data["err3"]
pLbreak, pLbreak_err  = data["L_s"],    data["err3"]
pz=data['z']
zpoints=np.array(pz)

fit_evolve=np.genfromtxt("../../Fit_parameters/codes/zevolution_fit_global.dat",names=True)
paraid, pglobal, pglobal_err = fit_evolve['paraid'], fit_evolve['value'], (fit_evolve['uperr']+fit_evolve['loerr'])/2.
zlist=np.linspace(0.1,7,100)

#load the shared object file
c_extenstion = CDLL(homepath+'codes/c_lib/convolve.so')
convolve_c = c_extenstion.convolve
convolve_c.restype = ctypes.POINTER(ctypes.c_double * N_bol_grid)

#
c_extenstion_old = CDLL(homepath+'codes/c_lib/convolve_old.so')
convolve_c_old = c_extenstion_old.convolve
convolve_c_old.restype = ctypes.POINTER(ctypes.c_double * N_bol_grid)

def cumulative_count(L_band,Phi_band,L_limit_low,L_limit_up,ABmag=False):
	if ABmag==False:
		logphi=inter.interp1d(L_band,Phi_band)
		def phi(x):
			return np.power(10.,logphi(x))
		return quad(phi,L_limit_low,L_limit_up)[0]
	else:
		Mband, Mlow, Mup = L_band, L_limit_low, L_limit_up
		logphi=inter.interp1d(Mband,Phi_band)
		def phi(x):
			return np.power(10.,logphi(x))
		return quad(phi,Mup,Mlow)[0]

def cumulative_lum(L_band,Phi_band,L_limit_low,L_limit_up,ABmag=False):
	if ABmag==False:
		logphi=inter.interp1d(L_band,Phi_band)
		def lum(x):
			return np.power(10.,logphi(x))*np.power(10.,x)
		return quad(lum,L_limit_low,L_limit_up)[0]
	else:
		Mband, Mlow, Mup = L_band, L_limit_low, L_limit_up
		logphi=inter.interp1d(Mband,Phi_band)
		def lum(x):
			nuLnu = np.power(10.,-0.4*x)*3631*1e-23*4*np.pi(10*con.pc.value*100)**2 * con.c.value/(1450.*1e-10)
			return np.power(10.,logphi(x))*nuLnu
		return quad(lum,Mup,Mlow)[0]

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

def get_model_lf_global(parameters,nu,redshift,magnitude=False):
	zref = 2.
        p=parameters[paraid==0]
        gamma1 = polynomial(redshift,p,2)
        p=parameters[paraid==1]
        gamma2 = doublepower(redshift,(p[0],zref,p[1],p[2]))
        p=parameters[paraid==2]
        logphi = polynomial(redshift,p,1)
        p=parameters[paraid==3]
        Lbreak = doublepower(redshift,(p[0],zref,p[1],p[2]))
        parameters_at_z = np.array([gamma1,gamma2,logphi,Lbreak])
        return get_model_lf(parameters_at_z,nu,redshift,magnitude=magnitude)

def get_model_lf(parameters,nu,redshift,magnitude=False):
	L_band = bolometric_correction(L_bol_grid,nu)
	nu_c = c_double(nu)
	redshift_c = c_double(redshift)
	dtg_c = c_double(return_dtg(redshift))
	input_c= np.power(10.,LF(L_bol_grid,parameters)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
	res = convolve_c(input_c,nu_c,redshift_c,dtg_c)
	res = [i for i in res.contents]
	PHI_band = np.array(res,dtype=np.float64)
	if magnitude==False:
		return L_band, PHI_band
	else:
		M_1450 = (M_sun_Bband_AB -2.5*L_band) + 0.706
		PHI_1450 = np.log10(PHI_band) - np.log10(2.5)
		return M_1450, PHI_1450

result=np.zeros((len(zlist),3))
for i in range(len(zlist)):
	L_band = bolometric_correction_old(L_bol_grid,-4)
        nu_c = c_double(-4)
        input_c= np.power(10.,LF_at_z_H07(L_bol_grid,parameters_init,zlist[i],"Fiducial")).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        res = convolve_c_old(input_c,nu_c)
        res = [j for j in res.contents]
        PHI_band = np.array(res,dtype=np.float64)

	result[i,0]=np.log10( cumulative_count(L_band+L_solar, np.log10(PHI_band),43.,44.))
	result[i,1]=np.log10( cumulative_count(L_band+L_solar, np.log10(PHI_band),44.,45.))
	result[i,2]=np.log10( cumulative_count(L_band+L_solar, np.log10(PHI_band),45.,46.))
ax.plot(zlist,result[:,0],'--',dashes=(25,15),c='crimson',label=r'$\rm Hopkins$ $\rm 2007$')
ax.plot(zlist,result[:,1],'--',dashes=(25,15),c='crimson')
ax.plot(zlist,result[:,2],'--',dashes=(25,15),c='crimson')

result=np.zeros((len(zlist),3))
for i in range(len(zlist)):
	L_HX, PHI_HX = get_model_lf_global(pglobal, -4, zlist[i])
	result[i,0]=np.log10(cumulative_count(L_HX+L_solar,np.log10(PHI_HX),43.,44.))
	result[i,1]=np.log10(cumulative_count(L_HX+L_solar,np.log10(PHI_HX),44.,45.))
	result[i,2]=np.log10(cumulative_count(L_HX+L_solar,np.log10(PHI_HX),45.,46.))
ax.plot(zlist,result[:,0],'-',c='darkorchid')
ax.plot(zlist,result[:,1],'-',c='darkorchid')
ax.plot(zlist,result[:,2],'-',c='darkorchid')

'''
result=np.zeros((len(zpoints),3))
for i in range(len(zpoints)):
        id=pz==zpoints[i]
        L_HX, PHI_HX = get_model_lf([pgamma1[id],pgamma2[id],plogphis[id],pLbreak[id]], -4)
	result[i,0]=np.log10( cumulative_count(L_HX+L_solar,np.log10(PHI_HX),43.5,44.5))
	result[i,1]=np.log10( cumulative_count(L_HX+L_solar,np.log10(PHI_HX),44.5,45.5))
        result[i,2]=np.log10( cumulative_count(L_HX+L_solar,np.log10(PHI_HX),45.5,46.5))
ax.plot(zpoints,result[:,0],'o',c='royalblue',mec='royalblue',ms=15)
ax.plot(zpoints,result[:,1],'o',c='royalblue',mec='royalblue',ms=15)
ax.plot(zpoints,result[:,2],'o',c='royalblue',mec='royalblue',ms=15)
'''
data = np.genfromtxt("../obdata/aird15.dat",names=True)
ax.plot( 10**data["logzplusone"][data["id"]==1]-1, data["logPhi"][data["id"]==1], '--', dashes=(25,10), c='seagreen', label=r'$\rm Aird+$ $\rm 2015$')
ax.plot( 10**data["logzplusone"][data["id"]==2]-1, data["logPhi"][data["id"]==2], '--', dashes=(25,10), c='seagreen')
ax.plot( 10**data["logzplusone"][data["id"]==3]-1, data["logPhi"][data["id"]==3], '--', dashes=(25,10), c='seagreen')

data = np.genfromtxt("../obdata/miyaji2015.dat",names=True)
ax.errorbar( data["zplusone"]-1, data["logPhi"], yerr=(data["logPhi"]-data["lo"],data["up"]-data["logPhi"]), c='royalblue', mec="royalblue", linestyle='none', marker='o', ms=15, capthick=4, capsize=0, label=r'$\rm Miyaji+$ $\rm 2015$')

data = np.genfromtxt("../obdata/ueda14.dat",names=True)
ax.errorbar( data["z"], data["logPhi"], yerr=(data["logPhi"]-data["lo"],data["up"]-data["logPhi"]), c='k', mec="k", linestyle='none', marker='o', ms=15, capthick=4, capsize=0, label=r'$\rm Ueda+$ $\rm 2014$')

prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=1,ncol=1,frameon=False)
ax.set_xlabel(r'$\rm z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\Phi[{\rm Mpc}^{-3}])}$',fontsize=40,labelpad=5)

ax.text(0.32, 0.3, r'$\rm 45-46$'  ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=30,color='gray')
ax.text(0.32, 0.64, r'$\rm 44-45$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=30,color='gray')
ax.text(0.32, 0.9, r'$\rm 43-44$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=30,color='gray')

ax.text(0.2, 0.1, r'$\rm Hard$ $\rm X-ray$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=40,color='navy')

ax.set_xlim(0,7)
ax.set_ylim(-9.6,-2.4)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/cumu_num_Xray.pdf",fmt='pdf')
#plt.show()

