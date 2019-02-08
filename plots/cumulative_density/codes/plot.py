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

fit_res=np.genfromtxt("../../fitresult/fit_at_z.dat",names=True)
pgamma1=fit_res["gamma1"]
pgamma2=fit_res["gamma2"]
plogphis=fit_res["phi_s"]
pLbreak=fit_res["L_s"]
pz=fit_res['z']
zpoints=np.array([0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6])

fit_evolve=np.genfromtxt("../../Fit_parameters/codes/zevolution_fit.dat",names=['gamma1','gamma2','phis','Lbreak'])
T0 = np.polynomial.chebyshev.Chebyshev((1,0,0,0))
T1 = np.polynomial.chebyshev.Chebyshev((0,1,0,0))
T2 = np.polynomial.chebyshev.Chebyshev((0,0,1,0))
T3 = np.polynomial.chebyshev.Chebyshev((0,0,0,1))
def polynomial(z,p):
        xsi=1.+z
        return p[0]*T0(xsi)+p[1]*T1(xsi)+p[2]*T2(xsi)+p[3]*T3(xsi)

def doublepower(z,p):
        xsi=1.+z
        zref=p[1]
        return 2*p[0]/(np.power(xsi/(1+zref),p[2]) + np.power(xsi/(1+zref),p[3]))
zlist=np.linspace(0,7,100)
p=fit_evolve['gamma1']
gamma1=polynomial(zlist,p)
p=fit_evolve['gamma2']
gamma2=doublepower(zlist,p)
p=fit_evolve['phis']
logphis=polynomial(zlist,p)
p=fit_evolve['Lbreak']
Lbreak=doublepower(zlist,p)

#load the shared object file
c_extenstion = CDLL(homepath+'codes/c_lib/convolve.so')
convolve_c = c_extenstion.convolve
convolve_c.restype = ctypes.POINTER(ctypes.c_double * N_bol_grid)

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

def get_model_lf(parameters,nu,magnitude=False):
	L_band = bolometric_correction(L_bol_grid,nu)
	nu_c = c_double(nu)
	input_c= np.power(10.,LF(L_bol_grid,parameters)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
	res = convolve_c(input_c,nu_c)
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
	result[i,0]=np.log10( cumulative_count(L_bol_grid+L_solar ,LF_at_z_H07(L_bol_grid,parameters_init,zlist[i],'Fiducial') ,45.5,46.5))
	result[i,1]=np.log10( cumulative_count(L_bol_grid+L_solar ,LF_at_z_H07(L_bol_grid,parameters_init,zlist[i],'Fiducial') ,46.5,47.5))
	result[i,2]=np.log10( cumulative_count(L_bol_grid+L_solar ,LF_at_z_H07(L_bol_grid,parameters_init,zlist[i],'Fiducial') ,47.5,48.5))
ax.plot(zlist,result[:,0],'--',dashes=(25,15),c='crimson',label=r'$\rm Hopkins$ $\rm 2007$')
ax.plot(zlist,result[:,1],'--',dashes=(25,15),c='crimson')
ax.plot(zlist,result[:,2],'--',dashes=(25,15),c='crimson')

result=np.zeros((len(zlist),3))
for i in range(len(zlist)):
        id=pz==zlist[i]
        result[i,0]=np.log10( cumulative_count(L_bol_grid+L_solar ,LF(L_bol_grid,[gamma1[i],gamma2[i],logphis[i],Lbreak[i]]) ,45.5,46.5))
	result[i,1]=np.log10( cumulative_count(L_bol_grid+L_solar ,LF(L_bol_grid,[gamma1[i],gamma2[i],logphis[i],Lbreak[i]]) ,46.5,47.5))
	result[i,2]=np.log10( cumulative_count(L_bol_grid+L_solar ,LF(L_bol_grid,[gamma1[i],gamma2[i],logphis[i],Lbreak[i]]) ,47.5,48.5))
ax.plot(zlist,result[:,0],'-',c='seagreen',label=r'$\rm Fit$ $\rm on$ $\rm local$ $\rm fits$')
ax.plot(zlist,result[:,1],'-',c='seagreen')
ax.plot(zlist,result[:,2],'-',c='seagreen')

result=np.zeros((len(zpoints),3))
for i in range(len(zpoints)):
	id=pz==zpoints[i]
	result[i,0]=np.log10( cumulative_count(L_bol_grid+L_solar ,LF(L_bol_grid,[pgamma1[id],pgamma2[id],plogphis[id],pLbreak[id]]) ,45.5,46.5))
	result[i,1]=np.log10( cumulative_count(L_bol_grid+L_solar ,LF(L_bol_grid,[pgamma1[id],pgamma2[id],plogphis[id],pLbreak[id]]) ,46.5,47.5))
	result[i,2]=np.log10( cumulative_count(L_bol_grid+L_solar ,LF(L_bol_grid,[pgamma1[id],pgamma2[id],plogphis[id],pLbreak[id]]) ,47.5,48.5))

ax.plot(zpoints,result[:,0],'o',c='royalblue',mec='royalblue',ms=15)
ax.plot(zpoints,result[:,1],'o',c='royalblue',mec='royalblue',ms=15)
ax.plot(zpoints,result[:,2],'o',c='royalblue',mec='royalblue',ms=15)

prop = matplotlib.font_manager.FontProperties(size=30.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=1,ncol=1,frameon=False)
ax.set_xlabel(r'$\rm z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\Phi[{\rm Mpc}^{-3}])}$',fontsize=40,labelpad=5)

ax.text(0.4, 0.33, r'$\rm 47.5-48.5$'  ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=30,color='gray')
ax.text(0.4, 0.58, r'$\rm 46.5-47.5$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=30,color='gray')
ax.text(0.4, 0.88, r'$\rm 45.5-46.5$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=30,color='gray')

ax.text(0.2, 0.1, r'$\rm Bolometric$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=40,color='navy')

ax.set_xlim(0,7)
ax.set_ylim(-10.3,-3.5)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/cumu_num_bol.pdf",fmt='pdf')
#plt.show()

