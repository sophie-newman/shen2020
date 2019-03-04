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

def Gamma(epsilon,z):
	alphaEUV=-1.70
	return ( 2.*(1+z)**(-1.5)*epsilon/1e24/(3.+np.abs(alphaEUV)) )

parameters_init = np.array([0.41698725, 2.17443860, -4.82506430, 13.03575300, 0.63150872, -11.76356000, -14.24983300, -0.62298947, 1.45993930, -0.79280099])

fit_res=np.genfromtxt("../../fitresult/fit_at_z.dat",names=True)
pgamma1, pgamma1_err  = fit_res["gamma1"], fit_res["err1"]
pgamma2, pgamma2_err  = fit_res["gamma2"], fit_res["err2"]
plogphis,plogphis_err = fit_res["phi_s"], fit_res["err3"]
pLbreak, pLbreak_err  = fit_res["L_s"], fit_res["err3"]
pz=fit_res['z']
#zpoints=np.array([0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6])
zpoints=np.array(pz)

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
zlist=np.linspace(0.4,7,100)
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

def get_model_lf(parameters,nu,magnitude=False,global_model=False):
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

def cumulative_emissivity(L_band,Phi_band,L_limit_low,L_limit_up):
	Mband, Mlow, Mup = L_band, L_limit_low, L_limit_up
	logphi=inter.interp1d(Mband,Phi_band)
	def emis(x):
		fnu = np.power(10.,-0.4*x)*3631*1e-23*4*np.pi*(10*con.pc.value*100)**2
		fnu912 = fnu * (912./1450.)**(0.61)
		return np.power(10.,logphi(x))*fnu
	return quad(emis,Mlow,Mup)[0]

def Gamma_err(parameters,errs,L_limit_low,L_limit_up,redshift):
	partials = 0.0 * parameters
	delta = 1e-6
	def fobjective(parameters):
		M_1450, PHI_1450 = get_model_lf(parameters, -1, magnitude=True)
        	return Gamma(cumulative_emissivity(M_1450, PHI_1450, L_limit_low, L_limit_up),redshift)
	for i in range(len(parameters)):
		parameters_add = parameters.copy()
		parameters_add[i] += delta
		partials[i] = ( fobjective(parameters_add) - fobjective(parameters))/delta
		partials[i] = np.abs(partials[i])

	print np.ravel(partials)
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
ax.plot(zlist,np.log10(result[:,0]),'--',dashes=(25,15),c='crimson',label=r'$\rm Hopkins$ $\rm 2007$')
#ax.plot(zlist,result[:,1],'--',dashes=(25,15),c='crimson')

result=np.zeros((len(zlist),2))
for i in range(len(zlist)):
	L_band = bolometric_correction(L_bol_grid,-1)
	M_1450 = (M_sun_Bband_AB -2.5*L_band) + 0.706
	PHI_1450 = return_kk18_lf_fitted(M_1450 ,zlist[i]) 
	result[i,0]= Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, -18),zlist[i])
	result[i,1]= Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, -21),zlist[i])
ax.plot(zlist,np.log10(result[:,0]),'--',dashes=(25,15),c='darkorchid',label=r'$\rm Kulkarni$ $\rm 2018$')
#ax.plot(zlist,result[:,1],'--',dashes=(25,15),c='cyan')

result=np.zeros((len(zlist),2))
for i in range(len(zlist)):
	id=pz==zlist[i]
	M_1450, PHI_1450 = get_model_lf([gamma1[i],gamma2[i],logphis[i],Lbreak[i]], -1, magnitude=True)
	result[i,0]= Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, -18),zlist[i])
	result[i,1]= Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, -21),zlist[i])
ax.plot(zlist,np.log10(result[:,0]),'-',c='seagreen',label=r'$\rm Fit$ $\rm on$ $\rm local$ $\rm fits$')
#ax.plot(zlist,result[:,1],'-',c='seagreen')

result=np.zeros((len(zpoints),2))
uncertainty=np.zeros((len(zpoints),2))
for i in range(len(zpoints)):
        id=pz==zpoints[i]
        M_1450, PHI_1450 = get_model_lf([pgamma1[id],pgamma2[id],plogphis[id],pLbreak[id]], -1, magnitude=True)
        result[i,0]= Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, -18),zpoints[i])
	uncertainty[i,0]= Gamma_err(np.array([pgamma1[id],pgamma2[id],plogphis[id],pLbreak[id]]), [pgamma1_err[id],pgamma2_err[id],plogphis_err[id],pLbreak_err[id]], lowlimit, -18, zpoints[i])
	result[i,1]= Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, -21),zpoints[i])
print uncertainty[:,0]
uperr = np.log10(result + uncertainty)-np.log10(result)
loerr = np.log10(result)-np.log10(result - uncertainty)
loerr[np.invert(np.isfinite(loerr))] = 100
ax.plot(zpoints,np.log10(result[:,0]),'o',c='royalblue',mec='royalblue',ms=15)
ax.errorbar(zpoints,np.log10(result[:,0]),yerr=(loerr[:,0],uperr[:,0]),linestyle='none',marker='o',c='royalblue',mec='royalblue',ms=15,capsize=9,capthick=4)
#ax.plot(zpoints,result[:,1],'o',c='royalblue',mec='royalblue',ms=15)

ax.errorbar([2.40,2.80,3.20,3.60,4.00,4.40,4.75],[0.015,-0.066,-0.103,-0.097,-0.072,-0.019,-0.029],yerr=([0.146, 0.131, 0.121, 0.118, 0.117, 0.122, 0.147],[0.132, 0.129, 0.130, 0.131, 0.135, 0.140, 0.156]),marker='s',linestyle='none',ms=15,color='k',mec='k',capsize=9,capthick=4,label=r'$\rm Becker$ $\rm &$ $\rm Bolton+$ $\rm 2013$')

ydata=np.array([0.58, 0.53, 0.48, 0.47, 0.45, 0.29])
lowerr=np.array([0.20, 0.19, 0.18, 0.18, 0.17, 0.11])
uperr=np.array([0.08, 0.09, 0.10, 0.12, 0.14, 0.11])
ax.errorbar([4.8,5.0,5.2,5.4,5.6,5.8], np.log10(ydata) ,yerr=(np.log10(ydata)-np.log10(ydata-lowerr),np.log10(ydata+uperr)-np.log10(ydata)),marker='s',linestyle='none',ms=15,color='navy',mec='navy',capsize=9,capthick=4,label=r'$\rm Aloisio+$ $\rm 2018$')

#data=np.genfromtxt("kk18.dat",names=['z','gamma'])
#ax.plot(data['z'],data['gamma'],'--',c='gray',alpha=0.3)

prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\rm z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\Gamma_{\rm -12}$ [$\rm s^{-1}\,atom^{-1}$]',fontsize=40,labelpad=5)

#ax.text(0.25, 0.64, r'$\rm <-21$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=30,color='gray')
#ax.text(0.25, 0.87, r'$\rm <-18$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=30,color='gray')

ax.set_xlim(0.3,7)
ax.set_ylim(-2.5,0.3)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/ionizing_photon.pdf",fmt='pdf')
#plt.show()

