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

def lamb_vs_Lbol(Lbol):
	median = 0.469 * Lbol - 22.46
	#median = -1.5 
	sigma = 0.4
	return median, sigma

def Lbol_to_mass(Lbol):
	C = np.log10(1.26e38)
	median, sigma = lamb_vs_Lbol(Lbol)
	lamb = np.random.normal(median, sigma, size=nsample)
	return Lbol - C -lamb

mock_volume = 100. #Mpc^3
nsample = 30000
def generate_ensemble(Lbol, LogPhi):
	Masses = np.array([])
	Weights = np.array([])
	binlength = Lbol[5]-Lbol[4]
	for i in range(len(Lbol)):
		Masses = np.append(Masses, Lbol_to_mass(Lbol[i]))
		Weights = np.append(Weights, np.ones(nsample)*10**LogPhi[i]*mock_volume*binlength/nsample)
	return Masses, Weights
	
def get_mass_function(Masses, Weights):
	bins = np.linspace(5., 12., 26)
	lenbin = bins[5]-bins[4]
	result,_ = np.histogram(Masses, weights=Weights, bins=bins)
	result = result/lenbin/mock_volume
	result = np.log10(result)
	return (bins[1:]+bins[:-1])/2., result

fduty=0.03

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

i=0
colors='royalblue'
redshift = 0.2
	
source = np.genfromtxt("../../Fit_parameters/codes/zevolution_fit_global.dat",names=True)
zref = 2.
p=source['value'][ source['paraid']==0 ]
gamma1 = polynomial(redshift,p,2)
p=source['value'][ source['paraid']==1 ]
gamma2 = doublepower(redshift,(p[0],zref,p[1],p[2]))
p=source['value'][ source['paraid']==2 ]
logphi = polynomial(redshift,p,1)
p=source['value'][ source['paraid']==3 ]
Lbreak = doublepower(redshift,(p[0],zref,p[1],p[2]))
parameters_global_2 = np.array([gamma1,gamma2,logphi,Lbreak])
catalog = {"mass":0, "weight":0}
L_bol_grid_extended = np.linspace(5,15,100)
catalog["mass"], catalog["weight"] = generate_ensemble(L_bol_grid_extended+L_solar, LF(L_bol_grid_extended,parameters_global_2) )
x,y = get_mass_function(catalog["mass"],catalog["weight"])
ax.plot(x,y-np.log10(fduty),'-',c=colors,label=r'$\rm Convolved$ ($\rm total$)')
ax.plot(x,y+np.log10(0.38),'--',dashes=(25,15),c=colors,label=r'$\rm Convolved$ ($\rm type-1$)')

def BHMF1(logM, logphi_s, logM_s, alpha, beta):
	a = 10.**(logM-logM_s)
	return 10.**logphi_s * a**(1.+alpha) * np.exp(1-a)

def BHMF2(logM, logphi_s, logM_s, alpha, beta):
	a = 10.**(logM-logM_s)
	return 10.**logphi_s / (a**alpha + a**beta)
xfit = np.linspace(5,12,100)
#pfit = np.log10(BHMF2(xfit, -3.997, 8.033, 1.707, 0.569))
pfit = np.log10(BHMF2(xfit, -1.39, 6.435, 1.7, -0.8423))

ax.plot(xfit , pfit-np.log10(fduty),'-', color='crimson', label=r'$\rm Deconvolved$ ($\rm total$)')
ax.plot(xfit , pfit+np.log10(0.38),'--', dashes=(25,15), color='crimson', label=r'$\rm Deconvolved$ ($\rm type-1$)')

data = np.genfromtxt(datapath + 'BHMF/Kelly2013/table1.dat', names=True)
id = data['z']==0.4
ax.errorbar( data['logM'][id], data['logPhi'][id], yerr=( data['logPhi'][id]-data['lo'][id], data['up'][id]-data['logPhi'][id]),
 c='gray', mec='gray', linestyle='none' ,marker='o',markersize=15, capthick=4, capsize=9, label=r'$\rm Kelly+$ $\rm 2013$')

data = np.genfromtxt(datapath + 'BHMF/Vika2009_local.dat', names=True)
scale = 1e-4 * (hubble/0.70)**3
ax.errorbar( data['logM'], np.log10(data['Phi']*scale), yerr=( np.log10(data['Phi']) - np.log10(data['Phi']-data['loerr'])
	, np.log10(data['Phi']+data['uperr']) - np.log10(data['Phi']) ), c='seagreen', mec='seagreen', linestyle='none' ,
	marker='o',markersize=15, capthick=4, capsize=9, label=r'$\rm Vika+$ $\rm 2009$')

data = np.genfromtxt(datapath + 'BHMF/Shankar2009_local.dat', names=True)
ax.fill_between( data['logM'], y1= data['lo'] - data['logM'], y2= data['up'] - data['logM'], color='chocolate', 
	edgecolor='none',alpha=0.3,label=r'$\rm Shankar+$ $\rm 2009$')

data = np.genfromtxt(datapath + 'BHMF/Marconi2004_local.dat', names=True)
ax.plot( data['logM'], data['logPhi'], c='darkorchid', mec='darkorchid', linestyle='none' ,
	marker='o',markersize=15, label=r'$\rm Marconi+$ $\rm 2004$')

prop = matplotlib.font_manager.FontProperties(size=22.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\log{(M_{\rm BH}[{\rm M}_{\odot}])}$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\phi[{\rm dex}^{-1}{\rm cMpc}^{-3}])}$',fontsize=40,labelpad=5)
ax.text(0.88, 0.92, r'${\rm local}$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=40)

ax.set_xlim(5.8, 12.)
ax.set_ylim(-8,-1)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/BHMF_local.pdf",fmt='pdf')
#plt.show()

