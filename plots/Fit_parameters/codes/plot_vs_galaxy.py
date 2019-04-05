from data import *
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

T0 = np.polynomial.chebyshev.Chebyshev((1,0,0,0))
T1 = np.polynomial.chebyshev.Chebyshev((0,1,0,0))
T2 = np.polynomial.chebyshev.Chebyshev((0,0,1,0))
T3 = np.polynomial.chebyshev.Chebyshev((0,0,0,1))
def polynomial(z,p,n=3):
	xsi=1.+z
	if n==1: return p[0]*T0(xsi)+p[1]*T1(xsi)
	elif n==2: return p[0]*T0(xsi)+p[1]*T1(xsi)+p[2]*T2(xsi)
	elif n==3:
		return p[0]*T0(xsi)+p[1]*T1(xsi)+p[2]*T2(xsi)+p[3]*T3(xsi)
	else: return False

def doublepower(z,p):
	xsi=1.+z
	zref=p[1]
	return 2*p[0]/(np.power(xsi/(1+zref),p[2]) + np.power(xsi/(1+zref),p[3]))

def bestfit(z,field):
	source=np.genfromtxt("zevolution_fit.dat",names=['gamma1','gamma2','phi_s','Lbreak'])
	p=source[field]
	if (field=='gamma1'): 
		return polynomial(z,p,2)
	elif (field=='phi_s'):
		return polynomial(z,p,1)
	else: return doublepower(z,p)

def bestfit_global(z,paraid):
	source=np.genfromtxt("zevolution_fit_global.dat",names=True)
	p=source['value'][ source['paraid']==paraid ]
	print p
	if (paraid==0): 
		return polynomial(z,p,2)
	elif (paraid==2):
		return polynomial(z,p,1)
	else: return doublepower(z,p)

def Hopkins07(z):
	parameters_init = np.array([0.41698725, 2.17443860, -4.82506430, 13.03575300, 0.63150872, -11.76356000, -14.24983300, -0.62298947, 1.45993930, -0.79280099])
	xsi_log	= np.log10((1.+z)/(1.+2.))
	gamma1_0	= parameters_init[0]	#faint-end slope
	gamma2_0	= parameters_init[1] 	#bright-end slope
	P0			= parameters_init[2]	#normalization in log
	L0			= parameters_init[3]	#break in log
	k1, k2, k3 = parameters_init[4], parameters_init[5], parameters_init[6]
	k_gamma1 = parameters_init[7]
	k_gamma2_1 = parameters_init[8]
	k_gamma2_2 = parameters_init[9]
	gamma1   = gamma1_0 * np.power(10., k_gamma1*xsi_log)
	gamma2   = 2.*gamma2_0 / (np.power(10., xsi_log*k_gamma2_1) + np.power(10., xsi_log*k_gamma2_2))
	Lbreak  = L0 + k1*xsi_log + k2*xsi_log**2 + k3*xsi_log**3
	return gamma1,gamma2,P0*np.ones(len(z)),Lbreak
z_a=np.linspace(0.01,8,1000)
gamma1_a, gamma2_a, phi_s_a, Lbreak_a = Hopkins07(z_a)

fig=plt.figure(figsize = (15,10))
ax = fig.add_axes([0.11,0.12,0.79,0.83])

data=np.genfromtxt("../galaxyLFdata/finkelstein2016.dat",names=True)
ax.errorbar(data["z"],-data["alpha"]-1,yerr=(data['alphalowerr'],data['alphauperr']),capsize=10,linestyle='',c='crimson',mec='crimson',marker='o', ms=15,capthick=4,alpha=0.6)
data=np.genfromtxt("../galaxyLFdata/bouwens2015.dat",names=True)
ax.errorbar(data["z"],-data["alpha"]-1,yerr=data['alphaerr'],linestyle='',c='seagreen',mec='seagreen',marker='o', ms=15,lw=4,capsize=10,capthick=4,alpha=0.6)
data=np.genfromtxt("../galaxyLFdata/atek2015.dat",names=True)
ax.errorbar(data["z"],-data["alpha"]-1,yerr=([data['alphalowerr']],[data['alphauperr']]),linestyle='',c='saddlebrown',mec='saddlebrown',marker='o', ms=15,lw=4,capsize=10,capthick=4,alpha=0.6)
data=np.genfromtxt("../galaxyLFdata/atek2018.dat",names=True)
ax.errorbar(data["z"],-data["alpha"]-1,yerr=([data['alphalowerr']],[data['alphauperr']]),linestyle='',c='saddlebrown',mec='saddlebrown',marker='o', ms=15,lw=4,capsize=10,capthick=4,alpha=0.6)
data=np.genfromtxt("../galaxyLFdata/bowler2015.dat",names=True)
ax.errorbar(data["z"],-data["alpha"]-1,yerr=([data['alphalowerr']],[data['alphauperr']]),linestyle='',c='yellow',mec='yellow',marker='o', ms=15,lw=4,capsize=10,capthick=4,alpha=0.6)
data=np.genfromtxt("../galaxyLFdata/ish2018.dat",names=True)
ax.errorbar(data["z"],-data["alpha"]-1,yerr=(data['alphalowerr'],data['alphauperr']),linestyle='',c='darkcyan',mec='darkcyan',marker='o', ms=15,lw=4,capsize=10,capthick=4,alpha=0.6)
data=np.genfromtxt("../galaxyLFdata/parsa2016.dat",names=True)
ax.errorbar(data["z"],-data["alpha"]-1,yerr=data['alphaerr'],linestyle='',c='silver',mec='silver',marker='o', ms=15,lw=4,capsize=10,capthick=4,alpha=0.6)
data=np.genfromtxt("../galaxyLFdata/duncan2014.dat",names=True)
ax.errorbar(data["z"],-data["alpha"]-1,yerr=(data['alphalowerr'],data['alphauperr']),linestyle='',c='purple',mec='purple',marker='o', ms=15,lw=4,capsize=10,capthick=4,alpha=0.6)
data=np.genfromtxt("../galaxyLFdata/metha2017.dat",names=True)
ax.errorbar(data["z"],-data["alpha"]-1,yerr=(data['alphalowerr'],data['alphauperr']),linestyle='',c='pink',mec='pink',marker='o', ms=15,lw=4,capsize=10,capthick=4,alpha=0.6)

'''
data=np.genfromtxt("../galaxyLFdata/yung2018.dat",names=True)
ax.errorbar(data["z"],-data["alpha"]-1,linestyle='',c='chocolate',mec='chocolate',marker='s', ms=15,lw=4,capsize=10,alpha=0.6)
data=np.genfromtxt("../galaxyLFdata/mason2015.dat",names=True)
ax.errorbar(data["z"],-data["alpha"]-1,yerr=data['alphaerr'],linestyle='',c='darkorchid',mec='darkorchid',marker='s', ms=15,lw=4,capsize=10,capthick=4,alpha=0.6)
data=np.genfromtxt("../galaxyLFdata/tachella2013.dat",names=True)
ax.errorbar(data["z"],-data["alpha"]-1,yerr=(data['alphalowerr'],data['alphauperr']),linestyle='',c='royalblue',mec='royalblue',marker='s', ms=15,lw=4,capsize=10,capthick=4,alpha=0.6)
data=np.genfromtxt("../galaxyLFdata/tachella2018.dat",names=True)
ax.errorbar(data["z"],-data["alpha"]-1,yerr=(data['alphalowerr'],data['alphauperr']),linestyle='',c='slateblue',mec='slateblue',marker='s', ms=15,lw=4,capsize=10,capthick=4,alpha=0.6)
data=np.genfromtxt("../galaxyLFdata/wilkins2017.dat",names=True)
ax.errorbar(data["z"],-data["alpha"]-1,linestyle='',c='olive',mec='olive',marker='s', ms=15,lw=4,capsize=10,capthick=4,alpha=0.6)
data=np.genfromtxt("../galaxyLFdata/jaacks2012.dat",names=True)
ax.errorbar(data["z"],-data["alpha"]-1,yerr=(data['alphalowerr'],data['alphauperr']),linestyle='',c='tan',mec='tan',marker='s', ms=15,lw=4,capsize=10,capthick=4,alpha=0.6)
'''

data=np.genfromtxt("../../fitresult/fit_at_z.dat",names=True)
ax.errorbar(data["z"]+1,data["gamma1"],yerr=data['err1'],linestyle='',marker='o',
	c='black',mec='black',ms=18,capsize=10,capthick=4,label=r'$\rm Local$ $\rm fits$ ($\phi_{\ast}$ $\rm fixed$)')

ax.plot(z_a+1,bestfit_global(z_a,0),'-',c='gray',label=r'$\rm Global$ $\rm fit$')

prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=2,ncol=1,frameon=False)
ax.set_xlabel(r'$\rm 1+z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\rm Faint-end$ $\rm Slope$ $\rm \gamma_{\rm 1}$',fontsize=40,labelpad=5)

ax.set_xlim(1,9.2)
ax.set_ylim(-0.1,1.9)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.show()
#plt.savefig("../figs/gamma1_vs_galaxy.pdf",fmt='pdf')
