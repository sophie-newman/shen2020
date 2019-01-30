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
def polynomial(z,p):
	xsi=1.+z
	return p[0]*T0(xsi)+p[1]*T1(xsi)+p[2]*T2(xsi)+p[3]*T3(xsi)

def doublepower(z,p):
	xsi=1.+z
	zref=p[1]
	return 2*p[0]/(np.power(xsi/(1+zref),p[2]) + np.power(xsi/(1+zref),p[3]))

def bestfit(z,field):
	source=np.genfromtxt("zevolution_fit.dat",names=['gamma1','gamma2','phi_s','Lbreak'])
	p=source[field]
	if (field=='gamma1') or (field=='phi_s'):
		return polynomial(z,p)
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

data=np.genfromtxt("../../fitresult/fit_at_z.dat",names=True)
ax.errorbar(data["z"]+1,data["gamma1"],yerr=data['err1'],xerr=0.25*np.ones(len(data["z"])),linestyle='',marker='o',c='royalblue',mec='royalblue',ms=15,capsize=9,capthick=3,label=r'$\rm This$ $\rm work$')

ax.plot(z_a+1,gamma1_a,'--',dashes=(25,15),c='crimson',label=r'$\rm Hopkins+$ $\rm 2007$')
ax.plot(z_a+1,bestfit(z_a,'gamma1'),'-',c='seagreen',label=r'$\rm This$ $\rm work$')

prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=2,ncol=1,frameon=False)
ax.set_xlabel(r'$\rm 1+z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\rm Faint-end$ $\rm Slope$ $\rm \gamma_{\rm 1}$',fontsize=40,labelpad=5)

ax.set_xlim(1,8.)
ax.set_ylim(0.1,1.5)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/gamma1.pdf",fmt='pdf')

fig=plt.figure(figsize = (15,10))
ax = fig.add_axes([0.11,0.12,0.79,0.83])

data=np.genfromtxt("../../fitresult/fit_at_z.dat",names=True)
ax.errorbar(data["z"]+1,data["gamma2"],yerr=data['err2'],xerr=0.25*np.ones(len(data["z"])),linestyle='',marker='o',c='royalblue',mec='royalblue',ms=15,capsize=9,capthick=3,label=r'$\rm This$ $\rm work$')

ax.plot(z_a+1,gamma2_a,'--',dashes=(25,15),c='crimson')
ax.plot(z_a+1,bestfit(z_a,'gamma2'),'-',c='seagreen')

#prop = matplotlib.font_manager.FontProperties(size=25.0)
#ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\rm 1+z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\rm Bright-end$ $\rm Slope$ $\rm \gamma_{\rm 2}$',fontsize=40,labelpad=5)

ax.set_xlim(1,8.)
ax.set_ylim(0.9,3.2)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/gamma2.pdf",fmt='pdf')

fig=plt.figure(figsize = (15,10))
ax = fig.add_axes([0.11,0.12,0.79,0.83])

data=np.genfromtxt("../../fitresult/fit_at_z.dat",names=True)
ax.errorbar(data["z"]+1,data["phi_s"],yerr=data['err3'],xerr=0.25*np.ones(len(data["z"])),linestyle='',marker='o',c='royalblue',mec='royalblue',ms=15,capsize=9,capthick=3,label=r'$\rm This$ $\rm work$')

ax.plot(z_a+1,phi_s_a,'--',dashes=(25,15),c='crimson')
ax.plot(z_a+1,bestfit(z_a,'phi_s'),'-',c='seagreen')

#prop = matplotlib.font_manager.FontProperties(size=25.0)
#ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\rm 1+z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\phi_{\ast})}$ $[\rm Mpc]$',fontsize=40,labelpad=5)

ax.set_xlim(1,8.)
ax.set_ylim(-6.3,-4.1)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/phi_s.pdf",fmt='pdf')

fig=plt.figure(figsize = (15,10))
ax = fig.add_axes([0.11,0.12,0.79,0.83])

data=np.genfromtxt("../../fitresult/fit_at_z.dat",names=True)
ax.errorbar(data["z"]+1,data["L_s"],yerr=data['err4'],xerr=0.25*np.ones(len(data["z"])),linestyle='',marker='o',c='royalblue',mec='royalblue',ms=15,capsize=9,capthick=3,label=r'$\rm This$ $\rm work$')

ax.plot(z_a+1,Lbreak_a,'--',dashes=(25,15),c='crimson')
ax.plot(z_a+1,bestfit(z_a,'Lbreak'),'-',c='seagreen')

#prop = matplotlib.font_manager.FontProperties(size=25.0)
#ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\rm 1+z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(L_{\ast})}$ $[{\rm L}_{\odot}]$',fontsize=40,labelpad=5)

ax.set_xlim(1,8.)
ax.set_ylim(11.3,14.0)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/Lbreak.pdf",fmt='pdf')
