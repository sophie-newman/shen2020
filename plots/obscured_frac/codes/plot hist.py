import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager
from scipy.interpolate import interp1d

matplotlib.style.use('classic')
matplotlib.rc('xtick.major', size=15, width=3)
matplotlib.rc('xtick.minor', size=7.5, width=3)
matplotlib.rc('ytick.major', size=15, width=3)
matplotlib.rc('ytick.minor', size=7.5, width=3)
matplotlib.rc('lines',linewidth=4,markersize=12)
matplotlib.rc('axes', linewidth=4)

###################################   Ueda2014
phi_4375_0=0.43
phi_min=0.2
phi_max=0.84
a1=0.48
beta=0.24
fCTK=1.
eps=1.7
def phi_4375(z):
	if z<2:	return phi_4375_0*np.power(1.+z,a1)
	else: return phi_4375_0*np.power(1.+2,a1)

def phi(logLx,z):
	return np.minimum( phi_max, np.maximum( phi_4375(z) - beta*(logLx-43.75), phi_min ) )

def f_nh(logNH,logLx,z):
	if phi(logLx,z)<(1.+eps)/(3.+eps):
		f1 = 1-phi(logLx,z)*(2.+eps)/(1.+eps)
		f2 = phi(logLx,z)/(1.+eps)
		f3 = phi(logLx,z)/(1.+eps)
		f4 = phi(logLx,z)*eps/(1.+eps)
		f5 = phi(logLx,z)*fCTK/2.
	else:
		f1 = 2./3-phi(logLx,z)*(3.+2*eps)/(3.+3*eps)
		f2 = 1./3-phi(logLx,z)*eps/(3.+3*eps)
		f3 = phi(logLx,z)/(1.+eps)
		f4 = phi(logLx,z)*eps/(1.+eps)
		f5 = phi(logLx,z)*fCTK/2.
	
	f1 = f1/(1.+2.*f5)
	f2 = f2/(1.+2.*f5)
	f3 = f3/(1.+2.*f5)
	f4 = f4/(1.+2.*f5)
	f5 = f5/(1.+2.*f5)
	
	results = 0.0*logNH
	select = (logNH>=20) & (logNH<21)
	results[select] = f1
	select = (logNH>=21) & (logNH<22)
	results[select] = f2
	select = (logNH>=22) & (logNH<23)
	results[select] = f3
	select = (logNH>=23) & (logNH<24)
	results[select] = f4
	select = (logNH>=24) & (logNH<26)
	results[select] = f5
	
	return results

def fraction(logLx,z):
	return (1-phi(logLx,z))/(1+phi(logLx,z)*fCTK)

'''
Ls=np.linspace(43.7,47.0,10)
f_fit=0.0*Ls
for i in range(len(Ls)):
	f_fit[i]= fraction(Ls[i], 2)
print Ls
print f_fit
exit()
'''

###########################################################################
###################################   Ueda2003
phi_44_U03=0.47
beta_U03=0.10
eps_U03=1.7
phi_max_U03= (1.+eps_U03)/(3.+eps_U03)

def phi_U03(logLx):
	return np.minimum( phi_max_U03, np.maximum( phi_44_U03 - beta_U03*(logLx-44.), 0 ) )

def f_nh_U03(logNH,logLx):
	f1 = 2.-phi_U03(logLx)*(5.+2.*eps_U03)/(1.+eps_U03)
	f2 = phi_U03(logLx)/(1.+eps_U03)
	f3 = phi_U03(logLx)*eps_U03/(1.+eps_U03)
	f4 = f3
	
	f1 = f1/(1.+f4)
	f2 = f2/(1.+f4)
	f3 = f3/(1.+f4)
	f4 = f4/(1.+f4)
	
	results = 0.0*logNH
	select = (logNH>=20)   & (logNH<20.5)
	results[select] = f1
	select = (logNH>=20.5) & (logNH<23)	 
	results[select] = f2
	select = (logNH>=23)   & (logNH<24)	 
	results[select] = f3
	select = (logNH>=24)   & (logNH<25)	 
	results[select] = f4

	return results

###########################################################################

redshift =0.05
logLx = 43.5

fig=plt.figure(figsize=(15,10))
ax  = fig.add_axes([0.11,0.12,0.79,0.83])

xgrid = np.linspace(20,26,1000)

ax.plot(xgrid, f_nh_U03(xgrid, logLx), '-', alpha=0.5, c='crimson'  , label=r'$\rm Ueda+$ $\rm 2003$ ($\rm H07$)')

ai15 = np.genfromtxt("Aird2015_hist.dat", names=True)
ax.plot( ai15["logNH"],ai15["f"], '-', alpha=0.5, color="seagreen", label=r'$\rm Aird+$ $\rm 2015$')

tr09 = np.genfromtxt("Treister09_hist.dat", names=True)
ax.plot( tr09["logNH"],tr09["f"], '-', alpha=0.5, color="chocolate", label=r'$\rm Treister+$ $\rm 2009$')

gi07 = np.genfromtxt("Gilli07_hist.dat", names=True)
ax.plot( gi07["logNH"],gi07["f"], '-', alpha=0.5, color="darkorchid", label=r'$\rm Gilli+$ $\rm 2007$')

ud14 = np.genfromtxt("Ueda2014_hist.dat", names=True)
ax.errorbar( ud14["logNH"], ud14["f"], yerr=np.array([ud14["f"]-ud14["lo"],ud14["up"]-ud14["f"]]), 
			xerr=(ud14["logNH"]-ud14["left"],ud14["right"]-ud14["logNH"]), color='navy', mec='navy', 
			linestyle='', marker='o', capsize=8, capthick=3, lw=4,
			label=r'$\rm Ueda+$ $\rm 2014$ ($\rm Swift/BAT$ $\rm sample$)' )

ax.plot(xgrid, f_nh(xgrid, logLx, redshift),    lw=6,  c='royalblue', label=r'$\rm Ueda+$ $\rm 2014$ ($\bf fid.$)')

ax.text(0.10, 0.95, r'$\rm z=0.05$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=30)
ax.text(0.14, 0.88, r'$\log{L_{\rm X}}=43.5$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=30)


prop = matplotlib.font_manager.FontProperties(size=20)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=1,ncol=1)
ax.set_xlabel(r'$\log{{\rm N}_{\rm H}}$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'${\rm N}_{\rm H}$ $\rm function$',fontsize=40,labelpad=5)
ax.set_xlim(20,26)
ax.set_ylim(0,0.58)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
#plt.show()
plt.savefig("../figs/NH_hist.pdf",fmt='pdf')