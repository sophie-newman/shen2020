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
		if (logNH>=20) and (logNH<21):
			return 1-phi(logLx,z)*(2.+eps)/(1.+eps)
		elif (logNH>=21) and (logNH<22):
			return phi(logLx,z)/(1.+eps)
		elif (logNH>=22) and (logNH<23):
			return phi(logLx,z)/(1.+eps)
		elif (logNH>=23) and (logNH<24):
			return phi(logLx,z)*eps/(1.+eps)
		elif (logNH>=24) and (logNH<26):
			return phi(logLx,z)*fCTK/2.
		else: return False
	else:
		if (logNH>=20) and (logNH<21):
			return 2./3-phi(logLx,z)*(3.+2*eps)/(3.+3*eps)
		elif (logNH>=21) and (logNH<22):
			return 1./3-phi(logLx,z)*eps/(3.+3*eps)
		elif (logNH>=22) and (logNH<23):
			return phi(logLx,z)/(1.+eps)
		elif (logNH>=23) and (logNH<24):
			return phi(logLx,z)*eps/(1.+eps)
		elif (logNH>=24) and (logNH<26):
			return phi(logLx,z)*fCTK/2.
		else: return False
###########################################################################

redshift =0.05
logLx = 43.5

fig=plt.figure(figsize=(15,10))
ax  = fig.add_axes([0.11,0.12,0.79,0.83])

x_fit=np.linspace(0.,5.1,1000)
data=np.genfromtxt("buchner2015.dat",names=True)
id=data["id"]==0
yup=interp1d(data["x"][id],data["y"][id])(x_fit)
id=data["id"]==1
ymid=interp1d(data["x"][id],data["y"][id])(x_fit)
id=data["id"]==2
ydown=interp1d(data["x"][id],data["y"][id])(x_fit)
ax.fill_between( x_fit ,yup,ydown,color='orange', edgecolor='white' ,alpha=0.4)
ax.plot( x_fit, ymid, c='orange', label=r'$\rm Buchner+$ $\rm 2015$ ($\log{L_{\rm X}}=44.0-44.4$)')

x_fit=np.linspace(0.0,0.8,1000)
data=np.genfromtxt("aird2015.dat",names=True)
id=data["id"]==0
ymid=interp1d(data["x"][id],data["y"][id])(x_fit)
id=data["id"]==1
ydown=interp1d(data["x"][id],data["y"][id])(x_fit)
id=data["id"]==2
yup=interp1d(data["x"][id],data["y"][id])(x_fit)
ax.fill_between( 10**x_fit-1 ,yup,ydown,color='royalblue', edgecolor='white' ,alpha=0.3)
ax.plot( 10**x_fit-1, ymid, c='royalblue', label=r'$\rm Aird+$ $\rm 2015$ ($\log{L_{\rm X}}=44.5$)')

x_fit=np.linspace(0,5,1000)
f_fit=0.0*x_fit
for i in range(len(x_fit)):
	f_fit[i]= fraction(44.5, x_fit[i])
ax.plot(x_fit, f_fit, lw=5, c='crimson', label=r'$\rm Ueda+$ $\rm 2014$ ($\bf fid.$ $\log{L_{\rm X}}=44.5$)')

for i in range(len(x_fit)):
	f_fit[i]= fraction(43.5, x_fit[i])
ax.plot(x_fit, f_fit, '--', dashes=(15,9), lw=5, c='crimson', label=r'$\rm Ueda+$ $\rm 2014$ ($\bf fid.$ $\log{L_{\rm X}}=43.5$)')

ax.errorbar( (ge17["zmin"]+ge17["zmax"])/2., ge17["f"], yerr=np.array([ge17["fup"]-ge17["f"],ge17["f"]-ge17["fdown"]]), 
			xerr=(ge17["zmax"]-ge17["zmin"])/2., color='seagreen', mec='seagreen', linestyle='', marker='o', 
			capsize=8, capthick=3, lw=4, label=r'$\rm Georgakakis+$ $\rm 2017$ ($\log{L_{\rm X}}=44-45$)' )

ax.errorbar( (ml14["zmin"]+ml14["zmax"])/2., ml14["f"], yerr=np.array([ml14["fup"]-ml14["f"],ml14["f"]-ml14["fdown"]]), 
			xerr=(ml14["zmax"]-ml14["zmin"])/2., color='gray', mec='gray', linestyle='', marker='o', capsize=8, capthick=3, lw=4,
			label=r'$\rm Merloni+$ $\rm 2014$ ($\log{L_{\rm X}}=44.3-44.7$)' )

ax.errorbar( (ud03["zmin"]+ud03["zmax"])/2., ud03["f"], yerr=np.array([ud03["fup"]-ud03["f"],ud03["f"]-ud03["fdown"]]), 
			xerr=(ud03["zmax"]-ud03["zmin"])/2., color='magenta', mec='magenta', linestyle='', marker='o', capsize=8, capthick=3, lw=4,
			label=r'$\rm Ueda+$ $\rm 2003$ ($\log{L_{\rm X}}=43-44.5$)' )

#ax.errorbar( 0.03, 0.20, yerr= ([0.06],[0.09]), marker='o', c='gray', mec='gray',
#			 linestyle='', capsize=8, capthick=3, lw=4, label=r'$\rm Burlon+$ $\rm 2011$ ($\log{L_{\rm X}}=43.6$)')
#ax.errorbar( 2, 0.36, yerr= 0.12, marker='o', c='purple', mec='purple',
#			 linestyle='', capsize=8, capthick=3, lw=4, label=r'$\rm Del$ $\rm Moro+$ $\rm 2015$ ($\log{L_{\rm X}}\sim44$)')

#ax.errorbar( 0.055, 0.27, yerr= 0.06, lolims=True, marker='o', c='magenta',
#			 linestyle='', capsize=10, mew=0, lw=4, label=r'$\rm Ricci+$ $\rm 2015$')

#ax.errorbar( 1.1, 0.115, yerr= 0.06, lolims=True, marker='o', c='cyan',
#			 linestyle='', capsize=10, mew=0, lw=4, label=r'$\rm Masini+$ $\rm 2018$')

ax.errorbar( (la18["zmin"][:-1]+la18["zmax"][:-1])/2., la18["f"][:-1], yerr=la18["ferr"][:,:-1], xerr=(la18["zmax"][:-1]-la18["zmin"][:-1])/2.,
			color='navy', mec='navy', linestyle='', marker='o', capsize=8, capthick=3, lw=4,
			label=r'$\rm Lanzuisi+$ $\rm 2018$ ($\log{L_{\rm X}}=44.5$)' )


prop = matplotlib.font_manager.FontProperties(size=18.8)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=2,ncol=2)
ax.set_xlabel(r'$\rm z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$f_{\rm abs}$',fontsize=40,labelpad=5)
ax.set_xlim(-0.05,5)
ax.set_ylim(0,1.1)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/fabs_vs_z.pdf",fmt='pdf')