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

###### Brightman2012
br12={
	"zmin":np.array([0.1,1.0,2.5]),
	"zmax":np.array([1.0,2.5,4.0]),
	"f":np.array([35.5,32.0,49.5]),
	"ferr" :np.array([7.5,5.2,12.4]),
}
##########################################################################
# Lanzuisi 2018
la18={
	"zmin":np.array([0.04,1.0,2.0]),
	"zmax":np.array([1.0,2.0,3.5]),
	"f":np.array([19.,30.0,48.0]),
	"ferr" : np.array([[0.06,0.08,0.11],[0.07,0.10,0.12]])
}

####### Ueda2014
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

def fraction(logLx,z):
	if phi(logLx,z)< (1+eps)/(3+eps):
		return fCTK/2.*phi(logLx,z)*2/(fCTK/2.*phi(logLx,z)*2+1.)
	else:
		return fCTK/2.*phi(logLx,z)*2/(fCTK/2.*phi(logLx,z)*2+1.)
###########################################################################

fig=plt.figure(figsize=(15,10))
ax  = fig.add_axes([0.11,0.12,0.79,0.83])

x_fit=np.linspace(0.,5.1,1000)
data=np.genfromtxt("buchner2015.dat",names=True)
id=data["id"]==1
yup=interp1d(data["x"][id],data["y"][id])(x_fit)
id=data["id"]==0
ymid=interp1d(data["x"][id],data["y"][id])(x_fit)
id=data["id"]==2
ydown=interp1d(data["x"][id],data["y"][id])(x_fit)
ax.fill_between( x_fit ,yup,ydown,color='orange', edgecolor='white' ,alpha=0.4)
ax.plot( x_fit, ymid, c='orange', label=r'$\rm Buchner+$ $\rm 2015$ ($\log{L_{\rm X}}=43.2-43.6$)')

x_fit=np.linspace(0.005,0.78,1000)
data=np.genfromtxt("aird2015.dat",names=True)
id=data["id"]==0
yup=interp1d(data["x"][id],data["y"][id])(x_fit)
id=data["id"]==1
ymid=interp1d(data["x"][id],data["y"][id])(x_fit)
id=data["id"]==2
ydown=interp1d(data["x"][id],data["y"][id])(x_fit)
ax.fill_between( 10**x_fit-1 ,yup,ydown,color='royalblue', edgecolor='white' ,alpha=0.3)
ax.plot( 10**x_fit-1, ymid, c='royalblue', label=r'$\rm Aird+$ $\rm 2015$ ($\log{L_{\rm X}}=43.5$)')


data=np.genfromtxt("aird2015_44_5.dat",names=True)
ax.plot( 10**data["x"]-1, data["y"], '--', dashes=(15,9) ,c='royalblue', label=r'$\rm Aird+$ $\rm 2015$ ($\log{L_{\rm X}}=44.5$)')


x_fit=np.linspace(0,5,1000)
f_fit=0.0*x_fit
for i in range(len(x_fit)):
	f_fit[i]= fraction(43.5, x_fit[i])
ax.plot(x_fit, f_fit, c='crimson', lw=5, label=r'$\rm Ueda+$ $\rm 2014$ ($\bf fid.$ $\log{L_{\rm X}}=43.5$)')

for i in range(len(x_fit)):
	f_fit[i]= fraction(44.5, x_fit[i])
ax.plot(x_fit, f_fit, '--', dashes=(15,9), lw=5, c='crimson', label=r'$\rm Ueda+$ $\rm 2014$ (${\bf fid.}$ $\log{L_{\rm X}}=44.5$)')


ax.errorbar( (br12["zmin"]+br12["zmax"])/2., br12["f"]/100., yerr=br12["ferr"]/100., xerr=(br12["zmax"]-br12["zmin"])/2.,
			color='seagreen', mec='seagreen', linestyle='', marker='o', capsize=8, capthick=3, lw=4, fillstyle='none', mew=3,
			label=r'$\rm Brightman+$ $\rm 2012$ ($\log{L_{\rm X}}=43.5$)' )

ax.errorbar( 0.03, 0.20, yerr= ([0.06],[0.09]), marker='o', c='gray', mec='gray', fillstyle='none', mew=3,
			 linestyle='', capsize=8, capthick=3, lw=4, label=r'$\rm Burlon+$ $\rm 2011$ ($\log{L_{\rm X}}=43.6$)')

ax.errorbar( 2, 0.36, yerr= 0.12, xerr=1, marker='o', c='purple', mec='purple',
			 linestyle='', capsize=8, capthick=3, lw=4, label=r'$\rm Del$ $\rm Moro+$ $\rm 2016$ ($\log{L_{\rm X}}\sim44$)')

ax.errorbar( 0.055, 0.27, yerr= 0.06, lolims=True, marker='o', c='magenta',
			 linestyle='', capsize=10, mew=0, lw=4, label=r'$\rm Ricci+$ $\rm 2015$')

ax.errorbar( 1.1, 0.115, yerr= 0.06, lolims=True, marker='o', c='cyan',
			 linestyle='', capsize=10, mew=0, lw=4, label=r'$\rm Masini+$ $\rm 2018$')

ax.errorbar( (la18["zmin"]+la18["zmax"])/2., la18["f"]/100., yerr=la18["ferr"], xerr=(la18["zmax"]-la18["zmin"])/2.,
			color='navy', mec='navy', linestyle='', marker='o', capsize=8, capthick=3, lw=4,
			label=r'$\rm Lanzuisi+$ $\rm 2018$ ($\log{L_{\rm X}}=44.5$)' )


prop = matplotlib.font_manager.FontProperties(size=18.8)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=2,ncol=2)
ax.set_xlabel(r'$\rm z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$F_{\rm CTK}$',fontsize=40,labelpad=5)
ax.set_xlim(-0.05,5)
ax.set_ylim(0,0.9)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/fCTK_vs_z.pdf",fmt='pdf')