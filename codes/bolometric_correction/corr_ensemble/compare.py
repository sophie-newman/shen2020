import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
from scipy.stats import binned_statistic
from scipy.optimize import curve_fit
from data import *

def std(x):
	#return (np.percentile(x[np.isfinite(x)],84)-np.percentile(x[np.isfinite(x)],16))/2.
	return np.std(x)

def H07(logLbol,band):
	if band=='B':
		c1,k1,c2,k2 = 6.25, -0.37, 9.00, -0.012
	elif band=='HX':
		c1,k1,c2,k2 = 10.83, 0.28, 6.08, -0.020
	elif band=='SX':
		c1,k1,c2,k2 = 17.87, 0.28, 10.03, -0.020
	elif band=='IR':
		c1,k1,c2,k2 = 7.40, -0.37, 10.66, -0.014
	bolcorr = c1*np.power(10**(logLbol-10-L_solar),k1) + c2*np.power(10**(logLbol-10-L_solar),k2)
	return np.log10(bolcorr)

f = np.genfromtxt('./ensemble.dat',names=['BOL', 'HX', 'SX', 'B', 'IR'])
print f.shape

bins=np.linspace(42,48,21)
bincenter=(bins[1:]+bins[:-1])/2.
x_fit = np.linspace(39,52,100)

matplotlib.style.use('classic')
matplotlib.rc('xtick.major', size=15, width=3)
matplotlib.rc('xtick.minor', size=7.5, width=3)
matplotlib.rc('ytick.major', size=15, width=3)
matplotlib.rc('ytick.minor', size=7.5, width=3)
matplotlib.rc('lines',linewidth=4)
matplotlib.rc('axes', linewidth=4)
fig=plt.figure(figsize = (15,10))
ax = fig.add_axes([0.11,0.12,0.79,0.83])

numbers,bedges,_ = binned_statistic(f['BOL'], f['HX'], statistic="count", bins=bins)
medians,_,_ = binned_statistic(f['BOL'], f['HX'], statistic="median", bins=bins)
sigmas,_,_ = binned_statistic(f['BOL'], f['HX'], statistic=std, bins=bins)
ax.plot(bincenter, medians,'.', c='royalblue', mec='royalblue', marker='o', markersize=15)
ax.plot( x_fit, H07(x_fit, 'HX'), '--', dashes=(25,15), c='royalblue'  )

prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\log{(L_{\rm bol}[{\rm erg}\,{\rm s}^{-1}])}$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\sigma_{\rm corr}$',fontsize=40,labelpad=5)
#ax.text(0.88, 0.92, r'${\rm z\sim'+str(redshift)+'}$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=40)

ax.set_xlim(42.1, 47.8)
#ax.set_ylim(0.0,0.25)
#ax.set_yscale('log')
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
#plt.savefig("corr_sigma.pdf",fmt='pdf')
plt.show()