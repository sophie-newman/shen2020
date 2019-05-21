import numpy as np
import matplotlib.pyplot as plt 
import matplotlib
from scipy.stats import binned_statistic
from scipy.stats import norm
from scipy.optimize import curve_fit
from data import *
'''
def weighted_std(values, weights):
	average = np.average(values, weights=weights)
	variance = np.average((values-average)**2, weights=weights)
	return np.sqrt(variance)

def my_binned_statistic(x, y, w, func, bins):
	result = np.zeros(bins.shape[0]-1)
	for i in range(bins.shape[0]-1):
		id_sel = (x>bins[i]) & (x<bins[i+1])
		if len(x[id_sel])!=0:
			w_sel, y_sel = w[id_sel], y[id_sel]
			result[i] = weighted_std(y_sel, w_sel)
		else: result[i] = np.nan
	return result
'''

def std(x):
	#return (np.percentile(x[np.isfinite(x)],84)-np.percentile(x[np.isfinite(x)],16))/2.
	return np.std(x)

def sigma_H07(x,s1,beta,s2):
	return s1 * 10**( beta*(x-L_solar-9)) + s2

def new_fit_func(x,k1,k2,b):
	return k1 * (x-L_solar) + k2 * (x-L_solar)**2 + b

def new_fit_func2(x,A,B,x0,sig):
	return B+A*norm.cdf(x,loc=x0,scale=sig)

f = np.genfromtxt('./ensemble.dat',names=['BOL', 'HX', 'SX', 'B', 'IR'])
print f.shape
#f2 = np.genfromtxt('./ensemble_gaussian.dat',names=['BOL', 'HX', 'SX', 'B', 'IR'])

bins=np.linspace(38,48,31)
bincenter=(bins[1:]+bins[:-1])/2.
x_fit = np.linspace(37,50,100)

matplotlib.style.use('classic')
matplotlib.rc('xtick.major', size=15, width=3)
matplotlib.rc('xtick.minor', size=7.5, width=3)
matplotlib.rc('ytick.major', size=15, width=3)
matplotlib.rc('ytick.minor', size=7.5, width=3)
matplotlib.rc('lines',linewidth=4)
matplotlib.rc('axes', linewidth=4)
fig=plt.figure(figsize = (15,10))
ax = fig.add_axes([0.11,0.12,0.79,0.83])

#numbers,bedges,_ = binned_statistic(f['BOL'], f['HX'], statistic="count", bins=bins)
#medians,_,_ = binned_statistic(f['BOL'], f['HX'], statistic="median", bins=bins)
sigmas,_,_ = binned_statistic(f['BOL'], f['HX'], statistic=std, bins=bins)
ax.plot(bincenter,sigmas,'.', c='crimson', mec='crimson', marker='o', markersize=15)
ax.plot( x_fit, sigma_H07(x_fit, 0.06, 0.10, 0.08), '--', dashes=(25,15), c='crimson'  )

pbest,_=curve_fit(new_fit_func2, bincenter, sigmas, p0=(0.06, 0.05, 42, 2), maxfev=10000)
print pbest
ax.plot(x_fit, new_fit_func2(x_fit,*pbest),'-',c='crimson' )

######################
sigmas,_,_ = binned_statistic(f['BOL'], f['SX'], statistic=std, bins=bins)
ax.plot(bincenter,sigmas,'.', c='seagreen', mec='seagreen',  marker='o', markersize=15)
ax.plot( x_fit, sigma_H07(x_fit, 0.046, 0.10, 0.08), '--',dashes=(25,15),  c='seagreen'  )

pbest,_=curve_fit(new_fit_func2, bincenter, sigmas, p0=(0.06, 0.05, 42, 2), maxfev=10000)
print pbest
ax.plot(x_fit, new_fit_func2(x_fit,*pbest),'-',c='seagreen' )

######################
sigmas,_,_ = binned_statistic(f['BOL'], f['B'], statistic=std, bins=bins)
ax.plot(bincenter,sigmas,'.', c='royalblue', mec='royalblue', marker='o', markersize=15)
ax.plot( x_fit, sigma_H07(x_fit, 0.08, -0.25, 0.06), '--', dashes=(25,15), c='royalblue'  )

pbest,_=curve_fit(new_fit_func2, bincenter, sigmas, p0=(0.06, 0.05, 42, 2), maxfev=10000)
print pbest
ax.plot(x_fit, new_fit_func2(x_fit,*pbest),'-',c='royalblue' )

######################
sigmas,_,_ = binned_statistic(f['BOL'], f['IR'], statistic=std, bins=bins)
ax.plot(bincenter,sigmas,'.', c='chocolate', mec='chocolate', marker='o', markersize=15)
ax.plot( x_fit, sigma_H07(x_fit, 0.07, -0.17, 0.086), '--', dashes=(25,15), c='chocolate'  )

pbest,_=curve_fit(new_fit_func2, bincenter, sigmas, p0=(0.06, 0.05, 42, 2), maxfev=10000)
print pbest
ax.plot(x_fit, new_fit_func2(x_fit,*pbest),'-',c='chocolate' )

######################
ax.plot([],[],'.', c='k', mec='k', marker='o', markersize=15, label=r'$\rm Binned$ $\rm estimation$')
ax.plot([],[], '--', dashes=(25,15), c='k', label=r'$\rm Hopkins+$ $\rm 2007$')
ax.plot([],[], '-', c='k', label=r'$\rm Best-fit$')

prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=1,ncol=1,frameon=False)
ax.set_xlabel(r'$\log{(L_{\rm bol}[{\rm erg}\,{\rm s}^{-1}])}$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\sigma_{\rm corr}$',fontsize=40,labelpad=5)
#ax.text(0.88, 0.92, r'${\rm z\sim'+str(redshift)+'}$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=40)

ax.set_xlim(38.1, 47.8)
#ax.set_xlim(39, 50)
ax.set_ylim(0.0,0.5)
#ax.set_yscale('log')
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("corr_sigma.pdf",fmt='pdf')
#plt.show()