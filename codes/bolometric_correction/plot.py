from data import *
import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit 
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

def sigma_H07(x,s1,beta,s2):
	return s1 * 10**( beta*(x-L_solar-9)) + s2

def sigma_new(x,A,B,x0,sig):
	return B+A*norm.cdf(x,loc=x0,scale=sig)

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

def fit_func(logLbol,c1,k1,c2,k2):
	bolcorr = c1*np.power(10**(logLbol-10-L_solar),k1) + c2*np.power(10**(logLbol-10-L_solar),k2)
	return np.log10(bolcorr)

fig=plt.figure(figsize = (15,10))
ax = fig.add_axes([0.11,0.12,0.79,0.83])

data = np.genfromtxt("bolcorr.dat",names=['Lbol','LHX','LSX','LB','L1450','LIR'])

ax.plot(data['Lbol'], data['Lbol']-data['LB']   , lw=6, c='royalblue' ,label=r'$\rm B$ $\rm band$')
ax.plot(data['Lbol'], data['Lbol']-data['L1450'], lw=6, c='darkorchid',label=r'$\rm UV$ $\rm 1450\AA$')
ax.plot(data['Lbol'], data['Lbol']-data['LIR']  , lw=6, c='crimson'   ,label=r'$\rm Mid$ $\rm IR$ ($\rm 15\mu m$)')
ax.plot(data['Lbol'], data['Lbol']-data['LHX']  , lw=6, c='seagreen'  ,label=r'$\rm Hard$ $\rm X-ray$')
ax.plot(data['Lbol'], data['Lbol']-data['LSX']  , lw=6, c='chocolate' ,label=r'$\rm Soft$ $\rm X-ray$')

#data2 = np.genfromtxt("bolcorr_2.dat",names=['Lbol','LHX','LSX','LB','L1450','LIR'])
#ax.plot(data2['Lbol'], data2['Lbol']-data2['LHX'] , lw=6, c='cyan')

xs, ys = data['Lbol'], data['Lbol']-data['LB']
sig = sigma_new(xs,-0.38269949 , 0.4052673 , 42.38666388 , 2.3775969)
ax.fill_between(xs, y1=ys-sig, y2=ys+sig, color='royalblue', edgecolor='none', alpha=0.2 )
args,_=curve_fit(fit_func, xs, ys)
print 'B:', args

xs, ys = data['Lbol'], data['Lbol']-data['L1450']
sig = sigma_new(xs,-0.37199551 , 0.40486931, 42.30731163 , 2.30978248)
ax.fill_between(xs, y1=ys-sig, y2=ys+sig, color='darkorchid', edgecolor='none', alpha=0.2 )
args,_=curve_fit(fit_func, xs, ys)
print '1450:', args

xs, ys = data['Lbol'], data['Lbol']-data['LIR']
sig = sigma_new(xs,-0.33807144 , 0.40716262, 42.15882924 , 2.19283454)
ax.fill_between(xs, y1=ys-sig, y2=ys+sig, color='crimson', edgecolor='none', alpha=0.2 )
args,_=curve_fit(fit_func, xs, ys)
print 'IR:', args

xs, ys = data['Lbol'], data['Lbol']-data['LHX']
sig = sigma_new(xs,0.19262562 , 0.06592312, 42.9876632  , 1.88296393)
ax.fill_between(xs, y1=ys-sig, y2=ys+sig, color='seagreen', edgecolor='none', alpha=0.2 )
args,_=curve_fit(fit_func, xs, ys)
print 'HX:', args

xs, ys = data['Lbol'], data['Lbol']-data['LSX']
sig = sigma_new(xs,0.07969176 , 0.18037277, 44.16411558 , 1.49648227)
ax.fill_between(xs, y1=ys-sig, y2=ys+sig, color='chocolate', edgecolor='none', alpha=0.2 )
args,_=curve_fit(fit_func, xs, ys)
print 'SX:', args

ax.plot(data['Lbol'],H07(data['Lbol'],'B'),'--',dashes=(25,15), c='royalblue')
ax.plot(data['Lbol'],H07(data['Lbol'],'IR'),'--',dashes=(25,15), c='crimson')
ax.plot(data['Lbol'],H07(data['Lbol'],'HX'),'--',dashes=(25,15), c='seagreen')
ax.plot(data['Lbol'],H07(data['Lbol'],'SX'),'--',dashes=(25,15), c='chocolate')
ax.plot([],[],'--',dashes=(25,15), c='k',label=r'$\rm Hopkins+$ $\rm 2007$')

prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=2,ncol=1,frameon=False)
ax.set_xlabel(r'$\log{(L_{\rm bol}[{\rm erg}/{\rm s}])}$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(L_{\rm bol}/L_{\rm band})}$',fontsize=40,labelpad=5)

ax.set_xlim(42.1,47.8)
ax.set_ylim(0.55,2.45)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("./bolcorr.pdf",fmt='pdf')
#plt.show()

