from data import *
import numpy as np 
from lf_shape import *
import scipy.interpolate as inter
import lmfit
import matplotlib.pyplot as plt
# fit the luminosity function based on datasets at a given redshift
import sys

redshift=0.2

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

L_bol_grid_shrinked = np.linspace(10,14,50)
bolLF_x = L_bol_grid_shrinked + L_solar
bolLF_y = LF(L_bol_grid_shrinked,parameters_global_2)

######################################################

C = np.log10(1.26e38)
shift = C - 1.5 # the center of lambda grid

def kernel_func1(x, knee=np.log10(1.5), alpha=-0.6):
	#referencde point is lamb=0
	P = 10**((x-knee)*(alpha+1))*np.exp(-10**(x-knee))
	P[np.invert(np.isfinite(P))]=1e-30
	P[P<=1e-30]=1e-30
	return P

def kernel_func2(x,x0,sigma):
	P = np.exp(-0.5* (x-x0)**2/sigma**2 )
	P[np.invert(np.isfinite(P))]=1e-30
	P[P<=1e-30]=1e-30
	return P

x=np.linspace(-4.,1.,51)
kernel1 = kernel_func1(x, knee=np.log10(1.5), alpha=-0.6) 
kernel1 = kernel1/np.sum(kernel1)

kernel2 =kernel_func2(x, x0=-1.9+0.45*redshift, sigma=1.03-0.15*redshift)
kernel2 = kernel2/np.sum(kernel2)

kernel = 0.62*kernel1+0.38*kernel2
kernel = np.flip(kernel)
'''
plt.plot(np.flip(x),kernel)
plt.yscale('log')
plt.show()
exit()
'''

def BHMF1(logM, logphi_s, logM_s, alpha, beta):
	a = 10.**(logM-logM_s)
	return 10.**logphi_s * a**(1.+alpha) * np.exp(1-a)

def BHMF2(logM, logphi_s, logM_s, alpha, beta):
	a = 10.**(logM-logM_s)
	return 10.**logphi_s / (a**alpha + a**beta)

func_for_fit = BHMF2

def residual(pars):
	parvals  = pars.valuesdict()
	logphi_s = parvals['logphi_s']
	logM_s   = parvals['logM_s']
	alpha    = parvals['alpha']
	beta     = parvals['beta']

	x= np.linspace(6,12,1000)
	original = func_for_fit(x, logphi_s, logM_s, alpha, beta)
	convolved = np.convolve(original, kernel, 'valid')

	discard=int((len(kernel)-1)/2.)
	totnum_original = np.sum(original[discard:-discard])
	totnum_convolved= np.sum(convolved)
	lbol_mod = np.log10( convolved/totnum_convolved*totnum_original )

	res = (inter.interp1d(x[discard:-discard]+shift, lbol_mod)(bolLF_x) - bolLF_y)/0.3
	return res[np.isfinite(res)]

def plot_bestfit(logphi_s, logM_s, alpha, beta):
	x= np.linspace(6,12,1000)
	original = func_for_fit(x, logphi_s, logM_s, alpha, beta)
	convolved = np.convolve(original, kernel, 'valid')

	discard=int((len(kernel)-1)/2.)
	totnum_original = np.sum(original[discard:-discard])
	totnum_convolved= np.sum(convolved)
	lbol_mod = np.log10( convolved/totnum_convolved*totnum_original )

	plt.plot(x[discard:-discard]+shift, lbol_mod)
	plt.plot(bolLF_x, bolLF_y)
	plt.show()

params = lmfit.Parameters()
params.add_many(('logphi_s' , -3,  True, None, None, None, None),
                ('logM_s' ,   9.,  True, None, None, None, None),
                ('alpha',     0.5, True, None, None, None, None),
                ('beta',      1.,  True, None, None, None, None))

fitter = lmfit.Minimizer(residual, params, scale_covar=True,nan_policy='raise',calc_covar=True)
result=fitter.minimize(method='leastsq')
print "bestfit:"
result.params.pretty_print()

#print chisq(result.params['logphi_s'].value, result.params['logM_s'].value, result.params['alpha'].value ,result.params['beta'].value)
plot_bestfit(result.params['logphi_s'].value, result.params['logM_s'].value, result.params['alpha'].value, result.params['beta'].value)




