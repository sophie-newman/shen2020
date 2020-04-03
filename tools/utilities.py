### a python file that includes all the useful functions/models of this project
# pieces of cdoes from other places of the repository are copied here

from data import *
from astropy.cosmology import FlatLambdaCDM
import astropy.constants as con
import numpy as np
from scipy.optimize import fsolve

import numpy as np 
from lf_shape import *
import scipy.interpolate as inter
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import least_squares
from ctypes import *
import ctypes

##Section 1: the bolometric correction
#codes to load the bolometric corrections and their dispersions as a function of bolometric luminosity

def bolometric_correction(L_bol,nu):  #L_bol is in log10 and is in unit of solar luminosity  
	x = L_bol - 10.
	if nu==0.: return L_bol
	elif nu < 0:
		if (nu== -1.): P0, P1, P2, P3 =3.75876730, -0.36057087, 9.82954354, -6.29078046e-3   #B band
		if (nu== -2.): P0, P1, P2, P3 =4.36089518, -0.36057086, 11.4041656, -6.29077790e-3   #IR 15micron
		if (nu== -3.): P0, P1, P2, P3 =5.71209149, -0.02553868, 17.6679057,     0.27750559   #soft Xray
		if (nu== -4.): P0, P1, P2, P3 =4.07349835, -0.02553868, 12.5996205,     0.27750559   #hard Xray
		if (nu== -5.): P0, P1, P2, P3 =4.86961646, -0.00629077, 1.86211735,    -0.36057082   #FUV
		lband = P0 * np.power(10.,P1*x) + P2*np.power(10.,P3*x)
		return L_bol - np.log10(lband) 
	#if not one of the specified bands, then take advantage of the fact that our model 
	#   spectrum is not l-dependent below 500 angstroms or above 50 angstroms, so 
	#   just take the appropriate ratios to renormalize to those luminosities
	nu_angstrom = con.c.value/1e-10
	nu500 = nu_angstrom/500.
	nu50 = nu_angstrom/50.
	if (nu <= nu500): #just take the ratio relative to B-band
		P = (3.75876730,-0.36057087,9.82954354,-6.29078046e-3)
		lband = P[0]*np.power(10.,P[1]*x) + P[2]*np.power(10.,P[3]*x)
		return np.log10(return_ratio_to_b_band(nu)) + L_bol - np.log10(lband)
	elif (nu >= nu50): #just take the ratio relative to the hard X-rays
		P = (4.07349835,-0.02553868,12.5996205,0.27750559)
		lband = P[0]*np.power(10.,P[1]*x) + P[2]*np.power(10.,P[3]*x)
		return np.log10(return_ratio_to_hard_xray(nu)) + L_bol - np.log10(lband)
	elif (nu > nu500) and (nu < nu50): #interpolate between both regimes
		P = (3.75876730,-0.36057087,9.82954354,-6.29078046e-3)
		L500 = return_ratio_to_b_band(nu500)/(P[0]*np.power(10.,P[1]*x)+P[2]*np.power(10.,P[3]*x))
		Q = (4.07349835,-0.02553868,12.5996205,0.27750559)
		L50  = return_ratio_to_hard_xray(nu50)/(Q[0]*np.power(10.,Q[1]*x)+Q[2]*np.power(10.,Q[3]*x))
		L00  = np.log10(L500) + np.log10(L50/L500) * (np.log10(nu/nu500)/np.log10(nu50/nu500))
		return L00 + L_bol

# inverse of the bolometric correction, used to move data points to the bolometric plane
def bolometric_correction_inverse(L_obs,nu):    #L_obs is in log10 and is in unit of solar luminosity
	def fobjective(L_bol):
		return bolometric_correction(L_bol,nu)-L_obs
	return fsolve(fobjective,x0=L_obs)[0]

#return the lognormal dispersion in bolometric corrections for a given band and luminosity
def fit_func_disp(x,A,B,x0,sig):
	return B+A*norm.cdf(x,loc=x0,scale=sig)

def band_dispersion(L_bol,nu): #L_bol is in log10 and is in unit of solar luminosity 
	x = L_bol + L_solar

	sx_floor = 0.010 #minimum value of dispersion to enforce
	if nu == 0.:  return 0.01
	if nu < 0.:
		if nu==(-1.): P0, P1, P2, P3=-0.3826995, 0.4052673, 42.3866639, 2.3775969
		if nu==(-2.): P0, P1, P2, P3=-0.3380714, 0.4071626, 42.1588292, 2.1928345
        if nu==(-3.): P0, P1, P2, P3=0.07969176, 0.1803728, 44.1641156, 1.4964823
        if nu==(-4.): P0, P1, P2, P3=0.19262562, 0.0659231, 42.9876632, 1.8829639
        if nu==(-5.): P0, P1, P2, P3=-0.3719955, 0.4048693, 42.3073116, 2.3097825
        return fit_func_disp(x, P0, P1, P2, P3)
	# interpolate between the known ranges, and (conservatively) hold constant 
	#   outside of them. roughly consistent with Richards et al. 2006 dispersions, 
	#   but uncertainty in how large the dispersions should be yields ~10% uncertainties 
	#   between 15microns and 10keV, and larger outside those ranges (since the 
	#   dispersions there are poorly determined) -- still, the lognormal dispersions 
	#   vary relatively weakly over observed ranges, so these are probably the 
	#   maximal uncertainties due to this effect
	nu15 = 2.00e13
	nuBB = 6.81818e14
	nuSX = 0.5 * 2.418e17
	nuHX = 10. * 2.418e17
	if (nu < nu15): 
		P0, P1, P2, P3 = -0.3380714, 0.4071626, 42.1588292, 2.1928345 
		return fit_func_disp(x, P0, P1, P2, P3)
	if (nu >=nuHX): 
		P0, P1, P2, P3 = 0.19262562, 0.0659231, 42.9876632, 1.8829639
		return fit_func_disp(x, P0, P1, P2, P3)
	if (nu >= nu15) and (nu< nuBB): 
		P0, P1, P2, P3 = -0.3380714, 0.4071626, 42.1588292, 2.1928345
		sf1 = fit_func_disp(x, P0, P1, P2, P3)
		P0, P1, P2, P3 = -0.3826995, 0.4052673, 42.3866639, 2.3775969
		sf2 = fit_func_disp(x, P0, P1, P2, P3)
		sx = sf1 + (sf2-sf1) * (np.log10(nu/nu15)/np.log10(nuBB/nu15))
		if (sx<=sx_floor): sx=sx_floor
		return sx
	if (nu >= nuBB) and (nu< nuSX): 
		P0, P1, P2, P3 = -0.3826995, 0.4052673, 42.3866639, 2.3775969
		sf1 = fit_func_disp(x, P0, P1, P2, P3)
		P0, P1, P2, P3 = 0.07969176, 0.1803728, 44.1641156, 1.4964823
		sf2 = fit_func_disp(x, P0, P1, P2, P3)
		sx = sf1 + (sf2-sf1) * (np.log10(nu/nu15)/np.log10(nuBB/nu15))
		if (sx<=sx_floor): sx=sx_floor
		return sx
	if (nu >= nuSX) and (nu< nuHX): 
		P0, P1, P2, P3 = 0.07969176, 0.1803728, 44.1641156, 1.4964823
		sf1 = fit_func_disp(x, P0, P1, P2, P3)
		P0, P1, P2, P3 = 0.19262562, 0.0659231, 42.9876632, 1.8829639
		sf2 = fit_func_disp(x, P0, P1, P2, P3)
		sx = sf1 + (sf2-sf1) * (np.log10(nu/nu15)/np.log10(nuBB/nu15))
		if (sx<=sx_floor): sx=sx_floor
		return sx

# load the optical-IR template, based on the observations in text and specifically 
def return_ratio_to_b_band(nu):  #return the luminosity at a given frequency with respect to B band
	log_nu=np.array([    12.50, 12.52, 12.54, 12.56, 12.58, 12.60, 12.62, 12.64, 12.66, 12.68, 
	12.70, 12.72, 12.74, 12.76, 12.78, 12.80, 12.82, 12.84, 12.86, 12.88, 12.90, 12.92, 12.94, 
	12.96, 12.98, 13.00, 13.02, 13.04, 13.06, 13.08, 13.10, 13.12, 13.14, 13.16, 13.18, 13.20, 
	13.22, 13.24, 13.26, 13.28, 13.30, 13.32, 13.34, 13.36, 13.38, 13.40, 13.42, 13.44, 13.46, 
	13.48, 13.50, 13.52, 13.54, 13.56, 13.58, 13.60, 13.62, 13.64, 13.66, 13.68, 13.70, 13.72, 
	13.74, 13.76, 13.78, 13.80, 13.82, 13.84, 13.86, 13.88, 13.90, 13.92, 13.94, 13.96, 13.98, 
	14.00, 14.02, 14.04, 14.06, 14.08, 14.10, 14.12, 14.14, 14.16, 14.18, 14.20, 14.22, 14.24, 
	14.26, 14.28, 14.30, 14.32, 14.34, 14.36, 14.38, 14.40, 14.42, 14.44, 14.46, 14.48, 14.50, 
	14.52, 14.54, 14.56, 14.58, 14.60, 14.62, 14.64, 14.66, 14.68, 14.70, 14.72, 14.74, 14.76, 
	14.78, 14.80, 14.82, 14.84, 14.86, 14.88, 14.90, 14.92, 14.94, 14.96, 14.98, 15.00, 15.02, 
	15.04, 15.06, 15.08, 15.10, 15.12, 15.14, 15.16, 15.18, 15.20, 15.22, 15.24, 15.26, 15.28, 
	15.30, 15.32, 15.34, 15.36, 15.38, 15.40, 15.42, 15.44, 15.46, 15.48, 15.50, 15.52, 15.54, 
	15.56, 15.58, 15.60, 15.62, 15.64, 15.66, 15.68])
    
	log_nuLnu=np.array([      43.77927983, 43.82927983, 43.88927983, 43.93927983, 
	43.98927983, 44.03927983, 44.08927983, 44.12927983, 44.16927983, 44.20927983, 
	44.24927983, 44.27927983, 44.30927983, 44.33927983, 44.35927983, 44.38927983, 
	44.40927983, 44.42927983, 44.44927983, 44.46927983, 44.48927983, 44.50927983, 
	44.52927983, 44.54927983, 44.55927983, 44.57927983, 44.59927983, 44.60927983, 
	44.61927983, 44.62927983, 44.64927983, 44.65927983, 44.66927983, 44.67927983, 
	44.67927983, 44.68927983, 44.69927983, 44.70927983, 44.71927983, 44.71927983, 
	44.72827983, 44.73727983, 44.73727983, 44.74527983, 44.74527983, 44.75227983, 
	44.75227983, 44.75827983, 44.75727983, 44.76227983, 44.76127983, 44.76527983, 
	44.76327983, 44.76127983, 44.75927983, 44.76027983, 44.75727983, 44.75527983, 
	44.75427983, 44.75227983, 44.75127983, 44.75027983, 44.74927983, 44.74927983, 
	44.75027983, 44.75127983, 44.75327983, 44.75627983, 44.75827983, 44.76127983, 
	44.76427983, 44.76727983, 44.76927983, 44.77127983, 44.77327983, 44.77327983, 
	44.77227983, 44.77027983, 44.76627983, 44.76127983, 44.75527983, 44.74727983, 
	44.73827983, 44.72627983, 44.71427983, 44.70027983, 44.68627983, 44.67127983, 
	44.65727983, 44.64427983, 44.63227983, 44.62327983, 44.61527983, 44.61027983, 
	44.60727983, 44.60627983, 44.60727983, 44.61027983, 44.61427983, 44.61927983, 
	44.62627983, 44.63327983, 44.64227983, 44.65127983, 44.66227983, 44.67427983, 
	44.68727983, 44.70027983, 44.71127983, 44.71927983, 44.72527983, 44.73027983, 
	44.73727983, 44.74727983, 44.75927983, 44.77227983, 44.78627983, 44.79927983, 
	44.81227983, 44.82827983, 44.84727983, 44.86827983, 44.89127983, 44.91427983, 
	44.93627983, 44.95727983, 44.97427983, 44.98727983, 44.99527983, 45.00027983, 
	45.00227983, 45.00427983, 45.00927983, 45.01627983, 45.02627983, 45.03827983, 
	45.05027983, 45.06227983, 45.07327983, 45.08327983, 45.09227983, 45.10027983, 
	45.10627983, 45.10927983, 45.11927983, 45.12327983, 45.12527983, 45.12627983, 
	45.12927983, 45.13727983, 45.14727983, 45.13327983, 45.11927983, 45.10527983, 
	45.09127983, 45.07727983, 45.06327983, 45.04927983, 45.03527983, 45.02127983])

	#want the ratio with respect to the intrinsic B-band:
	nu_BB = 14.833657
	L_BB  = 44.793331
	log_nu_obs = np.log10(nu)
	nuLnu_obs = 0.

	if (log_nu_obs < log_nu[0]):   nuLnu_obs = log_nuLnu[0]
	if (log_nu_obs > log_nu[159]): nuLnu_obs = log_nuLnu[159]
	if ((log_nu_obs>=log_nu[0]) and (log_nu_obs<=log_nu[159])): 
		n0 = int((log_nu_obs-log_nu[0])/0.02)
		nuLnu_obs = log_nuLnu[n0] + (log_nuLnu[n0+1]-log_nuLnu[n0]) * ((log_nu_obs-log_nu[n0])/(log_nu[n0+1]-log_nu[n0]))

	return np.power(10.0,nuLnu_obs-L_BB)

##Section 2: Quasar luminosity functions
# codes to get the best-fit bolometric QLFs and to get predicted QLFs in different bands

def return_bolometric_qlf(redshift, model='A'):
	if model=='A':
		#load the global fit A
		source = np.genfromtxt("../plots/Fit_parameters/codes/zevolution_fit_global.dat",names=True)
		zref = 2.
		p=source['value'][ source['paraid']==0 ]
		gamma1 = polynomial(redshift,p,2)
		p=source['value'][ source['paraid']==1 ]
		gamma2 = doublepower(redshift,(p[0],zref, p[1], p[2]))
		p=source['value'][ source['paraid']==2 ]
		logphi = polynomial(redshift,p,1) 
		p=source['value'][ source['paraid']==3 ]
		Lbreak = doublepower(redshift,(p[0],zref, p[1], p[2]))
		parameters_global = np.array([gamma1,gamma2,logphi,Lbreak])
	elif model=='B':
		#load the global fit B with shallow faint end slopes at high z
		source = np.genfromtxt("../plots/Fit_parameters/codes/zevolution_fit_global_shallowfaint.dat",names=True)
		zref = 2.
		p=source['value'][ source['paraid']==0 ]
		gamma1 = powerlaw_gamma1(redshift,(p[0],zref,p[1]))
		p=source['value'][ source['paraid']==1 ]
		gamma2 = doublepower(redshift,(p[0],zref,p[1],p[2]))
		p=source['value'][ source['paraid']==2 ]
		logphi = polynomial(redshift,p,1)
		p=source['value'][ source['paraid']==3 ]
		Lbreak = doublepower(redshift,(p[0],zref,p[1],p[2]))
		parameters_global = np.array([gamma1,gamma2,logphi,Lbreak])
	else: 
		print "model not found"
		exit()
	
	x = L_bol_grid + L_solar     #luminosities in log10
	y = LF(L_bol_grid,parameters_global)   #number densities in log10
	return x, y

def return_qlf_in_band(redshift, nu, model='A'):  #nu here has the same definitions as in Section 1
	if model=='A':
		#load the global fit A
		source = np.genfromtxt("../plots/Fit_parameters/codes/zevolution_fit_global.dat",names=True)
		zref = 2.
		p=source['value'][ source['paraid']==0 ]
		gamma1 = polynomial(redshift,p,2)
		p=source['value'][ source['paraid']==1 ]
		gamma2 = doublepower(redshift,(p[0],zref, p[1], p[2]))
		p=source['value'][ source['paraid']==2 ]
		logphi = polynomial(redshift,p,1) 
		p=source['value'][ source['paraid']==3 ]
		Lbreak = doublepower(redshift,(p[0],zref, p[1], p[2]))
		parameters_global = np.array([gamma1,gamma2,logphi,Lbreak])
	elif model=='B':
		#load the global fit B with shallow faint end slopes at high z
		source = np.genfromtxt("../plots/Fit_parameters/codes/zevolution_fit_global_shallowfaint.dat",names=True)
		zref = 2.
		p=source['value'][ source['paraid']==0 ]
		gamma1 = powerlaw_gamma1(redshift,(p[0],zref,p[1]))
		p=source['value'][ source['paraid']==1 ]
		gamma2 = doublepower(redshift,(p[0],zref,p[1],p[2]))
		p=source['value'][ source['paraid']==2 ]
		logphi = polynomial(redshift,p,1)
		p=source['value'][ source['paraid']==3 ]
		Lbreak = doublepower(redshift,(p[0],zref,p[1],p[2]))
		parameters_global = np.array([gamma1,gamma2,logphi,Lbreak])
	else: 
		print "model not found"
		exit()

	# compile the c code before using it!!!
	c_extenstion = CDLL("./clib/convolve.so")
	convolve_c = c_extenstion.convolve
	convolve_c.restype = ctypes.POINTER(ctypes.c_double * N_bol_grid)

	if nu==-5:
		L_1450 = bolometric_correction(L_bol_grid,-5)
		nu_c = c_double(-5)
		redshift_c = c_double(redshift)
		dtg_c = c_double(return_dtg(redshift))
		input_c= np.power(10.,LF(L_bol_grid,parameters_global)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
		res = convolve_c(input_c,nu_c,redshift_c,dtg_c)
		res = [i for i in res.contents]
		PHI_1450 = np.array(res,dtype=np.float64)  
		x = -2.5*( L_1450 + L_solar - np.log10(Fab*(con.c.value/1450e-10)) )  # convert lum to mag
		xm = x.copy()
		ym = np.log10(PHI_1450) - np.log10(2.5) # convert to number density per mag
		return xm, ym
	else:
		L_band = bolometric_correction(L_bol_grid,nu)
		nu_c = c_double(nu)
		redshift_c = c_double(redshift)
		dtg_c = c_double(return_dtg(redshift))
		input_c= np.power(10.,LF(L_bol_grid,parameters_global)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
		res = convolve_c(input_c,nu_c,redshift_c,dtg_c)
		res = [i for i in res.contents]
		PHI_HX = np.array(res,dtype=np.float64)
		x = L_HX + L_solar
		y = np.log10(PHI_HX)
		return x,y

print return_qlf_in_band(5, -5)
