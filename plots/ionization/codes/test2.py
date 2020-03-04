import numpy as np 
import scipy.interpolate as inter
from scipy.integrate import quad
from scipy.integrate import romberg
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import least_squares
import astropy.constants as con
# fit the luminosity function based on datasets at a given redshift
import scipy.special as spe
import mpmath
def cumulative(M_lim, Phis, Ms, alpha):
	return 10**(-0.4*Ms) * Phis*mpmath.gammainc(alpha+2., a=np.power(10.,-0.4*(M_lim-Ms)))

def cumulative_emissivity(L_band,Phi_band,L_limit_low,L_limit_up):
	Mband, Mlow, Mup = L_band, L_limit_low, L_limit_up
	logphi=inter.interp1d(Mband,Phi_band)
	def emis(x):
		fnu = np.power(10.,-0.4*x)#*3631*1e-23*4*np.pi*(10*con.pc.value*100)**2
		fnu912 = fnu #* (912./1450.)**(0.61)
		return np.power(10.,logphi(x))*fnu912
	return romberg(emis,Mlow,Mup,divmax=20)
	#return quad(emis,Mlow,Mup)[0]

def single_sch(M,phi_s,M_s,alpha):
        return np.log10( 0.4*np.log(10.)*phi_s*np.power(10.,-0.4*(M-M_s)*(alpha+1.))*np.exp(-np.power(10.,-0.4*(M-M_s))) )

M = np.linspace(-30, -16, 101)
Phi = single_sch(M, 1e-6, -24, -1.5)

print Phi

print cumulative(-18, 1e-6, -24, -1.5)
print cumulative_emissivity(M, Phi, -30, -18)






