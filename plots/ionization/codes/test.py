from data import *
import numpy as np 
from lf_shape import *
from new_load_kk18_lf_shape import *
import scipy.interpolate as inter
from scipy.integrate import quad
from scipy.integrate import romberg
from convolve import *
from convolve_h07 import *
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import least_squares
import astropy.constants as con
# fit the luminosity function based on datasets at a given redshift
from ctypes import *
import ctypes
import sys

#dummy functions, no real meaning
def Gamma_old(epsilon,z):
	return epsilon

def Gamma(epsilon,z):
	return epsilon

fit_evolve=np.genfromtxt("../../Fit_parameters/codes/zevolution_fit_global.dat",names=True)
paraid, pglobal, pglobal_err = fit_evolve['paraid'], fit_evolve['value'], (fit_evolve['uperr']+fit_evolve['loerr'])/2.

fit_evolve_shallowfaint=np.genfromtxt("../../Fit_parameters/codes/zevolution_fit_global_shallowfaint.dat",names=True)
paraid_shallowfaint, pglobal_shallowfaint, pglobal_err_shallowfaint = fit_evolve_shallowfaint['paraid'], fit_evolve_shallowfaint['value'], (fit_evolve_shallowfaint['uperr']+fit_evolve_shallowfaint['loerr'])/2.

#load the shared object file
c_extenstion = CDLL(homepath+'codes/c_lib/convolve.so')
convolve_c = c_extenstion.convolve
convolve_c.restype = ctypes.POINTER(ctypes.c_double * N_bol_grid)

#old convolution code for H07
c_extenstion_old = CDLL(homepath+'codes/c_lib/convolve_old.so')
convolve_c_old = c_extenstion_old.convolve
convolve_c_old.restype = ctypes.POINTER(ctypes.c_double * N_bol_grid)

def get_model_lf_global(parameters,nu,redshift,magnitude=False,model="Fiducial"):
	if model=="Fiducial":
		zref = 2.
		p=parameters[paraid==0]
		gamma1 = polynomial(redshift,p,2) 
		p=parameters[paraid==1]
		gamma2 = doublepower(redshift,(p[0],zref,p[1],p[2]))
		p=parameters[paraid==2]
		logphi = polynomial(redshift,p,1) 
		p=parameters[paraid==3]	
		Lbreak = doublepower(redshift,(p[0],zref,p[1],p[2]))
		parameters_at_z = np.array([gamma1,gamma2,logphi,Lbreak])
        elif model=="Shallowfaint":
                zref = 2.
                p=parameters[paraid_shallowfaint==0]
                gamma1 = powerlaw_gamma1(redshift,(p[0],zref,p[1])) 
                p=parameters[paraid_shallowfaint==1]
                gamma2 = doublepower(redshift,(p[0],zref,p[1],p[2]))
                p=parameters[paraid_shallowfaint==2]
                logphi = polynomial(redshift,p,1)
                p=parameters[paraid_shallowfaint==3]
                Lbreak = doublepower(redshift,(p[0],zref,p[1],p[2]))	
		parameters_at_z = np.array([gamma1,gamma2,logphi,Lbreak])
	elif model=="2":
		zref = 2.
                p=parameters[paraid==0]
                gamma1 = polynomial(redshift,p,2) + 0.7
                p=parameters[paraid==1]
                gamma2 = doublepower(redshift,(p[0],zref,p[1],p[2]))
                p=parameters[paraid==2]
                logphi = polynomial(redshift,p,1)
                p=parameters[paraid==3]
                Lbreak = doublepower(redshift,(p[0],zref,p[1],p[2]))
                parameters_at_z = np.array([gamma1,gamma2,logphi,Lbreak])

	print model, ":", parameters_at_z
	return get_model_lf(parameters_at_z,nu,redshift,magnitude=magnitude)

def get_model_lf(parameters,nu,redshift,magnitude=False):
        L_band = bolometric_correction(L_bol_grid,nu)
        nu_c = c_double(nu)
	redshift_c = c_double(redshift)
	dtg_c = c_double(return_dtg(redshift))
        input_c= np.power(10.,LF(L_bol_grid,parameters)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        res = convolve_c(input_c,nu_c,redshift_c,dtg_c)
        res = [i for i in res.contents]
        PHI_band = np.array(res,dtype=np.float64)
        if magnitude==False:
                return L_band, np.log10(PHI_band)
        else:
		M_1450 = -2.5*(L_band+L_solar-np.log10(Fab*(con.c.value/1450e-10)))
                PHI_1450 = np.log10(PHI_band) - np.log10(2.5)
                return M_1450, PHI_1450

def cumulative_emissivity(L_band,Phi_band,L_limit_low,L_limit_up):
	Mband, Mlow, Mup = L_band, L_limit_low, L_limit_up
	logphi=inter.interp1d(Mband,Phi_band)
	def emis(x):
		fnu = np.power(10.,-0.4*x)*3631*1e-23*4*np.pi*(10*con.pc.value*100)**2
		fnu912 = fnu * (912./1450.)**(0.61)
		return np.power(10.,logphi(x))*fnu912
	return romberg(emis,Mlow,Mup,divmax=20)
	#return quad(emis,Mlow,Mup)[0]

def kk18(M_list,z):
        c0=np.array([-7.084 , 0.753, -0.096])
        c1=np.array([-15.423, -6.725, 0.737, -0.029])
        c2=np.array([-2.973 , -0.347])
        c3=np.array([-2.545 , 1.581, 2.102, 1.965, -0.641])
	phi_s=10**( c0[0]*T0(1.+z)+c0[1]*T1(1.+z)+c0[2]*T2(1.+z))
	M_s  =c1[0]*T0(1.+z)+c1[1]*T1(1.+z)+c1[2]*T2(1.+z)+c1[3]*T3(1+z)
	alpha=c2[0]*T0(1.+z)+c2[1]*T1(1.+z)

        ksi = np.log10((1.+z)/(1.+c3[2]))
        beta = c3[0] + c3[1]/(10**(c3[3]*ksi)+10**(c3[4]*ksi))

        PHI=phi_s/( 10**(0.4*(alpha+1)*(M_list-M_s)) + 10**(0.4*(beta+1)*(M_list-M_s)) )
	print M_s, phi_s, alpha, beta
        return np.log10(PHI)

lowlimit=-35
uplimit =-13
redshift = 6

M_1450, PHI_1450 = get_model_lf_global(pglobal, -5, redshift, magnitude=True)
result1 = Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, uplimit), redshift)

M_1450, PHI_1450 = get_model_lf_global(pglobal_shallowfaint, -5, redshift, magnitude=True, model="Shallowfaint")
result2 = Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, uplimit), redshift)

PHI_1450 = kk18(M_1450, redshift)
result3 = Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, uplimit), redshift)

M_1450, PHI_1450 = get_model_lf_global(pglobal, -5, redshift, magnitude=True, model="2")
result4 = Gamma(cumulative_emissivity(M_1450, PHI_1450, lowlimit, uplimit), redshift)

print np.log10(result4)
print np.log10(result1)
print np.log10(result2)


