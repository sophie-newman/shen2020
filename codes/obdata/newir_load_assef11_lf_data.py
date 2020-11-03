# 
#  Assef et al. 2011
#
from data import *
import numpy as np 
from scipy.interpolate import interp1d
from scipy.interpolate import lagrange

def load_assef11_lf_data(z): # M_J, PHI_J, DPHI_J
	if (z <= 0) or (z > 5.85): return False
	else:	
		filename = datapath+'Assef2011.dat'
		data = np.genfromtxt(filename,names=True)

		id = (data["z_min"]< z) & (data["z_max"]>=z) & (data['label']==1) & (data['lf']!=data['lf_sigma'])

		Vega_to_AB = 0.912
		M_J  = data['MJ'][id] + Vega_to_AB
		P_J  = np.log10(data['lf'][id]/hubble**3)
		P_J_up  = np.log10(data['lf'][id]/hubble**3 + data['lf_sigma'][id]/hubble**3)
		P_J_down= np.log10(data['lf'][id]/hubble**3 - data['lf_sigma'][id]/hubble**3)
		D_J = (P_J_up - P_J_down)/2.

		PHI_J = P_J
		DPHI_J= D_J

		z_c = data["z_mean"][id][0]
		M_J = M_J - 2.5*np.log10(lum_correct_cosmo_flexible(z_c, 0.73, 0.3))
		PHI_J = PHI_J + np.log10(phi_correct_cosmo_flexible(z_c, 0.73, 0.3))

		L_J = np.log10( (10**(-0.4*M_J) * Fab) * (con.c.value/1.25e-6)) - L_solar

		L_IR = L_J + 0.11928 # convert to 15 micron
		PHI_IR = PHI_J + np.log10(2.5)
		DPHI_IR = DPHI_J

		return L_IR, PHI_IR, DPHI_IR

def return_assef11_lf_fitted(Llist,z):
	x = Llist + L_solar
	x = x - 0.11928

	logphis = lagrange(np.array([0.25,0.5,1.,2.,4.]),np.array([-3.41,-3.73,-4.17,-4.65, -5.77]))(z)
	alpha = 3.35
	beta  = 0.37
	Ms    = lagrange(np.array([0.5,1.,2.,4.]),np.array([-23.51,-24.64,-26.10,-27.08]))(z)

	Vega_to_AB = 0.912
	Ms = Ms + 0.912
	logLs = np.log10( (10**(-0.4*Ms) * Fab) * (con.c.value/1.25e-6))

	PHI=10**logphis / ( 10**((x-logLs)*(alpha-1)) + 10**((x-logLs)*(beta-1)) )
	return np.log10(PHI)








