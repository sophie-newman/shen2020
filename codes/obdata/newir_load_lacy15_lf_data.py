# 
#  Lacy et al. 2011
#
from data import *
import numpy as np 

def load_lacy15_lf_data(z): # 5 micron
	if (z <= 0.05) or (z > 3.8): return False
	else:	
		filename = datapath+'Lacy2015.dat'
		data = np.genfromtxt(filename,names=True)

		id = (data["z_min"]< z) & (data["z_max"]>=z) 

		L_IR = data['nuLnu'][id] - L_solar
		PHI_IR  = data['logphi'][id]
		DPHI_IR = data['dlogphi'][id]

		L_IR = L_IR - 0.02139

		return L_IR, PHI_IR, DPHI_IR

def return_lacy15_lf_fitted(Llist,z):
	x = Llist + L_solar 

	logphis = -4.75
	gamma1 = 1.07
	gamma2 = 2.48
	k1, k2, k3 = 1.05, -4.71, 0.034
	logL0 = 31.92 + np.log10(con.c.value/5e-6)

	eps = np.log10((1.+z)/(1+2.5))
	logLs = logL0 + k1*eps + k2*eps**2 + k3*eps**3

	PHI=10**logphis / ( 10**((x-logLs)*gamma1) + 10**((x-logLs)*gamma2) )
	return np.log10(PHI)
