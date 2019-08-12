# Giallongo et al. 2015
from data import *
import numpy as np

def return_giallongo15_lf_fitted(M_list,z):
	if (z>5.) and (z<6.5):
		phi_s = 10**(-5.8)
		M_s = -23.4
		alpha = 1.66
		beta = 3.35
	else: return False
	
	PHI=phi_s/( 10**(-0.4*(alpha-1)*(M_list-M_s)) + 10**(-0.4*(beta-1)*(M_list-M_s)) )
	return np.log10(PHI)

