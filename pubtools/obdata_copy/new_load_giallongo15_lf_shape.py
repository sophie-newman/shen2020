# Giallongo et al. 2015
from data_copy import *
import numpy as np

def return_giallongo15_lf_fitted(M_list,z):
	if (z>=5.) and (z<6.5):
		phi_s = 10**(-5.8)
		M_s = -23.4
		alpha = 1.66
		beta = 3.35
	elif (z>=4.5) and (z<5.): 
		phi_s = 10**(-5.7)
                M_s = -23.6
                alpha = 1.81
                beta = 3.14
	elif (z>=4.) and (z<4.5):
		phi_s = 10**(-5.2)
                M_s = -23.2
                alpha = 1.52
                beta = 3.13
	else: return False
	
	PHI=phi_s/( 10**(-0.4*(alpha-1)*(M_list-M_s)) + 10**(-0.4*(beta-1)*(M_list-M_s)) )
	return np.log10(PHI)


