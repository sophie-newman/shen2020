# Cristiani 2004 faint at z=3.5-5.2
#
from data_copy import *
import numpy as np

def load_cristiani_lf_data(z): # L_BB, PHI_BB, DPHI_BB, z
	if ((z <= 4.0) or (z > 5.2)): return False
	else:
		M_1450_o = np.array([-22.3])			# z = 4.0 - 5.2
		phi_o    = np.array([10**(-5.9)])
		phi_p    = np.array([10**(-5.55)])
		phi_m    = np.array([10**(-6.4)])
		M_B_o = M_1450_o - 0.83	
		# note this is a different convention for M_1450 than other places, 
		#   but is directly from the paper (& used therein to define M_1450)
		dphi  = 0.5*(np.log10(phi_p) - np.log10(phi_m))
		L_BB  = (0.4*(M_sun_Bband_AB-M_B_o))
		PHI_BB= np.log10(2.5*phi_o)
		DPHI_BB=dphi	
		return L_BB, PHI_BB, DPHI_BB
