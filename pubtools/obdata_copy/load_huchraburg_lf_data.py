#
# Huchra & Burg local AGN luminosity function
#
from data_copy import *
import numpy as np

def load_huchraburg_lf_data(z): # L_BB, PHI_BB, DPHI_BB, z
	if (z > 0.2): return False
	else:
		M_BB   = np.array([-17.75,-18.25,-18.75,-19.25,-19.75,-20.25,-20.75,-21.25,-21.75,-22.25])
		log_phi= np.array([-3.83, -3.69, -4.08, -3.90, -4.32, -4.46, -4.65, -5.29, -6.45, -6.54 ])
		
		L_BB   = (0.4*(M_sun_Bband_AB-M_BB)) - np.log10(hubble**2)
		PHI_BB = log_phi + np.log10(2.5) + np.log10(hubble**3)
		DPHI_BB= 0.0*log_phi + 0.2
		return L_BB, PHI_BB, DPHI_BB
