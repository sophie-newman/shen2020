# Yang et al. 2016
# Measuring M_1450
from data import *
import numpy as np

def load_yang16_lf_data(z): 
	if (z < 4.7) or (z > 5.4): return False
	else:	
		z_c = 5.05
		M_1450=np.array([-28.99, -28.55, -28.05, -27.55, -27.05])
		PHI_1450=np.array([-9.48, -9.86, -9.36, -9.09, -8.70])
		sigma=np.array([0.33, 0.08, 0.15, 0.19, 0.32])*1e-9	

	M_1450 = M_1450 - 2.5*np.log10(lum_correct_cosmo_flexible(z_c, 0.7, 0.272))
	PHI_1450 = PHI_1450 + np.log10(phi_correct_cosmo_flexible(z_c, 0.7, 0.272))
	sigma = sigma * phi_correct_cosmo_flexible(z_c, 0.7, 0.272)

	DPHI_1450= ( (np.log10(10**PHI_1450+sigma)-PHI_1450)+(PHI_1450-np.log10(10**PHI_1450-sigma)) )/2.
	return M_1450, PHI_1450, DPHI_1450
