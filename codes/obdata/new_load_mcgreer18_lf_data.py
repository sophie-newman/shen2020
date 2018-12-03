# McGreer et al. 2018
# Measuring M_1450
from data import *
import numpy as np

def load_mcgreer18_lf_data(z): 
	if (z < 4.7) or (z > 5.4): return False
	else:	
		z_mean=4.9
		M_1450=np.array([-26.35, -25.25, -24.35, -23.65, -22.90])
		PHI_1450=np.array([-8.12, -7.56, -7.25, -7.32, -7.32])
		sigma=np.array([4.34, 12.70, 18.05, 23.77, 28.24])*1e-9

		M_1450 = M_1450 - 2.5*np.log10(lum_correct_cosmo_flexible(z_mean, 0.7, 0.272))
		PHI_1450 = PHI_1450 + np.log10(phi_correct_cosmo_flexible(z_mean, 0.7, 0.272))
		sigma = sigma * phi_correct_cosmo_flexible(z_mean, 0.7, 0.272)

		DPHI_1450= ( (np.log10(10**PHI_1450+sigma)-PHI_1450)+(PHI_1450-np.log10(10**PHI_1450-sigma)) )/2.

		return M_1450, PHI_1450, DPHI_1450
#not completed