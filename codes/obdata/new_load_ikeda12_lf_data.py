# Ikeda et al. 2012
# Measuring M_1450
from data import *
import numpy as np

def load_ikeda12_lf_data(z): 
	if (z < 4.57) or (z > 5.57): return False
	else:	
		z_mean = 5.07
		M_1450= (-23.52-22.52)/2.
		PHI_1450= 0.87e-7
		sigma_up= 2.01e-7
		sigma_down= 0.72e-7

		PHI_1450 = np.log10(PHI_1450)
		DPHI_1450= ( (np.log10(10**PHI_1450+sigma_up)-PHI_1450)+(PHI_1450-np.log10(10**PHI_1450-sigma_down)) )/2.

	return M_1450, PHI_1450, DPHI_1450