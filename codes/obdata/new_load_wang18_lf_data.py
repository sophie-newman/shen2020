# Wang et al. 2018
# Measuring M_1450
from data import *
import numpy as np

def load_wang18_lf_data(z): 
	if (z < 6.45) or (z > 7.05): return False
	else:
		z_c = 6.7
		M_1450=np.array([25.831933, 26.440126, 27.183824])
		M_up =np.array([26.096639, 26.581933, 27.385504])
		M_low=np.array([25.567227, 26.298319, 26.988445])
		PHI_1450=np.array([-9.083333, -9.446667, -9.823333])
		PHI_up =np.array([-8.936667, -9.360000, -9.703333])
		PHI_low=np.array([-9.293333, -9.550000, -9.990000])

	DPHI_1450= (PHI_up - PHI_low)/2.
	return M_1450, PHI_1450, DPHI_1450