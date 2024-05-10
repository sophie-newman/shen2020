# Giallongo et al. 2019
# Measuring M_1450
from data_copy import *
import numpy as np

def load_giallongo19_lf_data(z): 
	if (z <= 4) or (z > 6.1): return False
	elif (z > 4) and (z <= 5):	
		M_1450  = np.array([-19.,-20.,-21.,-22.])
		phi     = np.array([14.54,11.47,5.08,1.31])*1e-6 
		up      = np.array([8.72 ,6.88 ,3.45,1.74])*1e-6
		lo      = np.array([5.81 ,4.59 ,2.21,0.87])*1e-6
	elif (z > 5) and (z <= 6.1):			 
		M_1450  = np.array([-19.,-20.,-21.,-22.])
                phi     = np.array([7.27,4.77,0.69,0.62])*1e-6
                up      = np.array([7.12,3.79,1.61,1.44])*1e-6
                lo      = np.array([4.02,2.31,0.60,0.54])*1e-6

	PHI_1450 = np.log10(phi)
	DPHI_1450 = (np.log10(phi+up) - np.log10(phi-lo))/2.
	return M_1450, PHI_1450, DPHI_1450
