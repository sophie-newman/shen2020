# 
# Shinozaki et al. 2006 local data
#
from data_copy import *
import numpy as np 

def load_shinozaki_lf_data(z): # L_HX, PHI_HX, DPHI_HX, z
	# determine which redshift interval its in
	if (z > 0.20): return False
	else:
		l = np.array([ 42.5, 43.5, 44.25, 44.85, 45.5 - 0.25, 46.0 - 0.25])
		p = np.log10( np.array([1.2e-4, 7.8e-6, 7.0e-7, 1.8e-8, 2.0e-9, 2.9e-10]) )
		d = np.array([0.13, 0.13, 0.10, 0.30, 0.42, 0.42])

		L_HX = (l - L_solar) 
		PHI_HX  = p
		DPHI_HX = d
		return L_HX, PHI_HX, DPHI_HX
