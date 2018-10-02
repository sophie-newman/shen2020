# 
# Beckmann et al. 2006 local data
#
from data import *
import numpy as np

def load_beckmann_lf_data(z): # L_HX, PHI_HX, DPHI_HX, z
	# determine which redshift interval its in
	if (z > 0.044): return False 
	else:
		#beckmann -- they measure 20-40 keV :: from our model spectrum, 
 		#L(20-40) = 1.15 * L(2-10)
		l = np.array([41.25, 41.75, 42.25, 42.75, 43.25, 43.75, 44.25, 44.75, 45.25]) - np.log10(1.1522883)
		p = np.array([-3.4,  -3.7,  -4.35, -4.7,   -5.25, -6.0,  -7.0, -8.2, -9.1])
		d = np.array([0.3,    0.18, 0.17, 0.14,     0.1,   0.11,  0.15, 0.25, 0.26])
		
		L_HX = (l - L_solar) 
		PHI_HX  = p
		DPHI_HX = d
		#still slightly low, but this sample is *really* local (z_bar = 0.022), so it's ok
		return L_HX, PHI_HX, DPHI_HX

