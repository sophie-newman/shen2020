# 
# Nandra et al. 2005 (LBG-selected) data
#
from data_copy import *
import numpy as np 

def load_nandra_lf_data(z): # L_HX, PHI_HX, DPHI_HX, z
	# determine which redshift interval its in
	WHICH_BLOCK = 2
	if ((z > 2.75) and (z <= 3.25)): WHICH_BLOCK = 1
	
	if (WHICH_BLOCK == 2): return False
	else:
		L_HX = np.array([0.5*(43.0+44.5)])
		P_HX = np.log10(np.array([4.2e-5])/1.5)
		D_HX = np.array([0.15])
		L_HX = (L_HX - L_solar) 
		PHI_HX  = P_HX
		DPHI_HX = D_HX
		return L_HX, PHI_HX, DPHI_HX
