# 
# La Franca et al. 2005 data
#
from data import *
import numpy as np

def load_lafranca_lf_data(z): # L_HX, PHI_HX, DPHI_HX, z
	# determine which redshift interval its in
	if ((z >= 0.0) and (z < 0.5)): WHICH_BLOCK = 1
	if ((z >= 0.5) and (z < 1.0)): WHICH_BLOCK = 2
	if ((z >= 1.0) and (z < 1.5)): WHICH_BLOCK = 3
	if ((z >= 1.5) and (z < 2.5)): WHICH_BLOCK = 4
	if ((z >= 2.5) and (z <=3.5)): WHICH_BLOCK = 5
	if (z > 3.5): WHICH_BLOCK = 6
	
	if (WHICH_BLOCK == 6): return False		
	else:
		if (WHICH_BLOCK == 1):
			L_HX = np.array([42.5, 43.5, 44.5, 45.5])
			P_HX = np.array([-3.7, -4.6, -6.2, -9.05])
			D_HX = np.array([ 0.2,  0.2,  0.2,  0.8])
		if (WHICH_BLOCK == 2):
			L_HX = np.array([42.5, 43.5, 44.5, 45.5])
			P_HX = np.array([-3.25, -4., -5.45, -8.0])
			D_HX = np.array([ 0.2,  0.2,  0.2,  0.33])
		if (WHICH_BLOCK == 3):
			L_HX = np.array([42.5, 43.5, 44.5, 45.5])
			P_HX = np.array([-3.1, -4.05, -5., -7.0])
			D_HX = np.array([ 0.3,  0.2,  0.2,  0.3])
		if (WHICH_BLOCK == 4):
			L_HX = np.array([43.5, 44.5, 45.5, 46.5])
			P_HX = np.array([-4.1, -4.9, -6.9, -8.95])
			D_HX = np.array([ 0.2,  0.2,  0.3,  0.5])
		if (WHICH_BLOCK == 5):
			L_HX = np.array([44.25, 45.0])
			P_HX = np.log10(np.array([2.0 * 8.0e-6, 1.0 * 2.0e-6]))
			D_HX = np.array([ 0.2, 0.25])
		L_HX = (L_HX - L_solar)
		PHI_HX  = P_HX
		DPHI_HX = D_HX
		return L_HX, PHI_HX, DPHI_HX


