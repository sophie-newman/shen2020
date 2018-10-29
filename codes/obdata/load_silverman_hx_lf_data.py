from data import *
import numpy as np 

def load_silverman_hx_lf_data(z):  # L_HX, PHI_HX, DPHI_HX, z
	# determine which redshift interval its in
	WHICH_BLOCK = 8
	if ((z >= 0.2) and (z < 0.5)): WHICH_BLOCK = 1
	if ((z >= 1.5) and (z < 2.0)): WHICH_BLOCK = 4
	if ((z >= 3.0) and (z < 4.0)): WHICH_BLOCK = 6
	
	if (WHICH_BLOCK == 8): return False
	else:
		if (WHICH_BLOCK == 1):
			L_HX = np.array([42.25,42.75,43.25,43.75,44.25,44.75])
			P_HX = np.array([2.3e-4,1.7e-4,4.3e-5,2.3e-5,7.0e-6,2.0e-7])
			D_HX = np.array([0.115,  0.081, 0.145, 0.157, 0.146, 0.222])
		if (WHICH_BLOCK == 2):
			L_HX = np.array([  42.5,  43.5,  44.5])
			P_HX = np.array([2.2e-4,6.0e-5,2.5e-6])
			D_HX = np.array([ 0.125, 0.125, 0.215])
		if (WHICH_BLOCK == 3):
			L_HX = np.array([  43.5,  44.5])
			P_HX = np.array([5.0e-5,7.5e-6])
			D_HX = np.array([ 0.130, 0.204])
		if (WHICH_BLOCK == 4):
			L_HX = np.array([43.25,43.75,44.25,44.75,45.25])
			P_HX = np.array([2.1e-5,4.0e-5,2.8e-5,7.0e-6,1.2e-6])
			D_HX = np.array([ 0.226, 0.138, 0.101, 0.109, 0.221])
		if (WHICH_BLOCK == 5):
			L_HX = np.array([  43.5,  44.5])
			P_HX = np.array([6.8e-6,2.8e-6])
			D_HX = np.array([ 0.167, 0.196])
		if (WHICH_BLOCK == 6):
			L_HX = np.array([43.75,44.25,44.75,45.25])
			P_HX = np.array([9.0e-6,3.0e-6,1.6e-6,8.0e-8])
			D_HX = np.array([ 0.222, 0.204, 0.155, 0.574])
		if (WHICH_BLOCK == 7):
			L_HX = np.array([  43.5,  0.5*(44.5+45.5)])	# the latter is soft X-ray selected
			P_HX = np.array([1.2e-5,2.0e-7])
			D_HX = np.array([ 0.368, 0.398])
		L_HX = (L_HX - L_solar) + np.log10(1.19) # convert the Barger 2-8kev to 2-10kev (gamma=1.8)
		PHI_HX  = np.log10(P_HX)
		DPHI_HX = D_HX + 0.1
		return L_HX, PHI_HX, DPHI_HX
