# 
# Barger al. 2005 data
#
from data import *
import numpy as np

def load_barger_lf_data(z):  # L_HX, PHI_HX, DPHI_HX, z
	# determine which redshift interval its in
	if (z <= 0.1): WHICH_BLOCK = 0
	if (z > 0.1) and (z <= 0.4): WHICH_BLOCK = 1
	if (z > 0.4) and (z <= 0.8): WHICH_BLOCK = 2
	if (z > 0.8) and (z <= 1.2): WHICH_BLOCK = 3
	if (z > 1.2) and (z <= 1.5): WHICH_BLOCK = 6
	if (z > 1.5) and (z <= 3.0): WHICH_BLOCK = 4
	if (z > 3.0) and (z <= 5.0): WHICH_BLOCK = 5
	if (z > 5.0) and (z <= 6.5): WHICH_BLOCK = 6
	if (z >  6.5): WHICH_BLOCK = 7
	
	if (WHICH_BLOCK == 0) or (WHICH_BLOCK == 7): return False
	else:
		if (WHICH_BLOCK == 1):
			L_HX = np.array([42.25, 42.75, 43.25, 43.75])
			P_HX = np.array([-4.15, -4.10, -4.52, -5.16])
			D_HX = np.array([ 0.25,  0.20,  0.20,  0.30])
		if (WHICH_BLOCK == 2):
			L_HX = np.array([42.25, 42.75, 43.25, 43.75, 44.25, 44.75])
			P_HX = np.array([-3.70, -4.15, -4.30, -4.55, -5.30, -6.40])
			D_HX = np.array([ 0.20,  0.25,  0.35,  0.20,  0.25,  0.40])
		if (WHICH_BLOCK == 3):
			L_HX = np.array([42.25, 42.75, 43.25, 43.75, 44.25, 44.75, 45.25])
			P_HX = np.array([-3.52, -4.05, -4.15, -4.40, -4.74, -5.89, -6.66])
			D_HX = np.array([ 0.20,  0.25,  0.25,  0.25,  0.35,  0.35,  0.40])
		if (WHICH_BLOCK == 4):
			L_HX = np.array([42.25, 42.75, 43.25, 43.75, 44.25, 44.75, 45.25])
			P_HX = np.array([-4.70, -4.75, -5.00, -5.10, -5.22, -5.58, -6.30])
			D_HX = np.array([ 0.80,  0.70,  0.80,  0.25,  0.20,  0.35,  0.45])
		if (WHICH_BLOCK == 5):
			L_HX = np.array([42.75, 43.25, 43.75, 44.25, 44.75, 45.25])
			P_HX = np.array([-5.30, -5.52, -5.30, -6.05, -6.52, -6.75])
			D_HX = np.array([ 0.70,  0.45,  0.40,  0.30,  0.45,  0.60])
		
		# Barger et al. 2003 (highest-redshift) data
		if (WHICH_BLOCK == 6):
			L_HX = np.array([43.5])
			P_HX = np.array([-5.95])
			D_HX = np.array([ 0.60])
		
		L_HX = (L_HX - L_solar) + np.log10(1.25) # convert the Barger 2-8kev to 2-10kev (gamma=1.4)
			# (note if instead were to use a gamma=1.8 photon index would get 1.19, so very small difference
		PHI_HX  = P_HX
		DPHI_HX = D_HX
		return L_HX, P_HX, DPHI_HX

