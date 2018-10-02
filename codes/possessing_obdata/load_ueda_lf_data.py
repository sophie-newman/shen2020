# 
# Ueda et al. 2003 data
#
from data import *
import numpy as np

def load_ueda_lf_data(z): # L_HX, PHI_HX, DPHI_HX, z
	# determine which redshift interval its in
	if ((z >= 0.0) and (z < 0.2)): WHICH_BLOCK = 1
	if ((z >= 0.2) and (z < 0.4)): WHICH_BLOCK = 2
	if ((z >= 0.4) and (z < 0.8)): WHICH_BLOCK = 3
	if ((z >= 0.8) and (z < 1.6)): WHICH_BLOCK = 4
	if ((z >= 1.6) and (z <=3.0)): WHICH_BLOCK = 5
	if (z > 3.0): WHICH_BLOCK = 6
	
	if (WHICH_BLOCK == 6): return False
	else:
		BEGIN = False
		filename=datapath+'load_ueda_lf_data.dat'
		with open(filename, 'r') as f:
			i=0
			j=0
			for line in f.readlines():
				elements=line.split()
				if (BEGIN == True) and (j<N_LUM_BINS):
					L_HX[j] = float(elements[2])
					PHI_HX[j] = float(elements[3])
					DPHI_HX[j] = float(elements[4])
					j+=1
		 		if (len(elements[0]) == 1) and ( int(elements[0])==WHICH_BLOCK ):
					BEGIN = True
					N_LUM_BINS = int(elements[1])
					L_HX = np.zeros(N_LUM_BINS)
					PHI_HX  = np.zeros(N_LUM_BINS)
					DPHI_HX = np.zeros(N_LUM_BINS)
				i+=1

		DPHIm_HX = -1.0*DPHI_HX
		DPHIp_HX = DPHI_HX
		L_HX = (L_HX - L_solar)
		return L_HX, PHI_HX, DPHI_HX
