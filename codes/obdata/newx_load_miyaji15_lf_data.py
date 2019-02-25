# Miyaji et al. 2015
# measured 2-10 keV
from data import *
import numpy as np

def load_miyaji15_lf_data(z): 
	if (z <= 0.015) or (z > 4.60): return False
	else:	
		filename = datapath+'miyaji2015_table4.dat'
		data = np.genfromtxt(filename,names=['z_min', 'z_max', 'z_mean', 
			'L_min', 'L_max', 'L_c', 'nqso', 'phi', 'phierr'])

		id = (data["z_min"]< z) & (data["z_max"]>=z)

		L_X  = (data['L_min'][id]+data['L_max'][id])/2.
		phi  = data['phi'][id]
		phierr  = data['phierr'][id]
			 
		PHI_X = np.log10(phi)
		DPHI_X = ((np.log10(phi+phierr)-np.log10(phi)) + (np.log10(phi)-np.log10(phi-phierr)))/2.
		return L_X, PHI_X, DPHI_X
