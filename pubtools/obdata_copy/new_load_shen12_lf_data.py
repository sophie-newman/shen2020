# Shen et al. 2013, BOSS
# Measuring M_i(z=2) (converted to M_1450 in this code)
from data_copy import *
import numpy as np

def load_shen12_lf_data(z): 
	if (z <= 0.3) or (z > 5.0): return False
	else:	
		filename = datapath+'shen2012_table3.dat'
		data = np.genfromtxt(filename,names=["z_min", "z_max", "z_mid", "magz2", "logPhi", "elogPhi"])

		id = (data["logPhi"]>-20) & (data["z_min"]< z) & (data["z_max"]>=z)

		M_1450  = data['magz2'][id] + 1.486 #the paper
		logphi  = data['logPhi'][id]
		dphi = 10**data["elogPhi"][id]
			 
		PHI_1450 = logphi
		DPHI_1450 = (np.log10(10**logphi + dphi) - np.log10(10**logphi - dphi) )/2.

		return M_1450, PHI_1450, DPHI_1450
