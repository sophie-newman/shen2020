# Willott et al. 2010
# Measuring M_1450
from data import *
import numpy as np

def load_willott10_lf_data(z): 
	if (z <= 5.75) or (z > 6.45): return False
	else:	
		filename = datapath+'kk18_compilation.dat'
		data = np.genfromtxt(filename,names=['counter','sample', 'z_bin', 'z_min', 'z_max', 'z_mean', 'M1450', 
			'left', 'right', 'log_phi', 'uperr', 'lowerr', 'nqso', 'Veff', 'P'])

		id = (data["sample"]==10) & (data["z_min"]< z) & (data["z_max"]>=z)

		M_1450  = data['M1450'][id]
		logphi  = data['log_phi'][id]
		dphi = (data['uperr'][id]+data['lowerr'][id])/2.
			 
		PHI_1450 = logphi
		DPHI_1450 = dphi
		return M_1450, PHI_1450, DPHI_1450


