# McGreer et al. 2013
# Measuring M_1450
from data_copy import *
import numpy as np

def load_mcgreer13_lf_data(z): 
	if (z <= 4.7) or (z > 5.1): return False
	else:	
		filename = datapath+'kk18_compilation.dat'
		data = np.genfromtxt(filename,names=['counter','sample', 'z_bin', 'z_min', 'z_max', 'z_mean', 'M1450', 
			'left', 'right', 'log_phi', 'uperr', 'lowerr', 'nqso', 'Veff', 'P'])

		id = (data["sample"]==8) & (data["z_min"]< z) & (data["z_max"]>=z)

		M_1450  = data['M1450'][id]
		logphi  = data['log_phi'][id]
		dphi = (data['uperr'][id]+data['lowerr'][id])/2.
			 
		PHI_1450 = logphi
		DPHI_1450 = dphi
		return M_1450, PHI_1450, DPHI_1450
