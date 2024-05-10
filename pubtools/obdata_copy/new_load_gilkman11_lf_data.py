# Gilkman et al. 2011
# Measuring M_1450, also reevaluate the M1450 for SDSS quasars
from data_copy import *
import numpy as np

def load_gilkman11_lf_data(z): 
	if (z <= 3.8) or (z > 5.2): return False
	else:	
		filename = datapath+'kk18_compilation.dat'
		data = np.genfromtxt(filename,names=['counter','sample', 'z_bin', 'z_min', 'z_max', 'z_mean', 'M1450', 
			'left', 'right', 'log_phi', 'uperr', 'lowerr', 'nqso', 'Veff', 'P'])

		id = (data["sample"]==6) & (data["z_min"]< z) & (data["z_max"]>=z)

		M_1450  = data['M1450'][id]
		logphi  = data['log_phi'][id]
		dphi = (data['uperr'][id]+data['lowerr'][id])/2.
			 
		PHI_1450 = logphi
		DPHI_1450 = dphi
		return M_1450, PHI_1450, DPHI_1450