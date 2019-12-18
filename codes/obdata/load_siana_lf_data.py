#
# Siana et al. 2006 (IR + optical selection (SWIRE)) at z~3
#   Measuring M_1450 -- convert to B-band following DR3 corrections
#
from data import *
import numpy as np

def load_siana_lf_data(z): # L_BB, PHI_BB, DPHI_BB, z
	if ((z <= 2.9) or (z > 3.4)): return False
	else:
		M_1450 = np.array([-23.75, -24.25, -24.75, -25.25, -25.75, -26.25])
		PHI    = np.array([1.0e-6, 8.0e-7, 5.5e-7, 3.0e-7, 1.2e-7, 1.8e-7]) * 2.5
		DPHI   = np.array([ 0.097,  0.056,  0.105,  0.125,  0.176,  0.143])	# way too optimistic
		DPHI   = np.array([ 0.097,  0.125,  0.155,  0.176,  0.269,  0.234])	# use lower errors - much better log-errors
		M_BB   = M_1450 - 0.776
		L_BB   = 0.4*(M_sun_Bband_AB-M_BB)
		PHI_BB = np.log10(PHI)
		DPHI_BB= DPHI
		return L_BB, PHI_BB, DPHI_BB
