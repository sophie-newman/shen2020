#
# Hunt et al. (possibly being updated with new data)
#   Measuring M_1450 (converted to M_B following DR3 data)
#
from data import *
import numpy as np

def load_hunt_lf_data(z): # L_BB, PHI_BB, DPHI_BB, z
	if (z <= 2.5) or (z >= 3.5): return False
	else:
		M_1450 = np.array([ -20.4,  -21.4,  -22.4,  -23.4,  -24.4,  -25.4]) + 5.*np.log10(7./5.) - 2.5*np.log10(lum_correct_cosmology(3.0)) 						# correct from h=0.5, Om=1, Ol=0
		PHI    = np.array([1.2e-6, 2.8e-6, 3.5e-7, 1.5e-6, 5.0e-7, 9.1e-7]) * 2.5 * (7./5.)**3 * phi_correct_cosmology(3.0)
		DPHI   = np.array([   0.4,   0.25,    0.7,    0.2,    0.5,    0.3])
		#M_BB   = M_1450 + 1.75
		M_BB   = M_1450 - 0.776
		L_BB   = 0.4*(M_sun_Bband_AB-M_BB)
		PHI_BB = np.log10(PHI)
		DPHI_BB= DPHI
		return L_BB, PHI_BB, DPHI_BB
