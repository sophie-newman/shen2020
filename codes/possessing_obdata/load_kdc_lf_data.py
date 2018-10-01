# KDC -- rest B-band
#
from data import *
import numpy as np

def load_kdc_lf_data(z): # L_BB, PHI_BB, DPHI_BB, z
	if ((z < 4.0) or (z > 4.5)): return False
	else:	
		M_B_o = np.array([-27.0,-28.0]) + 5.*np.log10(7./5.) - 2.5*np.log10(lum_correct_cosmology(4.25))  # correct from h=0.5, Om=1, Ol=0
		phi_o  = np.array([5.7e-9,5.5e-10]) * 2.5 * (0.7/0.5)**3 * phi_correct_cosmology(4.25)
		dphi_o = np.array([2.6e-9,3.6e-10]) * 2.5 * (0.7/0.5)**3 * phi_correct_cosmology(4.25)

		L_BB   = 0.4*(M_sun_Bband_AB-M_B_o)
		PHI_BB = np.log10(phi_o)
		DPHI_BB= dphi_o/phi_o/np.log(10.)
		return L_BB, PHI_BB, DPHI_BB
