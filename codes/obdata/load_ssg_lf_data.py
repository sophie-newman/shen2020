# SSG -- broad range of redshift -- rest frame B-band
#
from data import * 
import numpy as np

def load_ssg_lf_data(z): # L_BB, PHI_BB, DPHI_BB, z
	if ((z <= 2.75) or (z > 4.75)): return False
	else:		
		pc = np.log10(2.5 * phi_correct_cosmology(3.75)*((7./5.)**3))
		M_B_o = np.array([-25.8,-26.3,-26.8,-27.5])+ 5.*np.log10(7./5.) - 2.5*np.log10(lum_correct_cosmology(3.75)) 
		 # correct from h=0.5, Om=1, Ol=0
		phi_o = np.array([-6.90, -7.30, -7.450, -8.20])  + pc
		phi_p = np.array([-6.80, -7.20, -7.375, -8.05])  + pc
		phi_m = np.array([-7.05, -7.40, -7.550, -8.40])  + pc
		dphi_o= 0.5*(phi_p-phi_m)
		
		L_BB   = 0.4*(M_sun_Bband_AB-M_B_o)
		PHI_BB = phi_o
		DPHI_BB= dphi_o
		return L_BB, PHI_BB, DPHI_BB
