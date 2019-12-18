#
#  Brown et al. 2005; 1<z<5 combined data set - but binned only given at z=2
#		8 micron data
#
#  potential issue -- the Brown et al. sample is only Type 1's. which even at the brightest 
#     luminosities will be affected sign. by obsc. since lower-lum heavily obscured 
#     sources can be up-scattered up here (see e.g. the matute et al. LFs)
#
from data import *
import numpy as np

def load_brown_lf_data(z): # L_IR, PHI_IR, DPHI_IR, z
	if (z <= 1.5) or (z > 2.5): return False
	else:
		M_8  = np.array([-27.2, -27.6, -28.0, -28.4, -28.8, -29.2, -29.6, -30.0, -30.4, -31.2])
		P_8  = np.array([4.3e-6,2.4e-6,3.0e-6,1.2e-6,6.1e-7,2.6e-7,2.4e-7,9.2e-8,4.7e-8,2.3e-8])
		D_8  = np.array([ 0.21,  0.18,  0.10,  0.15,  0.15,  0.18,  0.18,  0.35,  0.40,  0.45])

		# use below if converting to another band
		#M_B = M_8 + 3.08
		#L_B = 10^(DOUBLE(0.4*(M_sun_BB(0) - M_B)))		
		#L_HX   = alog10(L_B) + alog10(10.0)		

		L_8 = 10**(0.4*(-28.5 - M_8)) * 4.11e45	# conversion given therein
		L_8 = L_8 / np.power(10.,L_solar)
		L_15_over_L_8 = 10**(+0.060001373)	# from Hatziminaoglou model spectrum
		#L_15_over_L_8 = 10**(+0.0263996)	# from Hatziminaoglou model spectrum
		L_15 = L_15_over_L_8 * L_8
		L_IR = np.log10(L_15)
		PHI_IR = np.log10(P_8) + np.log10(2.5) #+ alog10(1./(1.-0.55))
		# correction for mean Type-2 fraction (NOTE: not 
		#   obscured in the IR, but not optical qso so 
		#   no redshift data in this limited sample
		DPHI_IR= D_8
		return L_IR, PHI_IR, DPHI_IR

def return_brown_lf_fitted(L0_list,z):	
	# given L in L_solar at 15 microns, convert it 
	L_15_over_L_8 = 10**(+0.060001373)	# from Hatziminaoglou model spectrum
		#L_15_over_L_8 = 10^(+0.0263996)	# from Hatziminaoglou model spectrum
	L_8 = (10.**L0_list * 10.**L_solar )/L_15_over_L_8
	M_8 = -28.5 - 2.5*np.log10(L_8/4.11e45)

	phi_29_2 = 4.46e-7 * 2.5
	alpha = -2.75
	k1    =  1.15
	k2    = -0.34
	k3    =  0.03
	zpeak =  2.56

	M_z   = -29.0 - 2.5*(k1*z + k2*z**2 + k3*z**3) + 2.5*(2.*k1 + 4.*k2 + 8.*k3)
	phi   = phi_29_2 * 10**(0.4*(alpha+1.)*(M_z-M_8))
	
	return np.log10(phi) #+ np.log10(1./(1.-0.55))	
		# correction for mean Type-2 fraction (NOTE: not 
		#   obscured in the IR, but not optical qso so 
		#   no redshift data in this limited sample
