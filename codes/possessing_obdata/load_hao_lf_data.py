# Hao et al. local SDSS -- choose their OIII measure, as its least 
#   contaminated by star formation -- optical technically, but really an isotropic measure
#   that seems to be correlated with L(2-10kev), and even better correlated with harder 
#   X-rays (Mulchaey 1994, Heckman et al. 2005 ApJ 634, 161)
#
#  remember, should ignore obscuration when doing the narrow line fits, 
#  as the Heckman study demonstrates
from data import *
import numpy as np

def load_hao_lf_data(z): # L_HX, PHI_HX, DPHI_HX, z
	if (z > 0.15): return False
	else:
		LO3 = np.array([2.0e5,  3.0e5,  4.4e5,  6.5e5,  9.5e5,  1.4e6,  2.1e6,  3.1e6,  4.8e6,  7.0e6,
			   1.0e7,  1.5e7,  2.0e7,  3.1e7,  5.0e7,  7.0e7,  1.0e8,  1.6e8,  2.2e8,  3.5e8,  5.0e8])
		PO31= np.array([1.1e-2, 5.0e-3, 2.5e-3, 1.4e-3, 1.3e-3, 5.1e-4, 3.0e-4, 2.6e-4, 2.8e-4, 1.9e-4,
		       1.3e-4, 8.0e-5, 5.5e-5, 2.9e-5, 1.4e-5, 4.5e-6, 7.0e-6, 1.4e-6, 5.0e-8, 1.4e-6, 2.0e-8])
		PO32= np.array([1.1e-2, 5.0e-3, 2.5e-3, 1.4e-3, 0.5*1.3e-3, 0.5*5.1e-4, 0.5*3.0e-4, 0.5*2.6e-4,
		       0.3*2.8e-4, 0.3*1.9e-4, 0.3*1.3e-4, 0.3*8.0e-5, 0.3*5.5e-5, 
		       0.5*2.9e-5, 0.5*1.4e-5, 4.5e-6, 2.0e-6, 9.0e-8, 1.8e-6, 0.1e-8, 0.1e-8])	
		       # Type 2 contribution to each bin
		
		PO3 = np.log10(PO31 + PO32)
		D00 = np.log10(3./2.)
		DO3 = np.array([   D00,    D00,    D00,    D00,    D00,    D00,    D00,    D00,    D00,    D00, 
		          D00,    D00,    D00,    D00,    0.2,    0.3,    0.3,   0.48,    0.4,   0.35,   0.8])
		# they adopt H0 = 100

		PO3 = PO3 + 3.*np.log10(hubble)			
		LO3 = LO3/hubble**2
	
		# Mulchaey calibration considerably older, for 2-10; 
		#   Heckman's more recent study has more (still only ~50), shows a problem -- 
		#   2-10 significantly affected by obscuration (compton-thick) at the lowest L, 
		#   and so very different ratio inferred from e.g. type 1 & type 2 ( & not clear 
		#   that the type 1's are completely unaffected). also forces splitting the 
		#   sample -> worse stats. But 3-20 seems clean - no type 1/2 dep't. Mean there 
		#   is a 2.15 dex offset, then from our model spectrum we can correct to 
		#   2-10 = (2/3) * (3-20); since we're interested in the intrinsic spectrum that 
		#   should be ok?
	
		L_O3_over_L_2to10keV = 0.01 #4 	# Mulchaey calibration
		#L_O3_over_L_2to10keV = 10**(-2.00015) / (0.661) # very similar -- at this point, just splitting hairs
		#L_O3_over_L_2to10keV = 10**(-1.59)  # direct for unobscured type-1s, less certain (could be ~0.2 dex larger)
		
		# Heckman calibration with our 2-10/3-20 intrinsic spectrum
		L_HX   = np.log10( LO3 / L_O3_over_L_2to10keV )
		PHI_HX = PO3	# assuming a linear proportionality for now
		DPHI_HX= DO3

		return L_HX, PHI_HX, DPHI_HX


