# 
# Sazonov & Revnivtsev al. 2004 local data
#
from data_copy import *
import numpy as np

def load_sazrev_lf_data(z): # L_HX, PHI_HX, DPHI_HX, z
	# determine which redshift interval its in
	if (z > 0.15): return False
	else:
		L_HX = np.array([41.25,  41.75,  42.25,  42.75,  43.25,  43.75,  44.25,  44.75,  45.25,  45.75])
		P_HX = np.array([1.0e-3, 2.9e-4, 1.0e-4, 3.5e-5, 7.3e-6, 1.5e-6, 1.7e-7, 2.2e-8, 4.8e-10,1.4e-10])
		D_HX = np.array([ 0.4,    0.35,   0.30,   0.20,   0.20,   0.15,   0.20,   0.25,   0.35,   0.4])

	L_HX = (L_HX - L_solar) + np.log10( 0.66676 * (0.75/0.7)**2 )
	#L_HX = (L_HX - L_solar) + np.log10( (0.75/0.7)**2 ) ;+ np.log10(1./1.25)
	# convert the 2-30kev to 2-10kev (using our spectrum model in the X-rays)
	#  and from H0 = 75 to our 70
	PHI_HX  = np.log10(P_HX *1.43* (0.75/0.7)**3)	# includes the completeness correction they estimate
	DPHI_HX = D_HX

 	# Sazonov et al. 2006 (17-60 keV, 0.0 <~ z <~ 0.1-1.0 (ideal 0.15))
	l = np.array([ 40.25,  41.75,  42.25,  42.75,  43.25,  43.75,  44.25,  44.75,  45.25]) - 2.*np.log10(0.7/0.75) - np.log10(2.0)
	L_HX = np.append(L_HX,l-L_solar)
	p = np.log10(np.array([3.2e-2, 7.5e-4, 3.0e-4, 8.0e-5, 3.5e-5, 4.3e-6, 5.0e-7, 8.0e-8, 1.3e-9]) * (0.7/0.75)**3)
	PHI_HX = np.append(PHI_HX,p)
	d = np.array([0.4,    0.3,    0.24,   0.16,   0.09,   0.15,   0.17,   0.22,   0.39])
	DPHI_HX = np.append(DPHI_HX,d)

	return L_HX, PHI_HX, DPHI_HX
