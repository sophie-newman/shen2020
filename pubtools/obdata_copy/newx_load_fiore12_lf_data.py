# 
# Fiore et al. 2012 data
# Chandra GOODS field, ERS field
#
from data_copy import *
import numpy as np 

def load_fiore12_lf_data(z): # L_HX, PHI_HX, DPHI_HX, z
	# determine which redshift interval its in
	if (z > 3.) and (z <= 4.):
		L_HX_high=np.array([43.50,44.00,44.50])
		L_HX_low =np.array([42.75,43.50,44.00])
		P_HX = np.log10(np.array([4.3e-5,4.0e-5,1.6e-5]))
		D_HX_high = np.array([2.9,1.7,1.3])*1e-5
		D_HX_low  = np.array([1.9,1.2,0.8])*1e-5
	elif (z > 4.) and (z <= 5.):
		L_HX_high=np.array([44.00])
		L_HX_low =np.array([43.00])
		P_HX = np.log10(np.array([2.1e-5]))
		D_HX_high = np.array([0.9])*1e-5
		D_HX_low  = np.array([0.7])*1e-5
	elif (z > 5.8) and (z < 7.5):
		L_HX_high=np.array([44.50])
		L_HX_low =np.array([43.50])
		P_HX = np.log10(np.array([0.66e-5]))
		D_HX_high = np.array([1.1])*1e-5
		D_HX_low  = np.array([0.5])*1e-5
	else: return False

	L_HX = (L_HX_low+L_HX_high)/2.
	L_HX = (L_HX - L_solar) 
	PHI_HX  = P_HX
	DPHI_HX = ( np.log10(10**P_HX + D_HX_high) - np.log10(10**P_HX - D_HX_low) )/2.
	return L_HX, PHI_HX, DPHI_HX