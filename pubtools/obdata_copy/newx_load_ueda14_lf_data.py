# Ueda et al. 2015
# measured 2-10 keV
from data_copy import *
import numpy as np

def load_ueda14_lf_data(z): 
	if (z <= 0.002) or (z > 5.0): return False
	else:	
		filepath = datapath+'ueda2014/'
		if (z > 0.002) and (z <= 0.2): annex = "-z0.002-0.2.dat"
		elif (z > 0.2) and (z <= 0.4): annex = "-z0.2-0.4.dat"
		elif (z > 0.4) and (z <= 0.6): annex = "-z0.4-0.6.dat"
		elif (z > 0.6) and (z <= 0.8): annex = "-z0.6-0.8.dat"
		elif (z > 0.8) and (z <= 1.0): annex = "-z0.8-1.0.dat"
		elif (z > 1.0) and (z <= 1.2): annex = "-z1.0-1.2.dat"
		elif (z > 1.2) and (z <= 1.6): annex = "-z1.2-1.6.dat"
		elif (z > 1.6) and (z <= 2.0): annex = "-z1.6-2.0.dat"
		elif (z > 2.0) and (z <= 2.4): annex = "-z2.0-2.4.dat"
		elif (z > 2.4) and (z <= 3.0): annex = "-z2.4-3.0.dat"
		elif (z > 3.0) and (z <= 4.0): annex = "-z3.0-4.0.dat"
		elif (z > 4.0) and (z <= 5.0): annex = "-z4.0-5.0.dat"

		data_HX = np.genfromtxt(filepath+"hard"+annex, names=['L', 'phi', 'loerr', 'uperr'])
		data_SX = np.genfromtxt(filepath+"soft"+annex, names=['L', 'phi', 'loerr', 'uperr'])

		L_X  = np.append(data_HX["L"], data_SX["L"])
		phi  = np.append(data_HX["phi"], data_SX["phi"])
		loerr  = np.append(data_HX["loerr"], data_SX["loerr"])
		uperr  = np.append(data_HX["uperr"], data_SX["uperr"])	

		L_X = L_X - L_solar		 
		PHI_X = np.log10(phi)

		DPHI_X = ((np.log10(phi+uperr)-np.log10(phi)) + (np.log10(phi)-np.log10(phi-loerr)))/2.

		PHI_X = PHI_X + absorption_correction_xray(L_X, z)
		return L_X, PHI_X, DPHI_X

