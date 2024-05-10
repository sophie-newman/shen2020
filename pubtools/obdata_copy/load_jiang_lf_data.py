#
# Jiang et al. 2006 (low-L SDSS quasars)
#
from data_copy import *
import numpy as np

def load_jiang_lf_data(z): # L_BB, PHI_BB, DPHI_BB, z
	if (z <= 0.5) or (z > 3.6): return False
	else:
		if ((z > 0.5) and (z <= 1.0)): WHICH_BLOCK = 0
		if ((z > 1.0) and (z <= 1.5)): WHICH_BLOCK = 1
		if ((z > 1.5) and (z <= 2.0)): WHICH_BLOCK = 2
		if ((z > 2.0) and (z <= 2.5)): WHICH_BLOCK = 3
		if ((z > 2.5) and (z <= 3.0)): WHICH_BLOCK = 4
		if ((z > 3.0) and (z <= 3.6)): WHICH_BLOCK = 5
		
		if (WHICH_BLOCK == 0):
			M_g = np.array([-21.25, -21.75, -22.50, -23.25, -24.25])
			phi = np.array([2.3e-6, 2.3e-6, 1.4e-6, 1.0e-6, 4.3e-7])
			dphi= np.array([0.3000, 0.3000, 0.3000, 0.3000, 0.3500])
		if (WHICH_BLOCK == 1):
			M_g = np.array([-22.00, -22.50, -23.00, -23.50, -24.50])
			phi = np.array([3.5e-6, 3.0e-6, 2.3e-6, 2.2e-6, 7.8e-7])
			dphi= np.array([0.2500, 0.2400, 0.2400, 0.2300, 0.2400])
		if (WHICH_BLOCK == 2):
			M_g = np.array([-22.60, -23.10, -23.60, -24.20, -24.75, -25.20])
			phi = np.array([4.0e-6, 3.0e-6, 2.9e-6, 1.7e-6, 2.25e-6, 1.05e-6])
			dphi= np.array([0.2000, 0.1800, 0.1800, 0.1800, 0.18000, 0.20000])
		if (WHICH_BLOCK == 3):
			M_g = np.array([-23.00, -23.80, -24.15, -24.50, -25.40])
			phi = np.array([3.3e-6, 1.8e-6, 2.0e-6, 1.0e-6, 6.0e-7])
			dphi= np.array([0.4300, 0.3200, 0.3000, 0.3000, 0.4000])
		if (WHICH_BLOCK == 4):
			M_g = np.array([-23.35, -23.70, -24.60, -26.00])
			phi = np.array([3.0e-6, 1.4e-6, 8.1e-7, 3.0e-7])
			dphi= np.array([0.4250, 0.3680, 0.3580, 0.4840])
		if (WHICH_BLOCK == 5):
			M_g = np.array([-23.90, -24.60, -25.20])
			phi = np.array([1.0e-6, 7.9e-7, 4.3e-7])
			dphi= np.array([0.3520, 0.3500, 0.4820])

		phi = np.log10(phi*2.5) # convert to per log(L)
		dphi= dphi/2. # numbers above for full +/- 1sigma range
		M_B = M_g + 0.05 # bJ = g + 0.05 from Fukugita et al. (R06 DR3 ppr)
		
		L_BB   = 0.4*(M_sun_Bband_AB-M_B)
		PHI_BB = phi
		DPHI_BB= dphi
		return L_BB, PHI_BB, DPHI_BB

#
# here are the fit params including this data in the fitting procedure :: 
#   although the fit is technically changed, looking at it, it's clear that 
#   it's more or less *exactly* the same (like 0.01 dex change at most) at 
#   all luminosities of interest for our fitting/comparisons, so 
#   we know adding this data won't change our results
# {P0=-4.80530; P1=13.0073; P2=0.629385; P3=-11.7977; P4=-14.3867; 
#  P5=0.406049; P6=-0.649790; P7=2.17013;  P8=1.46220;  P9=-0.800668; 
#  P10=0.;P11=0.;P12=0.;P13=0.;P14=0.;}
#
