# Matute et al. 2006 -- 15 micron local and z=1.2
#
#   ADD THE TYPE 2's? -- weak contribution, although up to factor ~2 at low-L
#
from data import *
import numpy as np 

def load_matute_lf_data(z): # L_IR, PHI_IR, DPHI_IR, z
	if (((z > 0.2) and (z < 0.5)) or (z >= 2.0)): return False
	else:
		if (z <= 0.2):
			L_15 = np.array([42.6,  43.0,  43.4,  43.8,  44.2,  44.6,  45.0,  45.4,  45.8,  46.2])
			P_15 = np.array([-4.10, -4.30, -4.90, -5.25, -5.74, -6.16, -7.12, -8.05, -8.70, -8.97]) 
			D_15 = np.array([ 0.25,  0.15,  0.25,  0.25,  0.25,  0.20,  0.40,  0.50,  0.55,  0.55])
		if ((z >= 0.5) and (z <= 2.0)):
			L_15 = np.array([43.8,  44.2,  44.6,  45.0,  45.4,  45.8,  46.2,  46.6])
			P_15 = np.array([-4.72, -5.25, -5.30, -5.30, -6.13, -6.70, -7.50, -8.52])
			D_15 = np.array([ 0.30,  0.30,  0.20,  0.25,  0.20,  0.20,  0.25,  0.55])	
		
		# use below if converting to other bands
		# log_L_15_over_L_R = 0.23
		# L_R = L_15 - log_L_15_over_L_R - (alog10(3.9)+33.0)
		# L_HX   = L_R + alog10(17.5) - 2.*alog10(7./7.5)
		
		L_IR   = (L_15-L_solar) - 2.5*np.log10(7./7.5) #+ 0.25
		PHI_IR = P_15 + np.log10(2.5) + 3.*np.log10(7./7.5)		 #- 0.2
		DPHI_IR= D_15
		return L_IR, PHI_IR, DPHI_IR
