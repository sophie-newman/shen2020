# 
# Silverman et al. 2008
#
from data import *
import numpy as np 

def load_silverman08_lf_data(z): # L_HX, PHI_HX, DPHI_HX
	if ((z <= 0.2) and (z > 5.5)): return False
	else:
		L_HX = np.array([42.74510])
		P_HX = np.array([-4.30889])
		P_HX_up  = np.array([-3.79011])
		P_HX_down= np.array([-5.05974])
		D_HX = (P_HX_up - P_HX_down)/2.
		L_HX = (L_HX - L_solar) 
		PHI_HX  = P_HX
		DPHI_HX = D_HX
		return L_HX, PHI_HX, DPHI_HX

#not completed