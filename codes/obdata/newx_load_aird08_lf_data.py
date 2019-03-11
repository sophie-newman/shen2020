# 
# Aird et al. 2008
#
from data import *
import numpy as np 

def load_aird08_lf_data(z): # L_HX, PHI_HX, DPHI_HX
	if ((z <= 2.5) or (z > 3.5)): return False
	else:
		L_HX = np.array([42.74510,43.24724,43.74163,44.24398,44.74645,45.49320,46.49968,47.49698])
		P_HX = np.array([-4.30889,-4.35098,-4.99375,-5.24062,-5.59672,-7.49603,-9.68262,-11.1047])
		P_HX_up  = np.array([-3.79011,-4.14621,-4.80262,-5.06315,-5.43289,-7.22299,-9.16383, -10.59958])
		P_HX_down= np.array([-5.05974,-4.56940,-5.18487,-5.45905,-5.76054,-7.82367,-10.41980,-11.85556])
		D_HX = (P_HX_up - P_HX_down)/2.
		L_HX = (L_HX - L_solar) 
		PHI_HX  = P_HX
		DPHI_HX = D_HX
		return L_HX, PHI_HX, DPHI_HX
