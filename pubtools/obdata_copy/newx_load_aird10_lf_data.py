# 
# Aird et al. 2010
#
from data_copy import *
import numpy as np 

def load_aird10_lf_data(z): # L_HX, PHI_HX, DPHI_HX
	if (z <= 0) or (z > 3.5): return False
	else:	
		filename = datapath+'Aird10_Fig9.dat'
		data = np.genfromtxt(filename,names=True)

		id = (data["z_min"]< z) & (data["z_max"]>=z) & (data['lf']!=0)

		L_HX  = data['loglx'][id]
		P_HX  = np.log10(data['lf'][id])
		P_HX_up  = np.log10(data['lf'][id] + data['lf_hi'][id])
		P_HX_down= np.log10(data['lf'][id] - data['lf_lo'][id])
		D_HX = (P_HX_up - P_HX_down)/2.

		
		PHI_HX = P_HX
		DPHI_HX= D_HX
		L_HX = (L_HX - L_solar)

		return L_HX, PHI_HX, DPHI_HX
