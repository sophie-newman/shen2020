# 
# Aird et al. 2015
#
from data_copy import *
import numpy as np 

def load_aird15_lf_data_softsel(z): # L_HX, PHI_HX, DPHI_HX
	if (z <= 0) or (z > 7.): return False
	else:	
		filename = datapath+'Aird15_Fig7.dat'
		data = np.genfromtxt(filename,names=True)

		id = (data["z_min"]< z) & (data["z_max"]>=z) & (data['LFsoft']!=0) & (data['LFsoft_lo']>0)

		L_HX  = (data['loglx_lo'][id] + data['loglx_hi'][id])/2.
		P_HX  = np.log10(data['LFsoft'][id])
		P_HX_up  = np.log10(data['LFsoft_hi'][id])
		P_HX_down= np.log10(data['LFsoft_lo'][id])
		D_HX = (P_HX_up - P_HX_down)/2.

		PHI_HX = P_HX
		DPHI_HX= D_HX
		L_HX = (L_HX - L_solar)

		return L_HX, PHI_HX, DPHI_HX

def load_aird15_lf_data_hardsel(z): # L_HX, PHI_HX, DPHI_HX
	if (z <= 0) or (z > 7.): return False
	else:	
		filename = datapath+'Aird15_Fig7.dat'
		data = np.genfromtxt(filename,names=True)

		id = (data["z_min"]< z) & (data["z_max"]>=z) & (data['LFhard']!=0) & (data['LFhard_lo']>0)

		L_HX  = (data['loglx_lo'][id] + data['loglx_hi'][id])/2.
		P_HX  = np.log10(data['LFhard'][id])
		P_HX_up  = np.log10(data['LFhard_hi'][id])
		P_HX_down= np.log10(data['LFhard_lo'][id])
		D_HX = (P_HX_up - P_HX_down)/2.

		PHI_HX = P_HX
		DPHI_HX= D_HX
		L_HX = (L_HX - L_solar)

		return L_HX, PHI_HX, DPHI_HX
	