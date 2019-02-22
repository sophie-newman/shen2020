# 
#  Khorunzhev et al. 2018
#
from data import *
import numpy as np 

def load_khorunzhev18_lf_data_vitobins(z): # L_HX, PHI_HX, DPHI_HX
	if (z <= 3) or (z > 5.): return False
	else:	
		filename = datapath+'Khorunzhev18_Vitobins.dat'
		data = np.genfromtxt(filename,names=True)

		id = (data["z_min"]< z) & (data["z_max"]>=z)

		L_HX  = np.log10(data['lx'][id])
		P_HX  = np.log10(data['lf'][id])
		P_HX_up  = np.log10(data['lf'][id] + data['lf_err_hi'][id])
		P_HX_down= np.log10(data['lf'][id] - data['lf_err_lo'][id])
		D_HX = (P_HX_up - P_HX_down)/2.

		PHI_HX = P_HX
		DPHI_HX= D_HX
		L_HX = (L_HX - L_solar)

		L_HX = L_HX + np.log10(lum_correct_cosmo_flexible(z_c, 0.7, 0.27))
		PHI_HX = PHI_HX + np.log10(phi_correct_cosmo_flexible(z_c, 0.7, 0.27))

		return L_HX, PHI_HX, DPHI_HX


def load_khorunzhev18_lf_data_z345bins(z): # L_HX, PHI_HX, DPHI_HX
	if (z <= 3) or (z > 5.): return False
	else:	
		filename = datapath+'Khorunzhev18_z345bins.dat'
		data = np.genfromtxt(filename,names=True)

		id = (data["z_min"]< z) & (data["z_max"]>=z)

		L_HX  = np.log10(data['lx'][id])
		P_HX  = np.log10(data['lf'][id])
		P_HX_up  = np.log10(data['lf'][id] + data['lf_err_hi'][id])
		P_HX_down= np.log10(data['lf'][id] - data['lf_err_lo'][id])
		D_HX = (P_HX_up - P_HX_down)/2.

		PHI_HX = P_HX
		DPHI_HX= D_HX
		L_HX = (L_HX - L_solar)

		L_HX = L_HX + np.log10(lum_correct_cosmo_flexible(z_c, 0.7, 0.27))
		PHI_HX = PHI_HX + np.log10(phi_correct_cosmo_flexible(z_c, 0.7, 0.27))

		return L_HX, PHI_HX, DPHI_HX
	