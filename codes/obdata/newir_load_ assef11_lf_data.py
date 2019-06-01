# 
#  Assef et al. 2011
#
from data import *
import numpy as np 

def load_assef11_lf_data(z): # M_J, PHI_J, DPHI_J
	if (z <= 0) or (z > 5.85): return False
	else:	
		filename = datapath+'Assef2011.dat'
		data = np.genfromtxt(filename,names=True)

		id = (data["z_min"]< z) & (data["z_max"]>=z) & (data['label']==1) & (data['lf']!=data['lf_sigma'])

		M_J  = data['MJ'][id]
		P_J  = np.log10(data['lf'][id])
		P_J_up  = np.log10(data['lf'][id] + data['lf_sigma'][id])
		P_J_down= np.log10(data['lf'][id] - data['lf_sigma'][id])
		D_J = (P_J_up - P_J_down)/2.

		PHI_J = P_J
		DPHI_J= D_J

		z_c = data["z_mean"][id][0]
		M_J = M_J - 2.5*np.log10(lum_correct_cosmo_flexible(z_c, 0.73, 0.3))
		PHI_J = PHI_J + np.log10(phi_correct_cosmo_flexible(z_c, 0.73, 0.3))

		L_J = (10**(-0.4*M_J) * L_AB) * (con.c.value/1.25e-4)

		L_IR = L_J + 0.11928 # convert to 15 micron
		PHI_IR = PHI_J + np.log10(2.5)
		DPHI_IR = DPHI_J

		return L_IR, PHI_IR, DPHI_IR
