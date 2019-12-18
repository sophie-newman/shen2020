# 
# La Franca et al. 2005 data
#
from data import *
import numpy as np

def load_lafranca_lf_data(z): # L_HX, PHI_HX, DPHI_HX, z
	# determine which redshift interval its in
	if ((z > 0.0) and (z <= 0.5)): WHICH_BLOCK = 1
	if ((z > 0.5) and (z <= 1.0)): WHICH_BLOCK = 2
	if ((z > 1.0) and (z <= 1.5)): WHICH_BLOCK = 3
	if ((z > 1.5) and (z <= 2.5)): WHICH_BLOCK = 4
	if ((z > 2.5) and (z <= 3.5)): WHICH_BLOCK = 5
	if (z > 3.5): WHICH_BLOCK = 6
	
	if (WHICH_BLOCK == 6): return False		
	else:
		if (WHICH_BLOCK == 1):
			L_HX = np.array([42.5, 43.5, 44.5, 45.5])
			P_HX = np.array([-3.7, -4.6, -6.2, -9.05])
			D_HX = np.array([ 0.2,  0.2,  0.2,  0.8])
		if (WHICH_BLOCK == 2):
			L_HX = np.array([42.5, 43.5, 44.5, 45.5])
			P_HX = np.array([-3.25, -4., -5.45, -8.0])
			D_HX = np.array([ 0.2,  0.2,  0.2,  0.33])
		if (WHICH_BLOCK == 3):
			L_HX = np.array([42.5, 43.5, 44.5, 45.5])
			P_HX = np.array([-3.1, -4.05, -5., -7.0])
			D_HX = np.array([ 0.3,  0.2,  0.2,  0.3])
		if (WHICH_BLOCK == 4):
			L_HX = np.array([43.5, 44.5, 45.5, 46.5])
			P_HX = np.array([-4.1, -4.9, -6.9, -8.95])
			D_HX = np.array([ 0.2,  0.2,  0.3,  0.5])
		if (WHICH_BLOCK == 5):
			L_HX = np.array([44.25, 45.0])
			P_HX = np.log10(np.array([2.0 * 8.0e-6, 1.0 * 2.0e-6]))
			D_HX = np.array([ 0.2, 0.25])
		L_HX = (L_HX - L_solar)
		PHI_HX  = P_HX
		DPHI_HX = D_HX
		return L_HX, PHI_HX, DPHI_HX

def return_lafranca_lf_fitted(L0_list,z):	
	L_X_star_0 = 10**(44.25 - L_solar) # (converts to solar luminosities)
	GAMMA_1 = 1.01
	GAMMA_2 = 2.38
	
	A_0 = 1.21e-6	# Normalization - in AGNs/Mpc^3/log(L)	
	p1 = 4.62		# Account for redshift evolution
	p2 = -1.15
	L_a = 10**(45.74 - L_solar)
	z_cut_star = 2.49
	alpha = 0.20

	z_cut = 0.0*L0_list + z_cut_star
	zcut_ldde = (10**(L0_list) < L_a)
	if np.count_nonzero(zcut_ldde)!=0: 
		z_cut[zcut_ldde] = z_cut_star * (10**(L0_list[zcut_ldde])/L_a)**alpha
 
	e_z = 0.0*L0_list + (1.0 + z)**p1
	zl = (z_cut <= z)
	if np.count_nonzero(zl)!=0:
		e_z[zl] = ((1.+z_cut[zl])**(p1)) * ((1.+z)/(1.+z_cut[zl]))**p2
	
	x = 10**(L0_list) / L_X_star_0
	PHI = A_0 * e_z / (x**GAMMA_1 + x**GAMMA_2)

	return np.log10(PHI)
