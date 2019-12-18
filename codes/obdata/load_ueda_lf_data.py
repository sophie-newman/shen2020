# 
# Ueda et al. 2003 data
#
from data import *
import numpy as np

def load_ueda_lf_data(z): # L_HX, PHI_HX, DPHI_HX, z
	# determine which redshift interval its in
	if ((z > 0.0) and (z <= 0.2)): WHICH_BLOCK = 1
	if ((z > 0.2) and (z <= 0.4)): WHICH_BLOCK = 2
	if ((z > 0.4) and (z <= 0.8)): WHICH_BLOCK = 3
	if ((z > 0.8) and (z <= 1.6)): WHICH_BLOCK = 4
	if ((z > 1.6) and (z <= 3.0)): WHICH_BLOCK = 5
	if (z > 3.0): WHICH_BLOCK = 6
	
	if (WHICH_BLOCK == 6): return False
	else:
		BEGIN = False
		filename=datapath+'load_ueda_lf_data.dat'
		with open(filename, 'r') as f:
			i=0
			j=0
			for line in f.readlines():
				elements=line.split()
				if (BEGIN == True) and (j<N_LUM_BINS):
					L_HX[j] = float(elements[2])
					PHI_HX[j] = float(elements[3])
					DPHI_HX[j] = float(elements[4])
					j+=1
		 		if (len(elements[0]) == 1) and ( int(elements[0])==WHICH_BLOCK ):
					BEGIN = True
					N_LUM_BINS = int(elements[1])
					L_HX = np.zeros(N_LUM_BINS)
					PHI_HX  = np.zeros(N_LUM_BINS)
					DPHI_HX = np.zeros(N_LUM_BINS)
				i+=1

		DPHIm_HX = -1.0*DPHI_HX
		DPHIp_HX = DPHI_HX
		L_HX = (L_HX - L_solar)
		return L_HX, PHI_HX, DPHI_HX

def return_ueda_lf_fitted(L0_list,z):
	OMEGA_MATTER = 0.3	# Cosmology - not important
	OMEGA_LAMBDA = 0.7  
	N_QSO = 141			# Number used in sample ??
	P_KS  = 0.9			# K-S probability of this fit

	GAMMA_1 = 0.86
	GAMMA_2 = 2.23
	
	#L_X_star_0 = 10**(43.94)
	L_X_star_0 = 2.1774e10		# (converts to solar luminosities)
	
	A_0 = 5.04e-6	# Normalization - in AGNs/Mpc^3/log(L)
	
	p1 = 4.23		# Account for redshift evolution
	p2 = -1.5
	#L_a = 10**44.6
	L_a = 1.0e11
	z_cut_star = 1.9
	alpha = 0.335

	z_cut = 0.0*L0_list + z_cut_star
	zcut_ldde = (10**L0_list < L_a)
	if np.count_nonzero(zcut_ldde)!=0:
		z_cut[zcut_ldde] = z_cut_star * (10**(L0_list[zcut_ldde])/L_a)**alpha
	
	e_z = 0.0*L0_list + (1.0 + z)**p1
	zl = (z_cut <= z)
	if np.count_nonzero(zl):
		e_z[zl] = ((1.+z_cut[zl])**(p1)) * ((1.+z)/(1.+z_cut[zl]))**p2
	
	x = 10**(L0_list) / L_X_star_0
	PHI = A_0 * e_z / (x**GAMMA_1 + x**GAMMA_2)

	return np.log10(PHI)
