# 2SLAQ (Richards et al. -- supplants Boyle et al. & Croom et al.)
# rest-frame B-band
from data_copy import *
import numpy as np

def load_2slaq_lf_data(z): 
	#determine which redshift interval its in (bins follow Croom 2004)
	if (z <= 0.4): WHICH_BLOCK = 0
	if (z > 0.40) and (z <= 0.68): WHICH_BLOCK = 1
	if (z > 0.68) and (z <= 0.97): WHICH_BLOCK = 2
	if (z > 0.97) and (z <= 1.25): WHICH_BLOCK = 3
	if (z > 1.25) and (z <= 1.53): WHICH_BLOCK = 4
	if (z > 1.53) and (z <= 1.81): WHICH_BLOCK = 5
	if (z > 1.81) and (z <= 2.10): WHICH_BLOCK = 6
	if (z > 2.1): WHICH_BLOCK = 7

	if (WHICH_BLOCK == 0) or (WHICH_BLOCK == 7):
		return False
	else:
		data=np.genfromtxt(datapath+"load_2slaq_lf_data.dat",names=['z','M_B','phi','dphi','nqso','bin_flags'])

		N_LUM_BINS = 34
		N_Z_BINS   = 6
		
		z_list = np.zeros((N_LUM_BINS,N_Z_BINS))
		M_B    = z_list * 0
		phi    = z_list * 0
		dphi   = z_list * 0
		nqso   = z_list * 0
		bin_flags = np.zeros((N_LUM_BINS,N_Z_BINS),dtype=np.int32)
	
		for i in range(N_Z_BINS):
				z_list[:,i] = data['z'][i*N_LUM_BINS:(i+1)*N_LUM_BINS]
				M_B[:,i] = data['M_B'][i*N_LUM_BINS:(i+1)*N_LUM_BINS]
				phi[:,i] = data['phi'][i*N_LUM_BINS:(i+1)*N_LUM_BINS]
				dphi[:,i] = data['dphi'][i*N_LUM_BINS:(i+1)*N_LUM_BINS]
				nqso[:,i] = data['nqso'][i*N_LUM_BINS:(i+1)*N_LUM_BINS]
				bin_flags[:,i] = data['bin_flags'][i*N_LUM_BINS:(i+1)*N_LUM_BINS]
		
		M_B = M_B[:,WHICH_BLOCK-1] + 0.05	# bJ = g + 0.05 from Fukugita et al. (R06 DR3 ppr)
		phi = phi[:,WHICH_BLOCK-1]
		dphi = dphi[:,WHICH_BLOCK-1]
		bin_flags = bin_flags[:,WHICH_BLOCK-1]
		ok = (bin_flags >= 0)
		
		L_BB = 0.4*(M_sun_Bband_AB-M_B[ok]) 
		PHI_BB = np.log10(phi[ok]) + np.log10(2.5)
		DPHI_BB = np.log10(1.+dphi[ok]/phi[ok])
		return L_BB, PHI_BB, DPHI_BB

# Function to return the analytical Richards et al. (2005) luminosity function 
#   for a list of B-band luminosities L0_list (in SOLAR luminosities) at redshift z
#   (for an Omega_M = 0.3, Omega_Lambd = 0.7 cosmology)
#
def return_2slaq_lf_fitted(L0_list,z):
	OMEGA_MATTER = 0.3	# Cosmology - not important
	OMEGA_LAMBDA = 0.7  
	N_QSO = 6180		# Number used in sample
	P_KS  = 0.075		# K-S probability of this fit

	ALPHA = 3.28
	BETA  = 1.78

	M_B_star_0 = -22.68	 + 0.05 # Break magnitude at z=0, with b/g conversion
	L_B_star_0 = 10**((M_sun_Bband_AB-M_B_star_0)/2.5)
	k1 =  1.37			# Determines the polynomial evolution of the break
	k2 = -0.32			#   magnitude, as below
	
	PHI_STAR = 5.96e-7	# Normalization - in AGNs/Mpc^3/mag
	# Convert this to AGNs/Mpc^3/log(L) for luminosity units, 
	#   simply using dmag = -2.5 d(log(L)), but then get cancelling
	#   factor of 2.5 from other side of the equation
	PHI_STAR_L = PHI_STAR * 2.5
	
	# Account for redshift evolution
	# M_B_star = M_B_star_0 - 2.5 * (k1*z + k2*z*z), so have 
	L_B_star = L_B_star_0 * 10**(k1*z + k2*z*z)

	x = 10**(L0_list) / L_B_star
	PHI = PHI_STAR_L / (x**(ALPHA-1.) + x**(BETA-1.))
	return np.log10(PHI)
