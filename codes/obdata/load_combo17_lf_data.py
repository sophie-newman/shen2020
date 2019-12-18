# COMBO-17 (Wolf et al.)
#   measuring M_1450
#
from data import *
import numpy as np

def load_combo17_lf_data(z): # L_BB, PHI_BB, DPHI_BB, z
	if ((z <= 1.2) or (z > 4.8)): return False
	else:
		if ((z > 1.2) and (z <= 1.8)): WHICH_BLOCK = 0
		if ((z > 1.8) and (z <= 2.4)): WHICH_BLOCK = 1
		if ((z > 2.4) and (z <= 3.0)): WHICH_BLOCK = 2
		if ((z > 3.0) and (z <= 3.6)): WHICH_BLOCK = 3
		if ((z > 3.6) and (z <= 4.2)): WHICH_BLOCK = 4
		if ((z > 4.2) and (z <= 4.8)): WHICH_BLOCK = 5
		
		N_Z_BINS = 6
		N_LUM_BINS = 6
		M_1450 = np.array([-23.0,-24.0,-25.0,-26.0,-27.0,-28.0])
		M_B = M_1450 + 1.75 
		
		n_qso = np.zeros((N_LUM_BINS,N_Z_BINS),dtype=np.int32)
		phi   = np.zeros((N_LUM_BINS,N_Z_BINS),dtype=np.float32)
		dphi  = np.zeros((N_LUM_BINS,N_Z_BINS),dtype=np.float32)
		
		n_qso[:,0] = np.array([7,16,13,6,5,1])
		n_qso[:,1] = np.array([9,22,17,10,5,1])
		n_qso[:,2] = np.array([0,17,14,7,3,3])
		n_qso[:,3] = np.array([0,5,5,3,3,0])
		n_qso[:,4] = np.array([0,3,6,1,1,0])
		n_qso[:,5] = np.array([0,0,0,3,1,0])
		
		phi[:,0] = np.array([-5.78,-5.58,-5.69,-6.02,-6.10,-6.80])
		phi[:,1] = np.array([-5.53,-5.36,-5.60,-5.85,-6.15,-6.85])
		phi[:,2] = np.array([ 0.00,-5.42,-5.69,-6.01,-6.38,-6.38])
		phi[:,3] = np.array([ 0.00,-5.85,-6.10,-6.36,-6.36, 0.00])
		phi[:,4] = np.array([ 0.00,-6.02,-6.02,-6.81,-6.81, 0.00])
		phi[:,5] = np.array([ 0.00, 0.00, 0.00,-6.33,-6.79, 0.00])
		
		dphi[:,0] = np.array([ 0.10, 0.10, 0.15, 0.17, 0.50, 1.00])
		dphi[:,1] = np.array([ 0.10, 0.10, 0.13, 0.17, 0.50, 0.50])
		dphi[:,2] = np.array([ 0.00, 0.10, 0.10, 0.15, 0.25, 0.30])
		dphi[:,3] = np.array([ 0.00, 0.20, 0.30, 0.30, 0.50, 0.00])
		dphi[:,4] = np.array([ 0.00, 0.17, 0.17, 0.50, 0.50, 0.00])
		dphi[:,5] = np.array([ 0.00, 0.00, 0.00, 0.35, 0.60, 0.00])
		
		P0  = phi[:,WHICH_BLOCK]
		N0  = n_qso[:,WHICH_BLOCK]
		D0  = dphi[:,WHICH_BLOCK]
		ok  = (N0 > 0)
		PHI_BB = P0[ok] + np.log10(2.5) +3.*np.log10(7./6.5) # converts to per magnitude
		M_BB   = M_B[ok] + 5.*np.log10(7./6.5)
		DPHI_BB= D0[ok]
		L_BB   = 0.4*(M_sun_Bband_AB - M_BB)
		return L_BB, PHI_BB, DPHI_BB

# Function to return the analytical Wolf et al. (2003) luminosity function 
#   for a list of B-band luminosities L0_list (in SOLAR luminosities) at redshift z
#   (for an Omega_M = 0.3, Omega_Lambd = 0.7 cosmology)
#
def return_combo17_lf_fitted(L0_list,z):
	A0 = -5.620
	A1 = 0.1845
	A2 = -0.02652
	M0_145 = -25.0
	# convert to a B-band M0
	M0_B = M0_145 + 1.75 +5.*np.log10(7./6.5)
	B1 = 1.3455
	B2 = -80.845
	B3 = 127.32
	C1 = 0.3599
	C2 = -15.574
	
	zeta = np.log10((1.+z)/(1.+2.))
	Mz_145 = M0_145 + B1*zeta + B2*zeta**2 + B3*zeta**3

	Mz_B   = M0_B   + B1*zeta + B2*zeta**2 + B3*zeta**3
	M_B_t  = M_sun_Bband_AB - 2.5*(L0_list)
	mu     = M_B_t - Mz_B
	logphi = A0 + A1*mu + A2*mu*mu		# PLE


	mu     = M_B_t - Mz_B
	logphi = A0 + A1*mu + A2*mu*mu + C1*zeta + C2*zeta*zeta	# PDE

	PHI    = logphi + np.log10(2.5) + 3.*np.log10(7./6.5) 
	return PHI
