# Richards et al. SDSS DR3
# Measuring M_i(z=2) (outputs as M_B)
from data import *
import numpy as np

def load_sdss_dr3_lf_data(z): 
	if (z < 0.4) or (z > 5.0):
		return False
	else:	
		if (z >= 0.40) and (z < 0.68): filename = 'z.0.49'
		if (z >= 0.68) and (z < 1.06): filename = 'z.0.87'
		if (z >= 1.06) and (z < 1.44): filename = 'z.1.25'
		if (z >= 1.44) and (z < 1.82): filename = 'z.1.63'
		if (z >= 1.82) and (z < 2.21): filename = 'z.2.01'
		if (z >= 2.21) and (z < 2.60): filename = 'z.2.40'
		if (z >= 2.60) and (z < 3.03): filename = 'z.2.80'
		if (z >= 3.03) and (z < 3.50): filename = 'z.3.25'
		if (z >= 3.50) and (z < 4.00): filename = 'z.3.75'
		if (z >= 4.00) and (z < 4.50): filename = 'z.4.25'
		if (z >= 4.50) and (z <= 5.00): filename = 'z.4.75'
		
		filename = datapath+'load_sdss_dr3_lf_data.dat/'+filename
		data = np.genfromtxt(filename,names=['M_i','logphi','dphi','fill','zbar','Nqso','Ncor'],skip_header=1)
		M_i  = data['M_i']
		logphi  = data['logphi']
		dphi = (data['dphi']*1e-9)/(10**logphi)/np.log(10.)
		
		#M_B = M_i + 0.66	# Correction from the paper from Mi(z=2) to M_Bj
		M_B = M_i + 0.71	# Correction from the paper from Mi(z=2) to M_Bj
			
		L_BB = 0.4*(M_sun_Bband_AB - M_B) 
		PHI_BB = logphi + np.log10(2.5)
		DPHI_BB = dphi
		return L_BB, PHI_BB, DPHI_BB

# Function to return the analytical Richards et al. (2006) luminosity function 
#   for a list of B-band luminosities L0_list (in SOLAR luminosities) at redshift z
#   (for an Omega_M = 0.3, Omega_Lambd = 0.7 cosmology)
#
def return_sdss_dr3_lf_fitted(L0_list,z):	
	z_ref  = 2.45
	M_star = -26.0
	xsi    = np.log10((1.+z)/(1.+z_ref))
	A1     = 0.84
	A2     = -0.12
	B1     = 1.40
	B2     = 36.28
	B3     = 33.77
	log_phi= -5.72

	M_B    = M_sun_Bband_AB - 2.5*(L0_list)
	#M_i    = M_B - 0.66		# converts to their Mi(z=2) convention
	M_i    = M_B - 0.71 	# converts to their Mi(z=2) convention

	A      = A1 + A2*(z-z_ref)
	if (z > z_ref):
		A 	= A1 + A2*(z-z_ref)
	else: A = A1
	mu     = M_i - (M_star + B1*xsi + B2*xsi*xsi + B3*xsi*xsi*xsi)

	PHI    = 2.5 * (10**(log_phi)) * (10**(A*mu))
	return np.log10(PHI)
