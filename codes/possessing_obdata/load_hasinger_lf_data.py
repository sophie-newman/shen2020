# 
# HMS 2005 LF data (supplants Miyaji LF data)
#
from data import *
import numpy as np

def load_hasinger_lf_data(z): # L_SX, PHI_SX, DPHI_SX, z
	# determine which redshift interval its in
	if ((z >= 0.0) and (z < 0.2)): WHICH_BLOCK = 1
	if ((z >= 0.2) and (z < 0.4)): WHICH_BLOCK = 2
	if ((z >= 0.4) and (z < 0.8)): WHICH_BLOCK = 3
	if ((z >= 0.8) and (z < 1.6)): WHICH_BLOCK = 4
	if ((z >= 1.6) and (z < 3.2)): WHICH_BLOCK = 5
	if ((z >= 3.2) and (z <=4.8)): WHICH_BLOCK = 6
	if (z > 4.8): WHICH_BLOCK = 7

	if (WHICH_BLOCK == 7): return False
	else:
		BEGIN = False
		filename=datapath+'load_hasinger_lf_data.dat'
		with open(filename, 'r') as f:
			i=0
			j=0
			for line in f.readlines():
				elements=line.split()
				if (BEGIN == True) and (j<N_LUM_BINS):
					L_SX_min[j] = float(elements[2])
					L_SX_max[j] = float(elements[3])
					PHI_SX[j] = float(elements[4])
					DPHIp_SX[j] = float(elements[5]) 
					DPHIm_SX[j] = float(elements[6])
					j+=1
		 		if (len(elements[0]) == 1) and ( int(elements[0])==WHICH_BLOCK ):
					BEGIN = True
					N_LUM_BINS = int(elements[1])
					L_SX_min = np.zeros(N_LUM_BINS)
					L_SX_max = np.zeros(N_LUM_BINS)
					PHI_SX   = np.zeros(N_LUM_BINS)
					DPHIp_SX = np.zeros(N_LUM_BINS)
					DPHIm_SX = np.zeros(N_LUM_BINS)
					N_QSO_SX   = np.zeros(N_LUM_BINS,dtype=np.int32)
				i+=1

		L_SX = 0.5 * (L_SX_min + L_SX_max)

		ok = (L_SX < 48.0)# & (L_SX > 42.)		# things go seriously wonky with the Ueda fits if these are included
		L_SX = L_SX[ok]
		PHI_SX = PHI_SX[ok]
		DPHIp_SX = DPHIp_SX[ok]
		DPHIm_SX = DPHIm_SX[ok]

		# want to correct for the type-2 fraction, since this is only a Type-1 
		#  luminosity function: do so following the data in Hasinger et al. 2004
		#  as compiled in Simpson et al. 2005
		#
		L0 = 10**42.37/0.015 * 0.606	# our spectrum HX->SX
		xi = 0.23
		f1 = 1. - ((1.+3.*((10**L_SX/L0)**(1.-2.*xi)))**(-0.5))
		f1 = 1./(1.+ 10**(-(L_SX - 45.0) * 0.4)) * 3.0
		PHI_SX = PHI_SX + np.log10(1./f1) 

		DPHI_SX = 0.5 * (DPHIp_SX + np.abs(DPHIm_SX))
		L_SX = L_SX-L_solar
		return L_SX, PHI_SX, DPHI_SX
