from data import *
import numpy as np 
import astropy.constants as con

def load_miyaji_lf_data(z): # L_SX, PHI_SX, DPHI_SX, z
	# determine which redshift interval its in
	if ((z >= 0.0) and (z < 0.2)): WHICH_BLOCK = 1
	if ((z >= 0.2) and (z < 0.4)): WHICH_BLOCK = 2
	if ((z >= 0.4) and (z < 0.8)): WHICH_BLOCK = 3
	if ((z >= 0.8) and (z < 1.6)): WHICH_BLOCK = 4
	if ((z >= 1.6) and (z < 2.3)): WHICH_BLOCK = 5
	if ((z >= 2.3) and (z <=4.6)): WHICH_BLOCK = 6
	if (z > 4.6): WHICH_BLOCK = 7

	if (WHICH_BLOCK == 7): return False
	else:
		BEGIN = False
		filename=datapath+'load_miyaji_lf_data.dat'
		with open(filename, 'r') as f:
			i=0
			j=0
			for line in f.readlines():
				elements=line.split()
				if (BEGIN == True) and (j<N_LUM_BINS):
					L_SX_min[j] = float(elements[2])
					L_SX_max[j] = float(elements[3])
					PHI_SX[j]   = float(elements[4])
					DPHIp_SX[j] = float(elements[5]) 
					DPHIm_SX[j] = float(elements[6])
					N_QSO_SX[j] = int(elements[7])
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
		L_SX = 10**(L_SX - np.log10(con.L_sun.value*1e7)) * (0.7/0.5)**(-2.)
		DPHI_SX = 0.5 * (np.abs(DPHIp_SX) + np.abs(DPHIm_SX))
		ok = ( DPHI_SX > 0.)
		L_SX = L_SX[ok]
		DPHI_SX = DPHI_SX[ok]
		PHI_SX = PHI_SX[ok] + 3.*np.log10(7./5.)
		ok = (np.log10(L_SX) < 47.0) # things go seriously wonky with the Ueda fits if these are included
		L_SX = np.log10(L_SX[ok])
		PHI_SX = PHI_SX[ok]
		DPHI_SX = DPHI_SX[ok]
		return L_SX, PHI_SX, DPHI_SX
