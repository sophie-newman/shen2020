# 
#  Lacy et al. 2011
#
from data import *
import numpy as np 

def load_assef11_lf_data(z): # 5 micron
	if (z <= 0) or (z > 3.8): return False
	else:	
		filename = datapath+'Lacy2015.dat'
		data = np.genfromtxt(filename,names=True)

		id = (data["z_min"]< z) & (data["z_max"]>=z) 

		L_IR = data['nuLnu'] - L_solar
		PHI_IR  = data['logphi']
		DPHI_IR = data['dlogphi']

		return L_IR, PHI_IR, DPHI_IR
