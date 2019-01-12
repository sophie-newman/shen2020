# Jiang et al. 2016
# Measuring M_1450
from data import *
import numpy as np

def load_jiang16_lf_data(z): 
	if (z < 5.7) or (z > 6.4): return False
	else:	
		z_mean = 6.05
		M_1450=np.array([-24.83054,-25.92282,-26.43960,-26.59228,-27.19128,-27.79027,-28.69463])
		PHI_1450=np.array([-8.14801,-8.50295,-8.97014,-9.27017,-9.45805,-9.86234,-10.81781])
		PHI_1450_up=np.array([-8.01538,-8.34937,-8.85845,-9.15847,-9.33240,-9.73669])
		PHI_1450_down=np.array([-8.33650,-8.71237,-9.13070,-9.43771,-9.65351,-10.05779])
		DPHI_1450=(PHI_1450_up - PHI_1450_down)/2.

	return M_1450, PHI_1450, DPHI_1450
