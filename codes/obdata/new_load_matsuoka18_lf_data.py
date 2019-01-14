# Matsuoka et al. 2018
# Measuring M_1450
from data import *
import numpy as np

def load_matsuoka18_lf_data(z): 
	if (z < 5.7) or (z > 6.5): return False
	else:	
		z_c=6.1
		M_1450=np.array([-22.75, -23.25, -23.75, -24.25, -24.75, -25.25, -25.75, -26.25, -26.75, -27.50])
		PHI_1450= np.log10( np.array([23.0, 10.9, 8.3, 6.6, 7.0, 4.6, 1.33, 0.90, 0.58, 0.242])*1e-9 )
		sigma=np.array([8.1, 3.6, 2.6, 2.0, 1.7, 1.2, 0.60, 0.32, 0.17, 0.061])*1e-9

		DPHI_1450= ( (np.log10(10**PHI_1450+sigma)-PHI_1450)+(PHI_1450-np.log10(10**PHI_1450-sigma)) )/2.

		return M_1450, PHI_1450, DPHI_1450

#bins -22 and -29 are removed from original dataset, since there is only one object in the bin and the uncertainty is infinite
