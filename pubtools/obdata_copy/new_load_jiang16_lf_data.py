# Jiang et al. 2016
# Measuring M_1450
from data_copy import *
import numpy as np

'''#digitized version
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
'''

#official source
def load_jiang16_lf_data(z): 
	if (z < 5.7) or (z > 6.4): return False
	else:	
		z_mean = 6.05
		M_1450_left = np.array([-24.7300, -25.7400, -26.2700, -26.7800, -27.1100, -27.6100, -29.1000])
		M_1450_right= np.array([-24.8290, -25.9290, -26.4490, -26.5990, -27.1990, -27.7990, -28.6990])
		M_1450 = (M_1450_left+M_1450_right)/2.
		PHI_1450=np.log10(np.array([7.09152,  3.16253,  1.05519, 0.526870, 0.343092, 0.135817,0.0151006])*1e-9)
		sigma = np.array([2.65354,   1.27819,  0.330344,  0.173867,  0.128380, 0.0508206, 0.0149496])*1e-9
		PHI_1450_up= np.log10(10**PHI_1450+sigma)
		PHI_1450_down= np.log10(10**PHI_1450-sigma)
		DPHI_1450=(PHI_1450_up - PHI_1450_down)/2.
	return M_1450, PHI_1450, DPHI_1450