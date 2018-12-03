# Yang et al. 2018
# Measuring M_i, converted to M_1450
from data import *
import numpy as np

def load_yang18_lf_data(z): 
	if (z < 0.5) or (z > 4.5): return False
	elif (z>=0.5) and (z<1.0):
		z_c = 0.75
		M_1450=np.array([-23.5, -22.5, -21.5, -20.5])
		M_mean=np.array([-23.410, -22.588, -21.363, -20.422])
		PHI_1450=np.array([-6.174, -6.169, -5.662, -5.457])
		sigma=np.array([0.474, 0.479, 0.902, 1.935])*1e-6	
	elif (z>=1.0) and (z<1.5):
		z_c = 1.25
		M_1450=np.array([-25.823, -24.5, -23.5, -22.5, -21.5])
		M_mean=np.array([-25.823, -24.212, -23.658, -22.539, -21.566])
		PHI_1450=np.array([-6.687, -6.393, -5.907, -5.889, -5.481])
		sigma=np.array([0.206, 0.286, 0.506, 0.527, 1.331])*1e-6
	elif (z>=1.5) and (z<2.0):
		z_c = 1.75
		M_1450=np.array([-25.485, -23.5, -22.5, -21.5])
		M_mean=np.array([-25.485, -23.491, -22.627, -21.904])
		PHI_1450=np.array([-6.772, -5.859, -5.948, -5.339])
		sigma=np.array([0.169, 0.533, 0.518, 2.629])*1e-6
	elif (z>=2.0) and (z<2.5):
		z_c = 2.25
		M_1450=np.array([-25.5, -24.5, -23.5, -22.5])
		M_mean=np.array([-25.436, -24.202, -23.511, -22.558])
		PHI_1450=np.array([-6.194, -6.094, -5.652, -5.323])
		sigma=np.array([0.320, 0.361, 0.687, 1.766])*1e-6
	elif (z>=2.5) and (z<3.0):
		z_c = 2.75
		M_1450=np.array([-26.765, -25.5, -24.5, -23.5])
		M_mean=np.array([-26.765, -25.669, -24.813, -23.69])
		PHI_1450=np.array([-6.792, -6.501, -6.191, -6.061])
		sigma=np.array([0.162, 0.223, 0.322, 0.444])*1e-6
	elif (z>=3.0) and (z<3.5):
		z_c = 3.25	
		M_1450=np.array([-26.5, -24.5, -23.5])
		M_mean=np.array([-26.483, -24.457, -23.607])
		PHI_1450=np.array([-6.478, -6.156, -5.98])
		sigma=np.array([0.235, 0.349, 0.667])*1e-6
	elif (z>=3.5) and (z<4.0):	
		z_c = 3.75
		M_1450=np.array([-26.998, -24.835])
		M_mean=np.array([-26.998, -24.835])
		PHI_1450=np.array([-6.771, -6.739])
		sigma=np.array([0.170, 0.183])*1e-6	
	elif (z>=4.0) and (z<=4.5):	
		z_c = 4.25
		M_1450=np.array([-26.342, -25.050])
		M_mean=np.array([-26.342, -25.050])
		PHI_1450=np.array([-6.738, -6.532])
		sigma=np.array([0.183, 0.294])*1e-6	

	M_1450 = M_1450 - 2.5*np.log10(lum_correct_cosmo_flexible(z_c, 0.7, 0.272))
	PHI_1450 = PHI_1450 + np.log10(phi_correct_cosmo_flexible(z_c, 0.7, 0.272))
	sigma = sigma * phi_correct_cosmo_flexible(z_c, 0.7, 0.272)

	DPHI_1450= ( (np.log10(10**PHI_1450+sigma)-PHI_1450)+(PHI_1450-np.log10(10**PHI_1450-sigma)) )/2.
	return M_1450, PHI_1450, DPHI_1450
