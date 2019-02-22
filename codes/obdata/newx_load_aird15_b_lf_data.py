# 
# Aird et al. 2015b NuSTAR
#
from data import *
import numpy as np 

def load_aird15_b_lf_data(z): # L_HX, PHI_HX, DPHI_HX
	if (z <= 0.1) or (z > 3.): return False
	elif (z > 0.1) and (z <= 0.5):	
		L_HX_lo = np.array([42.75 ,43.00 ,43.25 ,43.50 ,43.75 ,44.00 ,44.25])
		L_HX_hi = np.array([43.00 ,43.25 ,43.50 ,43.75 ,44.00 ,44.25 ,44.50])
		
		P_HX = np.array([1.81e-4 ,9.39e-5 ,5.69e-5 ,2.63e-5 ,8.48e-6 ,3.20e-6 ,1.81e-6])
		sigma_up  = np.array([2.14e-4 ,8.13e-5 ,3.50e-5 ,1.46e-5 ,6.33e-6 ,2.95e-6 ,1.71e-6])
		sigma_down= np.array([1.08e-4 ,4.67e-5 ,2.27e-5 ,0.98e-5 ,3.85e-6 ,1.66e-6 ,0.95e-6])
	elif (z > 0.5) and (z <= 1.0):
		L_HX_lo = np.array([43.50 ,43.75 ,44.00 ,44.25 ,44.50 ,44.75])
		L_HX_hi = np.array([43.75 ,44.00 ,44.25 ,44.50 ,44.75 ,45.00])
		
		P_HX = np.array([1.74e-4, 5.92e-5, 2.32e-5, 1.53e-5, 4.17e-6, 3.74e-7])
		sigma_up  = np.array([3.09e-4, 5.42e-5, 1.50e-5, 0.57e-5, 1.98e-6, 6.38e-7])
		sigma_down= np.array([1.27e-4, 3.05e-5, 0.96e-5, 0.42e-5, 1.39e-6, 2.68e-7])
	elif (z > 1.0) and (z <= 3.0):
		L_HX_lo = np.array([44.25 , 44.50 , 44.75 , 45.00 , 45.25 , 45.50])
		L_HX_hi = np.array([44.50 , 44.75 , 45.00 , 45.25 , 45.50 , 45.75])
		
		P_HX = np.array([4.21e-5, 1.61e-5, 5.69e-6, 2.00e-6, 3.76e-7, 1.56e-7])
		sigma_up  = np.array([4.73e-5, 1.16e-5, 3.09e-6, 1.00e-6, 3.19e-7, 1.76e-7])
		sigma_down= np.array([2.44e-5, 0.71e-5, 2.09e-6, 0.69e-6, 1.85e-7, 0.91e-7])

	L_HX = (L_HX_hi+L_HX_lo)/2.
	P_HX_up = np.log10( P_HX + sigma_up  )
	P_HX_down = np.log10( P_HX - sigma_down )
	D_HX = (P_HX_up - P_HX_down)/2.

	PHI_HX = np.log10(P_HX)
	DPHI_HX= D_HX
	L_HX = (L_HX - L_solar)

	return L_HX, PHI_HX, DPHI_HX
