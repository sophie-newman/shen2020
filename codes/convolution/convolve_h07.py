import numpy as np 
from data import *
import astropy.constants as con
#all the luminosities are in log10
#all the NHs are in log10
#all Phis are "not" in log10

#bolometric correction function of HRH07

def bolometric_correction_old(L_bol,nu):  #L_bol is in log10   
	x = L_bol - 10.
	if nu==0.: return L_bol
	elif nu < 0:
		if (nu== -1.): P0, P1, P2, P3 =8.99833, 6.24800, -0.370587, -0.0115970
		if (nu== -2.): P0, P1, P2, P3 =10.6615, 7.40282, -0.370587, -0.0115970
		if (nu== -3.): P0, P1, P2, P3 =10.0287, 17.8653, 0.276804,  -0.0199558
		if (nu== -4.): P0, P1, P2, P3 =6.08087, 10.8331, 0.276802,  -0.0199597
		lband = P0 * np.power(10.,P3*x) + P1*np.power(10.,P2*x)
		return L_bol - np.log10(lband) 
