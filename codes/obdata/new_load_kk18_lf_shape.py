# Kulkarni et al. 2018
from data import *
import numpy as np

T0 = np.polynomial.chebyshev.Chebyshev((1,0,0,0))
T1 = np.polynomial.chebyshev.Chebyshev((0,1,0,0))
T2 = np.polynomial.chebyshev.Chebyshev((0,0,1,0))
T3 = np.polynomial.chebyshev.Chebyshev((0,0,0,1))

# Function to return the analytical Kulkarni et al. (2006) luminosity function 
#   for a list of M1450 at redshift z
#   (for an Omega_M = 0.3, Omega_Lambd = 0.7 cosmology)
#
def return_kk18_lf_fitted(M_list,z):
	#model 1	
	#c0=np.array([-7.559 , 1.013, -0.113])
	#c1=np.array([-17.006, 5.548, 0.588, -0.023])
	#c2=np.array([-3.246 , -0.25])
	#c3=np.array([-2.350 , 0.647, 3.857, 27.534, -0.002])
	#model 2
	c0=np.array([-7.084 , 0.753, -0.096])
	c1=np.array([-15.423, -6.725, 0.737, -0.029])
	c2=np.array([-2.973 , -0.347])
	c3=np.array([-2.545 , 1.581, 2.102, 1.965, -0.641])

	phi_s=10**( c0[0]*T0(1.+z)+c0[1]*T1(1.+z)+c0[2]*T2(1.+z))

	M_s  =c1[0]*T0(1.+z)+c1[1]*T1(1.+z)+c1[2]*T2(1.+z)+c1[3]*T3(1+z)

	alpha=c2[0]*T0(1.+z)+c2[1]*T1(1.+z)

	ksi = np.log10((1.+z)/(1.+c3[2]))
	beta = c3[0] + c3[1]/(10**(c3[3]*ksi)+10**(c3[4]*ksi))

	PHI=phi_s/( 10**(0.4*(alpha+1)*(M_list-M_s)) + 10**(0.4*(beta+1)*(M_list-M_s)) )
	return np.log10(PHI)

