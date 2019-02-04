import numpy as np
from data import *

def LF(L_bol,P):
	gamma1	= P[0]	#faint-end slope
	gamma2	= P[1] 	#bright-end slope
	P0		= P[2]	#normalization
	Lbreak	= P[3]	#luminosity break
	x = np.power(10.,(L_bol-Lbreak))
	P_temp = np.power(10.,P0) / ( np.power(x,gamma1) + np.power(x,gamma2) )
	P_temp = np.log10(P_temp)		    
	return P_temp

def LF_at_z(L_bol,P,z,model):
	if model=='Fiducial':
		xsi_log	= np.log10((1.+z)/(1.+2))
		xsi_lin	= z - 2
		xsi = xsi_lin
		xsi = xsi_log
	
		gamma1_0	= P[0]	#faint-end slope
		gamma2_0	= P[1] 	#bright-end slope
		P0			= P[2]	#normalization in log
		L0			= P[3]	#break in log
		k1, k2, k3 = P[4], P[5], P[6]
		k_gamma1 = P[7]
		k_gamma2_1 = P[8]
		k_gamma2_2 = P[9]

		gamma1   = gamma1_0 * np.power(10., k_gamma1*xsi_log)
		gamma2   = 2.*gamma2_0 / (np.power(10., xsi_log*k_gamma2_1) + np.power(10., xsi_log*k_gamma2_2))
		BETA_MIN = 1.3
		if KEYS["LF_model"]=="Schechter": BETA_MIN = 0.02
		if gamma2<BETA_MIN: gamma2 = BETA_MIN

		Lbreak  = L0 + k1*xsi + k2*xsi**2 + k3*xsi**3

		P_temp = LF( L_bol, [gamma1,gamma2,P0,Lbreak])
		return P_temp
'''	
T0 = np.polynomial.chebyshev.Chebyshev((1,0,0,0))
T1 = np.polynomial.chebyshev.Chebyshev((0,1,0,0))
T2 = np.polynomial.chebyshev.Chebyshev((0,0,1,0))
T3 = np.polynomial.chebyshev.Chebyshev((0,0,0,1))
def doublepower(z,p):
        xsi=1.+z
        zref=p[1]
        return 2*p[0]/(np.power(xsi/(1+zref),p[2]) + np.power(xsi/(1+zref),p[3]))

def LF_at_z(L_bol,P,z,model):
        if model=='Fiducial':
                xsi = 1.+z
		gamma1=P[0]*T0(xsi)+P[1]*T1(xsi)+P[2]*T2(xsi)+P[3]*T3(xsi)
		gamma2=doublepower(z,[P[4],P[5],P[6],P[7]])
		Phis  =P[8]*T0(xsi)+P[9]*T1(xsi)
		Lbreak=doublepower(z,[P[10],P[11],P[12],P[13]])

                P_temp = LF( L_bol, [gamma1,gamma2,Phis,Lbreak])
                return P_temp
'''
