import numpy as np
from data_copy import *

def LF(L_bol,P):
	gamma1	= P[0]	#faint-end slope
	gamma2	= P[1] 	#bright-end slope
	P0	= P[2]	#normalization
	Lbreak	= P[3]	#luminosity break
	x = np.power(10.,(L_bol-Lbreak))
	P_temp = np.power(10.,P0) / ( np.power(x,gamma1) + np.power(x,gamma2) )
	P_temp = np.log10(P_temp)		    
	return P_temp

def LF_at_z_H07(L_bol,P,z,model):
	if model=='Fiducial':
		xsi_log	= np.log10((1.+z)/(1.+2))
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
	
T0 = np.polynomial.chebyshev.Chebyshev((1,0,0,0))
T1 = np.polynomial.chebyshev.Chebyshev((0,1,0,0))
T2 = np.polynomial.chebyshev.Chebyshev((0,0,1,0))
T3 = np.polynomial.chebyshev.Chebyshev((0,0,0,1))
def polynomial(z,p,n=3): #Chebyshev polynomials
        xsi=1.+z
        if n==1: return p[0]*T0(xsi)+p[1]*T1(xsi)
        elif n==2: return p[0]*T0(xsi)+p[1]*T1(xsi)+p[2]*T2(xsi)
        elif n==3:
                return p[0]*T0(xsi)+p[1]*T1(xsi)+p[2]*T2(xsi)+p[3]*T3(xsi)
        else: return False

def doublepower(z,p): #double power law, definition slightly different from the DPL LF
        xsi=1.+z
        zref=p[1]
        return 2*p[0]/(np.power(xsi/(1+zref),p[2]) + np.power(xsi/(1+zref),p[3]))

def powerlaw_gamma1(z,p): #powerlaw, defined here for the evolution of the faint end slope
	xsi=1.+z
	zref=p[1]
	return p[0] * np.power(xsi/(1+zref),p[2])

def LF_at_z(L_bol,P,z,model):
	zref = 2.
        if model=='Fiducial':
                xsi = 1.+z
		gamma1=P[0]*T0(xsi)+P[1]*T1(xsi)+P[2]*T2(xsi)
		gamma2=doublepower(z,[P[3],zref,P[4],P[5]])
		Phis  =P[6]*T0(xsi)+P[7]*T1(xsi)
		Lbreak=doublepower(z,[P[8],zref,P[9],P[10]])

                P_temp = LF( L_bol, [gamma1,gamma2,Phis,Lbreak])
		if len(P_temp[np.invert(np.isfinite(P_temp))])!=0:
			return np.ones(len(L_bol))*(-40.)
		else: 
			return P_temp
	elif model=='ShallowFaint': #models that has shallow faint end at high z
		xsi = 1.+z
		gamma1=powerlaw_gamma1(z,(P[0],zref,P[1]))
                gamma2=doublepower(z,[P[2],zref,P[3],P[4]])
                Phis  =P[5]*T0(xsi)+P[6]*T1(xsi)
                Lbreak=doublepower(z,[P[7],zref,P[8],P[9]])

                P_temp = LF( L_bol, [gamma1,gamma2,Phis,Lbreak])
                if len(P_temp[np.invert(np.isfinite(P_temp))])!=0:
                        return np.ones(len(L_bol))*(-40.)
                else:
                        return P_temp

def pars_at_z(fit_evolve,redshift,model="Fiducial"):
	zref = 2.
	p=fit_evolve['gamma1']
	if model=="Fiducial": gamma1=polynomial(redshift,p)
	elif model=="ShallowFaint": gamma1=powerlaw_gamma1(redshift,(p[0],zref,p[1]))
	p=fit_evolve['gamma2']
	gamma2=doublepower(redshift,(p[0],zref, p[1], p[2]))
	p=fit_evolve['phis']
	logphis=polynomial(redshift,p)
	p=fit_evolve['Lbreak']
	Lbreak=doublepower(redshift,(p[0],zref, p[1], p[2]))
	parameters=np.array([gamma1,gamma2,logphis,Lbreak])
	return parameters

