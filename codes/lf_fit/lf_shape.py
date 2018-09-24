import numpy as np
from data import *

def LF(L_bol,P):
	P=float(P)
	alpha	= P[0]	#faint-end slope
	beta	= P[1] 	#bright-end slope
	P0		= P[2]	#normalization
	L0		= P[3]	#luminosity break
	x = np.power(10.,(L_bol-L0))
	P_temp = np.power(10.,P0) / ( np.power(x,alpha) + np.power(x,beta) )
	P_temp = np.log10(P_temp)		    
	return P_temp

def LF_at_z(L_bol,P,z,model):
	P=float(P)
	if model=='LDDE':
		x = np.power(10., (L_bol + np.log10(3.9) + 33. - P[1]))
		
		phi_0 = np.power(10., P[0] /(x**P[2] + x**P[3]))	#z=0 LF

		p1 = P[4] + P[7] * (L_bol + np.log10(3.9) +33. - 46.)
		p2 = P[5] + P[8] * (L_bol + np.log10(3.9) +33. - 46.)
		xc = np.power(10., (L_bol+ np.log10(3.9) + 33. - P[9]))
		
		zc = 0.0 * xc
		zc[xc<=1.] = P[6]*(xc**P[10])
		zc[xc>1.] = P[6]

		ed = 0.0 * zc
		ed[z<=zc] = np.power(1.+z,p1)
		ed[z>zc] = np.power(1.+zc,p1-p2) * np.power(1.+z,p2)
		
		P_temp = np.log10(phi_0*ed)
	
	elif model=='Fiducial':
		xsi_log	= np.log10((1.+z)/(1.+Z_NORM_SET))
		xsi_lin	= z - Z_NORM_SET
		xsi = xsi_lin
		if (FIT_LOG_Z) then xsi = xsi_log
	
		alpha	= P[0]	#faint-end slope
		beta	= P[1] 	#bright-end slope
		P0		= P[2]	#normalization
		L0		= P[3]	#break
	
		alpha   = P[0] * np.power(10., P[6]*xsi_log + P[7]*xsi_log**2)
		beta    = 2.*P[1] / np.power(10., xsi_log*P[8] + np.power(10., xsi_log*P[9]))
		BETA_MIN = 1.3
		if KEYS["LF_model"]=="Schechter": BETA_MIN = 0.02
		if beta<BETA_MIN: beta = BETA_MIN
		L0  = P[3] + P[4]*xsi + P[5]*xsi**2 + P[11]*xsi**3 + P[12]*xsi**4
		if P[13] > 0.: P0 = P[2] - np.log10(1. + np.power( (1.+z)/(1.+P[13]),P[14]))
		if P[13] > 0.:
			if (z <= P[13]): P0 = P[2]
			if (z > P[13]): P0 = P[2] + P[14]*np.log10((1.+z)/(1.+P[13]))
	
		x = np.power(10.,(log_L-L0))
		P_temp = LF( np.log10(x),[alpha,beta,P0,0.0])
		return, P_temp