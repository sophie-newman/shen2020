import numpy as np 
from data import *
from utils import *
import scipy.interpolate as inter
from scipy.optimize import fsolve
import astropy.constants as con
import time
#all the luminosities are in log10
#all the NHs are in log10
#all Phis are "not" in log10

def bolometric_correction(L_bol,nu):  #L_bol is in log10   
	x = L_bol - 10.
	if nu==0.: return L_bol
	elif nu < 0:
		if (nu== -1.): P0, P1, P2, P3 =3.75876730, -0.36057087, 9.82954354, -6.29078046e-3
		if (nu== -2.): P0, P1, P2, P3 =4.36089518, -0.36057086, 11.4041656, -6.29077790e-3
		if (nu== -3.): P0, P1, P2, P3 =5.71209149, -0.02553868, 17.6679057,     0.27750559
		if (nu== -4.): P0, P1, P2, P3 =4.07349835, -0.02553868, 12.5996205,     0.27750559
		if (nu== -5.): P0, P1, P2, P3 =4.86961646, -0.00629077, 1.86211735,    -0.36057082
		lband = P0 * np.power(10.,P1*x) + P2*np.power(10.,P3*x)
		return L_bol - np.log10(lband) 
	#if not one of the specified bands, then take advantage of the fact that our model 
	#   spectrum is not l-dependent below 500 angstroms or above 50 angstroms, so 
	#   just take the appropriate ratios to renormalize to those luminosities
	nu_angstrom = con.c.value/1e-10
	nu500 = nu_angstrom/500.
	nu50 = nu_angstrom/50.
	if (nu <= nu500): #just take the ratio relative to B-band
		P = (3.75876730,-0.36057087,9.82954354,-6.29078046e-3)
		lband = P[0]*np.power(10.,P[1]*x) + P[2]*np.power(10.,P[3]*x)
		return np.log10(return_ratio_to_b_band(nu)) + L_bol - np.log10(lband)
	elif (nu >= nu50): #just take the ratio relative to the hard X-rays
		P = (4.07349835,-0.02553868,12.5996205,0.27750559)
		lband = P[0]*np.power(10.,P[1]*x) + P[2]*np.power(10.,P[3]*x)
		return np.log10(return_ratio_to_hard_xray(nu)) + L_bol - np.log10(lband)
	elif (nu > nu500) and (nu < nu50): #interpolate between both regimes
		P = (3.75876730,-0.36057087,9.82954354,-6.29078046e-3)
		L500 = return_ratio_to_b_band(nu500)/(P[0]*np.power(10.,P[1]*x)+P[2]*np.power(10.,P[3]*x))
		Q = (4.07349835,-0.02553868,12.5996205,0.27750559)
		L50  = return_ratio_to_hard_xray(nu50)/(Q[0]*np.power(10.,Q[1]*x)+Q[2]*np.power(10.,Q[3]*x))
		L00  = np.log10(L500) + np.log10(L50/L500) * (np.log10(nu/nu500)/np.log10(nu50/nu500))
		return L00 + L_bol

def bolometric_correction_inverse(L_obs,nu):
	def fobjective(L_bol):
		return bolometric_correction(L_bol,nu)-L_obs
	return fsolve(fobjective,x0=L_obs)[0]

