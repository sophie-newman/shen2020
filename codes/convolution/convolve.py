import numpy as np 
from data import *
from utils import *
import scipy.interpolate as inter
import astropy.constants as con
import time
#all the luminosities are in log10
#all the NHs are in log10
#all Phis are "not" in log10

def bolometric_correction(L_bol,nu):  #L_bol is in log10   
	x = L_bol - 10.
	if nu==0.: return L_bol
	elif nu < 0:
		if (nu== -1.): P0, P1, P2, P3 =8.99833, 6.24800, -0.370587, -0.0115970
		if (nu== -2.): P0, P1, P2, P3 =10.6615, 7.40282, -0.370587, -0.0115970
		if (nu== -3.): P0, P1, P2, P3 =10.0287, 17.8653, 0.276804,  -0.0199558
		if (nu== -4.): P0, P1, P2, P3 =6.08087, 10.8331, 0.276802,  -0.0199597
		lband = P0 * np.power(10.,P3*x) + P1*np.power(10.,P2*x)
		return L_bol - np.log10(lband) 
	#if not one of the specified bands, then take advantage of the fact that our model 
	#   spectrum is not l-dependent below 500 angstroms or above 50 angstroms, so 
	#   just take the appropriate ratios to renormalize to those luminosities
	nu_angstrom = con.c.value/1e-10
	nu500 = nu_angstrom/500.
	nu50 = nu_angstrom/50.
	if (nu <= nu500): #just take the ratio relative to B-band
		P = (8.99833, 6.24800, -0.370587, -0.0115970)
		lband = P[0]*np.power(10.,P[3]*x) + P[1]*np.power(10.,P[2]*x)
		return np.log10(return_ratio_to_b_band(nu)) + L_bol - np.log10(lband)
	elif (nu >= nu50): #just take the ratio relative to the hard X-rays
		P = (6.08087, 10.8331, 0.276802, -0.0199597)
		lband = P[0]*np.power(10.,P[3]*x) + P[1]*np.power(10.,P[2]*x)
		return np.log10(return_ratio_to_hard_xray(nu)) + L_bol - np.log10(lband)
	elif (nu > nu500) and (nu < nu50): #interpolate between both regimes
		P = (8.99833, 6.24800, -0.370587, -0.0115970)
		L500 = return_ratio_to_b_band(nu500)/(P[0]*np.power(10.,P[3]*x)+P[1]*np.power(10.,P[2]*x))
		Q = (6.08087, 10.8331, 0.276802 , -0.0199597)
		L50  = return_ratio_to_hard_xray(nu50)/(Q[0]*np.power(10.,Q[3]*x)+Q[1]*np.power(10.,Q[2]*x))
		L00  = np.log10(L500) + np.log10(L50/L500) * (np.log10(nu/nu500)/np.log10(nu50/nu500))
		return L00 + L_bol

#return the appropriate jacobian factors for the above (dlogL/dlogL_band)
def bolometric_correction_jacobian(L_bol,nu):
	x = L_bol - 10.

	if nu == 0.:  return 1.0
	elif nu < 0.:
		if nu == -1.: P0, P1, P2, P3 = 8.99833, 6.24800, -0.370587, -0.0115970
		if nu == -2.: P0, P1, P2, P3 = 10.6615, 7.40282, -0.370587, -0.0115970
		if nu == -3.: P0, P1, P2, P3 = 10.0287, 17.8653, 0.276804 , -0.0199558
		if nu == -4.: P0, P1, P2, P3 = 6.08087, 10.8331, 0.276802 , -0.0199597
		D1 = P0*(1.+P3)*np.power(10.,P3*x) + P1*(1.+P2)*np.power(10.,P2*x);
		D2 = P0*np.power(10.,P3*x) + P1*np.power(10.,P2*x);
		return D1/D2

	nu_angstrom = con.c.value/1e-10
	nu500 = nu_angstrom/500.
	nu50 = nu_angstrom/50.
	if (nu <= nu500): #just take the ratio relative to B-band
		P = (8.99833, 6.24800, -0.370587, -0.0115970)
		D1 = P[0]*(1.+P[3])*np.power(10.,P[3]*x) + P[1]*(1.+P[2])*np.power(10.,P[2]*x);
		D2 = P[0]*np.power(10.,P[3]*x) + P[1]*np.power(10.,P[2]*x);
		return D1/D2
	elif (nu >= nu50): #just take the ratio relative to the hard X-rays
		P = (6.08087, 10.8331, 0.276802, -0.0199597)
		D1 = P[0]*(1.+P[3])*np.power(10.,P[3]*x) + P[1]*(1.+P[2])*np.power(10.,P[2]*x);
		D2 = P[0]*np.power(10.,P[3]*x) + P[1]*np.power(10.,P[2]*x);
		return D1/D2
	elif (nu > nu500) and (nu < nu50): #interpolate between both regimes
		P = (8.99833, 6.24800, -0.370587, -0.0115970)
		D1 = P[0]*(1.+P[3])*np.power(10.,P[3]*x) + P[1]*(1.+P[2])*np.power(10.,P[2]*x);
		D2 = P[0]*np.power(10.,P[3]*x) + P[1]*np.power(10.,P[2]*x);
		L500=D1/D2
		Q = (6.08087, 10.8331, 0.276802 , -0.0199597)
		D1 = Q[0]*(1.+Q[3])*np.power(10.,Q[3]*x) + Q[1]*(1.+Q[2])*np.power(10.,Q[2]*x);
		D2 = Q[0]*np.power(10.,Q[3]*x) + Q[1]*np.power(10.,Q[2]*x);
		L50=D1/D2
		L00  = np.log10(L500) + np.log10(L50/L500) * (np.log10(nu/nu500)/np.log10(nu50/nu500));	
		return np.power(10.,L00)

#######################################################################

def uncertainty_correction(Phi_bol,nu):
	L_band = bolometric_correction(L_bol_grid,nu)
	sig=0.0*L_bol_grid
	for i in range(len(L_bol_grid)):
		sig[i] = band_dispersion(L_bol_grid[i],nu)
	
	expfac = -0.5/sig**2
	prefac = 1./(sig*np.sqrt(2*np.pi))*Phi_bol*d_log_l_bol

	Phi_band = 0.0*L_band
	for i in range(len(L_bol_grid)):
		Phi_band[i] += np.sum( prefac*np.exp(expfac*(L_band-L_band[i])**2) )
	return Phi_band

def extinction_correction(Phi_band, nu):
	NHs = np.linspace(20, 25, 501)
	dNH = NHs[1]-NHs[0]
	
	taus=0.0*NHs
	for i in range(len(NHs)):
		taus[i] = return_tau(NHs[i], nu)

	eps = 1.7
	psi_44 = 0.47
	beta_L = 0.10
	psi_max = (1.+eps)/(3.+eps)

	L_band = bolometric_correction(L_bol_grid,nu)
	Phi_obs_corrected = Phi_band * 0.

	for i in range(len(NHs)):
		L_obs = L_band - taus[i]/np.log(10.)
		Phi_obs = Phi_band

		for j in range(len(L_band)):
			if L_band[j] <= L_obs[-1]: 
				p0 = inter.interp1d(L_obs, np.log10(Phi_obs))(L_band[j])
			else:
				k = (np.log10(Phi_obs[-1])-np.log10(Phi_obs[-2]))/(L_obs[-1]-L_obs[-2])
				p0 = np.log10(Phi_obs[-2]) + k*(L_band[j]-L_obs[-2])
			L_HX = bolometric_correction(L_bol_grid[j],-4)
			psi = psi_44 - beta_L * (L_HX + L_solar - 44.0)
			if psi < 0: psi = 0
			if psi > psi_max: psi = psi_max

			f_low = 2.0 - ((5.+2.*eps)/(1.+eps))*psi
			f_med = (1./(1.+eps))*psi
			f_hig = (eps/(1.+eps))*psi
			f_compton = f_hig
			f_low = f_low / (1. + f_compton)
			f_med = f_med / (1. + f_compton)
			f_hig = f_hig / (1. + f_compton)

			if NHs[i] <= 20.5: f_NH = f_low
			elif (NHs[i] > 20.5) and (NHs[i] <= 23.0): f_NH = f_med
			elif (NHs[i] > 23.0) and (NHs[i] <= 24.0): f_NH = f_hig
			elif NHs[i] > 24.0: f_NH = f_compton
			dN_NH = f_NH * dNH
			Phi_obs_corrected[j] += np.power(10., p0) * dN_NH

	return Phi_obs_corrected

def convolve(Phi_bol,nu):
	l_band=bolometric_correction(L_bol_grid,nu)
	P_1=Phi_bol * bolometric_correction_jacobian(L_bol_grid,nu)
	P_2=uncertainty_correction(P_1,nu)
	phi_obs= extinction_correction(P_2,nu)
	return l_band, phi_obs

def printall(redshift,nu):
	nu_eff=nu
	if nu ==  0.: nu_eff = 1.1992000e15 # 2500 angstrom
	if nu == -1.: nu_eff = 6.8136364e14 # 4400 angstrom
	if nu == -2.: nu_eff = 1.9986667e13 # 15 micron
	if nu == -3.: nu_eff = 2.4180000e17 # 'effective' 1 keV
	if nu == -4.: nu_eff = 1.2090000e18 # 'effective' 5 keV

	phi_bol= bolometricLF(L_bol_grid,redshift)
	l_band, phi_obs=convolve(phi_bol,nu)

	m_AB_band=-2.5*l_band+2.5*np.log10(nu_eff)-32.38265724887536
	m_AB_obs = m_AB_band - distance_modulus(redshift)
	S_nu = -0.4*(m_AB_obs-16.40) # returns log_{10}(S_nu/mJy)
	if (nu == 0.) or (nu == -3.) or (nu == -4.): 
		S_nu = l_band + 0.4*distance_modulus(redshift) - 6.486937100449856 
		# last factor for CGS conversion, 
		#  given that l_band is in L_sun
	S_nu = np.power(10.,S_nu)

	print "L band:",l_band + L_solar
	print "M_AB:  ",m_AB_band
	print "S_nu:  ",S_nu
	print "L_bol: ",L_bol_grid + L_solar
	print "Phiobs:",phi_obs
