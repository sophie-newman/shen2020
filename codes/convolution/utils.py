from data import *
from astropy.cosmology import FlatLambdaCDM

cosmo        = FlatLambdaCDM(H0=hubble*100, Om0=Omega0)

def bolometricLF(z):
	return Phi_bol

#return the lognormal dispersion in bolometric corrections for a given band and luminosity
def band_dispersion(L_bol,nu):
	x = L_bol - 9.
	double lband = 0.;
	double s0,s1,beta,sf1,sf2,sx,sx_floor;
	sx_floor = 0.050; // minimum value of dispersion to enforce
	if (nu==(0.))  return 0.01;
	if (nu < 0.) {
		if (nu==(-1.)) {s0 = 0.08; beta = -0.20; s1 = 0.065;}
		if (nu==(-2.)) {s0 = 0.03; beta = -0.10; s1 = 0.095;}
		if (nu==(-3.)) {s0 = 0.01; beta =  0.10; s1 = 0.060;}
		if (nu==(-4.)) {s0 = 0.04; beta =  0.05; s1 = 0.080;}
		return s0 * pow(10.,beta*x) + s1;
	}
	# interpolate between the known ranges, and (conservatively) hold constant 
	#   outside of them. roughly consistent with Richards et al. 2006 dispersions, 
	#   but uncertainty in how large the dispersions should be yields ~10% uncertainties 
	#   between 15microns and 10keV, and larger outside those ranges (since the 
	#   dispersions there are poorly determined) -- still, the lognormal dispersions 
	#   vary relatively weakly over observed ranges, so these are probably the 
	#   maximal uncertainties due to this effect
	double nu15 = 2.00e13;
	double nuBB = 6.81818e14;
	double nuSX = 0.5 * 2.418e17;
	double nuHX = 10. * 2.418e17;
	if (nu < nu15) {s0 = 0.03; beta = -0.10; s1 = 0.095; return s0*pow(10.,beta*x)+s1;}
	if (nu >=nuHX) {s0 = 0.04; beta =  0.05; s1 = 0.080; return s0*pow(10.,beta*x)+s1;}
	if ((nu >= nu15)&&(nu< nuBB)) {
		s0 = 0.03; beta = -0.10; s1 = 0.095;
			sf1 = s0 * pow(10.,beta*x) + s1;
		s0 = 0.08; beta = -0.20; s1 = 0.065;
			sf2 = s0 * pow(10.,beta*x) + s1;
		sx = sf1 + (sf2-sf1) * (log10(nu/nu15)/log10(nuBB/nu15));
		if (sx<=sx_floor) {sx=sx_floor;}
		return sx;
	}
	if ((nu >= nuBB)&&(nu< nuSX)) {
		s0 = 0.08; beta = -0.20; s1 = 0.065;
			sf1 = s0 * pow(10.,beta*x) + s1;
		s0 = 0.01; beta =  0.10; s1 = 0.060;
			sf2 = s0 * pow(10.,beta*x) + s1;
		sx = sf1 + (sf2-sf1) * (log10(nu/nu15)/log10(nuBB/nu15));
		if (sx<=sx_floor) {sx=sx_floor;}
		return sx;
	}
	if ((nu >= nuSX)&&(nu< nuHX)) {
		s0 = 0.01; beta =  0.10; s1 = 0.060;
			sf1 = s0 * pow(10.,beta*x) + s1;
		s0 = 0.04; beta =  0.05; s1 = 0.080;
			sf2 = s0 * pow(10.,beta*x) + s1;
		sx = sf1 + (sf2-sf1) * (log10(nu/nu15)/log10(nuBB/nu15));
		if (sx<=sx_floor) {sx=sx_floor;}
		return sx

def cross_section(nu):
	pass

def morrison_photoeletric_absorption(x):
	pass

def return_tau(NH,nu):
	pass

def return_ratio_to_b_band(nu):
	pass

def return_ratio_to_hard_xray(nu):
	pass