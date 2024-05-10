from data_copy import *
import numpy as np

def load_babbedge_lf_data(z): # L_IR, PHI_IR, DPHI_IR, z
	DLOGL = np.log10(8.)/5.
	mm= 1e-6
	nu_norm = np.log10(con.c.value/(15.0*mm))

	f_norm  = richards_return_spectrum(nu_norm, /BLUE)

	# 3.6 micron
	L = 9.55506 + DLOGL * np.arange(7)
	P = [2.3e-5, 2.0e-5, 1.3e-5, 7.0e-6, 4.5e-6, 3.2e-6, 1.9e-6]
	dp= [0.18,     0.17,   0.20,   0.43,   0.45,   0.43,   0.52]
	L = L - 2. * np.log10(hubble) + np.log10(f_norm/richards_return_spectrum(np.log10(con.c.value/(3.6*mm)),/BLUE))
	P = np.log10(P) + 3.*np.log10(hubble)
			
	# 4.5 micron
	L = 9.55506 + DLOGL/2. + DLOGL * np.arange(5)
	P = [2.2e-5, 1.5e-5, 9.0e-6, 6.0e-6, 3.1e-6]
	dp= [0.12,     0.33,   0.43,   0.43,   0.55]
	L = L - 2.*alog10(h) + alog10(f_norm/richards_return_spectrum(alog10(c/(4.5*mm)),/BLUE))
	P = alog10(P) + 3.*alog10(h)
	
	# 5.8 micron
	L = 9.55506 + 3.*DLOGL + DLOGL*findgen(10)
	P = [9.0d-6, 7.8d-6, 9.0d-6, 7.8d-6, 7.8d-6, 7.8d-6, 6.5d-6,  5.5d-6, 4.0d-6, 2.9d-6]
	dp= [0.46,     0.24,   0.35,   0.40,   0.38,   0.18,   0.33,    0.43,   0.43,   0.55]
	L = L - 2.*alog10(h) + alog10(f_norm/richards_return_spectrum(alog10(c/(5.8*mm)),/BLUE))
	P = alog10(P) + 3.*alog10(h)
			
	# 8.0 micron
	L = 9.55506 + 6.*DLOGL + DLOGL*findgen(5)
	P = [1.0d-5, 9.0d-6, 8.0d-6, 6.8d-6, 3.9d-6]
	dp= [0.30,     0.30,   0.60,   0.30,   0.62]
	L = L - 2.*alog10(h) + alog10(f_norm/richards_return_spectrum(alog10(c/(8.0*mm)),/BLUE))
	P = alog10(P) + 3.*alog10(h)

	# 24.0 micron
	L = 9.55506 + 4.*DLOGL + DLOGL*findgen(6)
	P = [8.0d-6, 8.5d-6, 8.0d-6, 6.5d-6, 5.5d-6, 4.0d-6]
	dp= [0.60,     0.42,   0.30,   0.40,   0.40,   0.65]
	L = L - 2.*alog10(h) + alog10(f_norm/richards_return_spectrum(alog10(c/(24.*mm)),/BLUE))
	P = alog10(P) + 3.*alog10(h)

#
# end story -- their photo-z's are *way* too inaccurate ( & their accuracy changes as a function 
#   of redshift, which could create severe bias), they ignore the hot dust when doing the 
#   photo-z's then fit it for the mid-IR excess -- their AGN template (at least for the photo-zs) 
#   looks nothing like a real AGN -- and so at the end of the day, their 
#   calculated 3.6 micron and 4.5 micron LFs are seriously discrepant from the longer wavelength LFs
#

	
	L_IR   = (DOUBLE(L_15-L_solar)) - 5.*alog10(7./7.5)
	PHI_IR = P_15 + alog10(2.5) + 3.*alog10(7./7.5)
	DPHI_IR= D_15
