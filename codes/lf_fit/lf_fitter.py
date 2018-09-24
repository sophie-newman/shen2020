from data import *
import numpy as np 
from lf_shape import *
import scipy.interpolate as inter

def load_sdss_dr3_LF_data(z):
	if (z >= 0.40) and (z <= 0.68): filename = 'z.0.49'
	if (z >= 0.68) and (z <= 1.06): filename = 'z.0.87'
	if (z >= 1.06) and (z <= 1.44): filename = 'z.1.25'
	if (z >= 1.44) and (z <= 1.82): filename = 'z.1.63'
	if (z >= 1.82) and (z <= 2.21): filename = 'z.2.01'
	if (z >= 2.21) and (z <= 2.60): filename = 'z.2.40'
	if (z >= 2.60) and (z <= 3.03): filename = 'z.2.80'
	if (z >= 3.03) and (z <= 3.50): filename = 'z.3.25'
	if (z >= 3.50) and (z <= 4.00): filename = 'z.3.75'
	if (z >= 4.00) and (z <= 4.50): filename = 'z.4.25'
	if (z >= 4.50) and (z <= 5.00): filename = 'z.4.75'

# fit the luminosity function based on datasets at a given redshift

# SDSS DR3 (Richards et al. 2006)
	z_min = [0.40, 0.68, 1.06, 1.44, 1.82, 2.21, 2.60, 3.03, 3.50, 4.00, 4.50]
	z_max = [0.68, 1.06, 1.44, 1.82, 2.21, 2.60, 3.03, 3.50, 4.00, 4.50, 5.00]
	z_lis = 0.5*(z_min + z_max)
	for iz in range(len(z_lis)):
		redshift = z_lis[iz]
		L_bol_tmp, phi_fit_tmp = load_sdss_dr3_LF_data(redshift)
			
		phi_fit_pts = inter.interp1d(L_bol_tmp, phi_fit_tmp, L_bol)

			if (RENORM_KEY) then PHI_BB = PHI_BB + (MEAN((phi_fit_pts))-MEAN((PHI_BB)))	
		if (n_elements(L_BB) GT 1) then begin
		convolved_LF, LF_PARAMS, LB, phi, /BB, $
			UEDA_OBSC=UEDA_KEY, SIMPSON_OBSC=SIMP_KEY, LAFRANCA_OBSC=LAFR_KEY
		phi_i = INTERPOL(phi,(LB),alog10(L_BB))
		ALL_P_PRED = [ALL_P_PRED , phi_i]
			ALL_L_OBS  = [ALL_L_OBS  , alog10(L_BB)]
			ALL_P_OBS  = [ALL_P_OBS  , PHI_BB]
			ALL_D_OBS  = [ALL_D_OBS  , DPHI_BB]
			ALL_Z_TOT  = [ALL_Z_TOT  , 0.0*phi_i + REDSHIFT]
			ALL_B      = [ALL_B      , 0*phi_i   + 0]
			ALL_ID     = [ALL_ID     , 0*phi_i   + 12]
		if (keyword_set(DUMP_SPEC_SURVEYS)) then begin
			if (iz EQ 0) then chi2tott = 0.
			if (iz EQ 0) then doftott  = 0
			
			chi2tott = chi2tott + TOTAL(((phi_i-PHI_BB)/DPHI_BB)^2)
			doftott  = doftott  + n_elements(PHI_BB)
			if (iz EQ (n_elements(z_lis)-1)) then print, 'DR3 ->',chi2tott,doftott
		endif
		endif
	endfor 

