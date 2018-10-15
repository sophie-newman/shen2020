function return_Brown_LF_fitted, L0_list, z	
	;; given L in L_solar at 15 microns, convert it 
	L_solar = 3.9d33
	L_15_over_L_8 = 10^(+0.060001373)	;; from Hatziminaoglou model spectrum
		;L_15_over_L_8 = 10^(+0.0263996)	;; from Hatziminaoglou model spectrum
	L_8 = (DOUBLE(10^(L0_list))*L_solar)/L_15_over_L_8
	M_8 = -28.5 - 2.5*alog10(L_8/4.11d45)

	phi_29_2 = 4.46d-7 * 2.5
	alpha = -2.75
	k1    =  1.15
	k2    = -0.34
	k3    =  0.03
	zpeak =  2.56

	M_z   = -29.0 - 2.5*(k1*z + k2*z*z + k3*z*z*z) + 2.5*(2.*k1 + 4.*k2 + 8.*k3)
	phi   = phi_29_2 * 10^(0.4*(alpha+1.)*(M_z-M_8))
	
	return, alog10(phi) ;+ alog10(1./(1.-0.55))	
		;; correction for mean Type-2 fraction (NOTE: not 
		;;   obscured in the IR, but not optical qso so 
		;;   no redshift data in this limited sample
end
