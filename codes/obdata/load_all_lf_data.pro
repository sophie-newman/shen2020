pro load_ALL_LF_data, z, L, P, D, B, ID, $
	RENORM_KEY=RENORM_KEY, ALL_RENORM_KEY=ALL_RENORM_KEY

	forward_function fitted_band_lum, M_sun_BB, return_SDSS_DR3_LF_fitted, $
		return_2SLAQ_LF_fitted, return_COMBO17_LF_fitted, return_Hasinger_LF_fitted, $
		return_Miyaji_LF_fitted, return_Ueda_LF_fitted, return_LaFranca_LF_fitted, $
		return_Brown_LF_fitted, return_Matute_LF_fitted

	all_L = [0.]	;; log(luminosity/solar) to return
	all_P = [0.]	;; log(phi/logL/Mpc^3)   to return
	all_D = [0.]	;; d(log(phi))
	all_B = [0]		;; BANDs : (0) = B, (1) = SX, (2) = HX, (3) = 15mu, (4) = narrow-line(HX)
	all_ID= [0]		;; ID numbers (different for each sample, for e.g. plotting, etc.)

	if (keyword_set(RENORM_KEY)) then begin
		RENORM_KEY = 1
	endif else begin
		RENORM_KEY = 0
	endelse
	if (keyword_set(ALL_RENORM_KEY)) then begin
		ALL_RENORM_KEY = 1
	endif else begin
		ALL_RENORM_KEY = 0
	endelse

	L_tmp = (7.0 + 0.1*findgen(111))
		L_BB_tmp = DOUBLE(fitted_band_lum(L_tmp,/BB))
		L_SX_tmp = DOUBLE(fitted_band_lum(L_tmp,/SX))
		L_HX_tmp = DOUBLE(fitted_band_lum(L_tmp,/HX))
		L_IR_tmp = DOUBLE(fitted_band_lum(L_tmp,/IR))


	;;;;;;
	;; 	OPTICAL (B-BAND)
	;;;;;;	

	;; SDSS DR3 (Richards et al. 2006)
	load_SDSS_DR3_LF_data, L_BB, PHI_BB, DPHI_BB, z
		phi_fit_tmp = return_SDSS_DR3_LF_fitted(L_BB_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_BB_tmp,L_BB)
		if (RENORM_KEY) then PHI_BB = PHI_BB + (MEAN((phi_fit_pts))-MEAN((PHI_BB)))		
		if (L_BB(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_BB]
			ALL_P  = [ALL_P  , PHI_BB]
			ALL_D  = [ALL_D  , DPHI_BB]
			ALL_B  = [ALL_B  , 0*L_BB   + 0]
			ALL_ID = [ALL_ID , 0*L_BB   + 12]
		endif	

	;; 2SLAQ (Richards et al. 2005)
	load_2SLAQ_LF_data, L_BB, PHI_BB, DPHI_BB, z
		phi_fit_tmp = return_2SLAQ_LF_fitted(L_BB_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_BB_tmp,L_BB)
		if ((z LE 0.4) OR (z GE 0.6)) then $
		if (RENORM_KEY) then PHI_BB = PHI_BB + (MEAN((phi_fit_pts))-MEAN((PHI_BB)))	
		if (L_BB(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_BB]
			ALL_P  = [ALL_P  , PHI_BB]
			ALL_D  = [ALL_D  , DPHI_BB]
			ALL_B  = [ALL_B  , 0*L_BB   + 0]
			ALL_ID = [ALL_ID , 0*L_BB   + 11]
		endif	

	;; Hunt (LBGs - to be updated shortly) (Hunt et al. 2003)
	load_Hunt_LF_data, L_BB, PHI_BB, DPHI_BB, z
		phi_fit_tmp = return_COMBO17_LF_fitted(L_BB_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_BB_tmp,L_BB)
		;if (RENORM_KEY) then PHI_BB = PHI_BB + (MEAN((phi_fit_pts))-MEAN((PHI_BB)))	
			if (ALL_RENORM_KEY) then PHI_BB = PHI_BB + (MEAN((phi_fit_pts))-MEAN((PHI_BB)))	
		if (L_BB(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_BB]
			ALL_P  = [ALL_P  , PHI_BB]
			ALL_D  = [ALL_D  , DPHI_BB]
			ALL_B  = [ALL_B  , 0*L_BB   + 0]
			ALL_ID = [ALL_ID , 0*L_BB   + 14]
		endif	

	;; Siana et al. 2006
	load_Siana_LF_data, L_BB, PHI_BB, DPHI_BB, z
		phi_fit_tmp = return_COMBO17_LF_fitted(L_BB_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_BB_tmp,L_BB)
		;if (RENORM_KEY) then PHI_BB = PHI_BB + (MEAN((phi_fit_pts))-MEAN((PHI_BB)))	
			if (ALL_RENORM_KEY) then PHI_BB = PHI_BB + (MEAN((phi_fit_pts))-MEAN((PHI_BB)))	
		if (L_BB(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_BB]
			ALL_P  = [ALL_P  , PHI_BB]
			ALL_D  = [ALL_D  , DPHI_BB]
			ALL_B  = [ALL_B  , 0*L_BB   + 0]
			ALL_ID = [ALL_ID , 0*L_BB   + 23]
		endif	

	;; Cristiani et al. 2004
	load_Cristiani_LF_data, L_BB, PHI_BB, DPHI_BB, z
		phi_fit_tmp = return_COMBO17_LF_fitted(L_BB_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_BB_tmp,L_BB)
		if (RENORM_KEY) then PHI_BB = PHI_BB + (MEAN((phi_fit_pts))-MEAN((PHI_BB)))	
			if (ALL_RENORM_KEY) then PHI_BB = PHI_BB + (MEAN((phi_fit_pts))-MEAN((PHI_BB)))	
		if (L_BB(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_BB]
			ALL_P  = [ALL_P  , PHI_BB]
			ALL_D  = [ALL_D  , DPHI_BB]
			ALL_B  = [ALL_B  , 0*L_BB   + 0]
			ALL_ID = [ALL_ID , 0*L_BB   + 15]
		endif	

	;; KDC LF
	load_KDC_LF_data, L_BB, PHI_BB, DPHI_BB, z
		phi_fit_tmp = return_SDSS_DR3_LF_fitted(L_BB_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_BB_tmp,L_BB)
		;if (RENORM_KEY) then PHI_BB = PHI_BB + (MEAN((phi_fit_pts))-MEAN((PHI_BB)))	
			if (ALL_RENORM_KEY) then PHI_BB = PHI_BB + (MEAN((phi_fit_pts))-MEAN((PHI_BB)))	
		if (L_BB(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_BB]
			ALL_P  = [ALL_P  , PHI_BB]
			ALL_D  = [ALL_D  , DPHI_BB]
			ALL_B  = [ALL_B  , 0*L_BB   + 0]
			ALL_ID = [ALL_ID , 0*L_BB   + 16]
		endif	

	;; Fan et al. 2001-2004
	load_SDSS_FAN_LF_data, L_BB, PHI_BB, DPHI_BB, z
		phi_fit_tmp = return_SDSS_FAN_LF_fitted(L_BB_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_BB_tmp,L_BB)
		if (RENORM_KEY) then PHI_BB = PHI_BB + (MEAN((phi_fit_pts))-MEAN((PHI_BB)))	
		if (L_BB(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_BB]
			ALL_P  = [ALL_P  , PHI_BB]
			ALL_D  = [ALL_D  , DPHI_BB]
			ALL_B  = [ALL_B  , 0*L_BB   + 0]
			ALL_ID = [ALL_ID , 0*L_BB   + 18]
		endif	

	;; SSG LF
	load_SSG_LF_data, L_BB, PHI_BB, DPHI_BB, z
		phi_fit_tmp = return_SDSS_DR3_LF_fitted(L_BB_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_BB_tmp,L_BB)
		if (RENORM_KEY) then PHI_BB = PHI_BB + (MEAN((phi_fit_pts))-MEAN((PHI_BB)))	
			if (ALL_RENORM_KEY) then PHI_BB = PHI_BB + (MEAN((phi_fit_pts))-MEAN((PHI_BB)))	
		if (L_BB(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_BB]
			ALL_P  = [ALL_P  , PHI_BB]
			ALL_D  = [ALL_D  , DPHI_BB]
			ALL_B  = [ALL_B  , 0*L_BB   + 0]
			ALL_ID = [ALL_ID , 0*L_BB   + 17]
		endif	

	;; COMBO-17 (Wolf et al.)
	load_COMBO17_LF_data, L_BB, PHI_BB, DPHI_BB, z
		phi_fit_tmp = return_COMBO17_LF_fitted(L_BB_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_BB_tmp,L_BB)
		if (RENORM_KEY) then PHI_BB = PHI_BB + (MEAN((phi_fit_pts))-MEAN((PHI_BB)))	
		if (L_BB(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_BB]
			ALL_P  = [ALL_P  , PHI_BB]
			ALL_D  = [ALL_D  , DPHI_BB]
			ALL_B  = [ALL_B  , 0*L_BB   + 0]
			ALL_ID = [ALL_ID , 0*L_BB   + 13]
		endif	

	;;;;;;
	;; 	SOFT X-RAYS
	;;;;;;	

	;; HMS 2005 (supplants Miyaji et al. LFs for Type I's, crudely corrected for Type-II's)
	load_Hasinger_LF_data, L_SX, PHI_SX, DPHI_SX, z
		phi_fit_tmp = return_Hasinger_LF_fitted(L_SX_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_SX_tmp,L_SX)
		if (RENORM_KEY) then PHI_SX = PHI_SX + (MEAN((phi_fit_pts))-MEAN((PHI_SX)))	
		if (L_SX(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_SX]
			ALL_P  = [ALL_P  , PHI_SX]
			ALL_D  = [ALL_D  , DPHI_SX]
			ALL_B  = [ALL_B  , 0*L_SX   + 1]
			ALL_ID = [ALL_ID , 0*L_SX   + 9]
		endif	

	;; Miyaji et al. 2001
	load_Miyaji_LF_data, L_SX, PHI_SX, DPHI_SX, z
		phi_fit_tmp = return_Miyaji_LF_fitted(L_SX_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_SX_tmp,L_SX)
		if (RENORM_KEY) then PHI_SX = PHI_SX + (MEAN((phi_fit_pts))-MEAN((PHI_SX)))	
		if (L_SX(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_SX]
			ALL_P  = [ALL_P  , PHI_SX]
			ALL_D  = [ALL_D  , DPHI_SX]
			ALL_B  = [ALL_B  , 0*L_SX   + 1]
			ALL_ID = [ALL_ID , 0*L_SX   + 21]
		endif	

	;; Silverman et al. 2005
	load_Silverman_SX_LF_data, L_SX, PHI_SX, DPHI_SX, z
		phi_fit_tmp = return_Hasinger_LF_fitted(L_SX_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_SX_tmp,L_SX)
		if (RENORM_KEY) then PHI_SX = PHI_SX + (MEAN((phi_fit_pts))-MEAN((PHI_SX)))	
		if (L_SX(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_SX]
			ALL_P  = [ALL_P  , PHI_SX]
			ALL_D  = [ALL_D  , DPHI_SX]
			ALL_B  = [ALL_B  , 0*L_SX   + 1]
			ALL_ID = [ALL_ID , 0*L_SX   + 10]
		endif	

	;;;;;;
	;; 	HARD X-RAYS
	;;;;;;	

	;; Ueda et al. 2003
	load_Ueda_LF_data, L_HX, PHI_HX, DPHI_HX, z
		phi_fit_tmp = return_Ueda_LF_fitted(L_HX_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_HX_tmp,L_HX)
		if (RENORM_KEY) then PHI_HX = PHI_HX + (MEAN((phi_fit_pts))-MEAN((PHI_HX)))	
		if (L_HX(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_HX]
			ALL_P  = [ALL_P  , PHI_HX]
			ALL_D  = [ALL_D  , DPHI_HX]
			ALL_B  = [ALL_B  , 0*L_HX   + 2]
			ALL_ID = [ALL_ID , 0*L_HX   + 1]
		endif	

	;; La Franca et al. 2005
	load_LaFranca_LF_data, L_HX, PHI_HX, DPHI_HX, z
		phi_fit_tmp = return_LaFranca_LF_fitted(L_HX_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_HX_tmp,L_HX)
		if (RENORM_KEY) then PHI_HX = PHI_HX + (MEAN((phi_fit_pts))-MEAN((PHI_HX)))	
		if (L_HX(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_HX]
			ALL_P  = [ALL_P  , PHI_HX]
			ALL_D  = [ALL_D  , DPHI_HX]
			ALL_B  = [ALL_B  , 0*L_HX   + 2]
			ALL_ID = [ALL_ID , 0*L_HX   + 1]
		endif	

	;; Silverman et al. 2005
	load_Silverman_HX_LF_data, L_HX, PHI_HX, DPHI_HX, z
		phi_fit_tmp = return_Ueda_LF_fitted(L_HX_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_HX_tmp,L_HX)
		;if (RENORM_KEY) then PHI_HX = PHI_HX + (MEAN((phi_fit_pts))-MEAN((PHI_HX)))	
			if (ALL_RENORM_KEY) then PHI_HX = PHI_HX + (MEAN((phi_fit_pts))-MEAN((PHI_HX)))	
		if (L_HX(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_HX]
			ALL_P  = [ALL_P  , PHI_HX]
			ALL_D  = [ALL_D  , DPHI_HX]
			ALL_B  = [ALL_B  , 0*L_HX   + 2]
			ALL_ID = [ALL_ID , 0*L_HX   + 3]
		endif	

	;; Barger et al. 2005
	load_Barger_LF_data, L_HX, PHI_HX, DPHI_HX, z
		phi_fit_tmp = return_Ueda_LF_fitted(L_HX_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_HX_tmp,L_HX)
		if (RENORM_KEY) then PHI_HX = PHI_HX + (MEAN((phi_fit_pts))-MEAN((PHI_HX)))	
		if (L_HX(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_HX]
			ALL_P  = [ALL_P  , PHI_HX]
			ALL_D  = [ALL_D  , DPHI_HX]
			ALL_B  = [ALL_B  , 0*L_HX   + 2]
			ALL_ID = [ALL_ID , 0*L_HX   + 4]
		endif	

	;; Nandra et al. 2005
	load_Nandra_LF_data, L_HX, PHI_HX, DPHI_HX, z
		if (L_HX(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_HX]
			ALL_P  = [ALL_P  , PHI_HX]
			ALL_D  = [ALL_D  , DPHI_HX]
			ALL_B  = [ALL_B  , 0*L_HX   + 2]
			ALL_ID = [ALL_ID , 0*L_HX   + 6]
		endif	

	;; Sazonov & Revnivtsev et al. 2004
	load_SazRev_LF_data, L_HX, PHI_HX, DPHI_HX, z
		phi_fit_tmp = return_Ueda_LF_fitted(L_HX_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_HX_tmp,L_HX)
		if (RENORM_KEY) then PHI_HX = PHI_HX + (MEAN((phi_fit_pts))-MEAN((PHI_HX)))	
		if (L_HX(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_HX]
			ALL_P  = [ALL_P  , PHI_HX]
			ALL_D  = [ALL_D  , DPHI_HX]
			ALL_B  = [ALL_B  , 0*L_HX   + 2]
			ALL_ID = [ALL_ID , 0*L_HX   + 7]
		endif	

	;;;;;;
	;; 	INFRARED
	;;;;;;	

	;; Brown et al. 2005
	load_Brown_LF_data, L_IR, PHI_IR, DPHI_IR, z
		phi_fit_tmp = return_Brown_LF_fitted(L_IR_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_IR_tmp,L_IR)
		if (RENORM_KEY) then PHI_IR = PHI_IR + (MEAN((phi_fit_pts))-MEAN((PHI_IR)))	
		if (L_IR(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_IR]
			ALL_P  = [ALL_P  , PHI_IR]
			ALL_D  = [ALL_D  , DPHI_IR]
			ALL_B  = [ALL_B  , 0*L_IR   + 3]
			ALL_ID = [ALL_ID , 0*L_IR   + 20]
		endif	

	;; Matute et al. 2006
	load_Matute_LF_data, L_IR, PHI_IR, DPHI_IR, z
		phi_fit_tmp = return_Matute_LF_fitted(L_IR_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_IR_tmp,L_IR)
		if (RENORM_KEY) then PHI_IR = PHI_IR + (MEAN((phi_fit_pts))-MEAN((PHI_IR)))	
		if (L_IR(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_IR]
			ALL_P  = [ALL_P  , PHI_IR]
			ALL_D  = [ALL_D  , DPHI_IR]
			ALL_B  = [ALL_B  , 0*L_IR   + 3]
			ALL_ID = [ALL_ID , 0*L_IR   + 19]
		endif	

	;;;;;;
	;; 	NARROW LINES
	;;;;;;	

	;; Hao et al. 2005
	load_Hao_LF_data, L_HX, PHI_HX, DPHI_HX, z
		phi_fit_tmp = return_Ueda_LF_fitted(L_HX_tmp,z)
		phi_fit_pts = INTERPOL(phi_fit_tmp,L_HX_tmp,L_HX)
		;if (RENORM_KEY) then PHI_HX = PHI_HX + (MEAN((phi_fit_pts))-MEAN((PHI_HX)))	
		;	if (ALL_RENORM_KEY) then PHI_HX = PHI_HX + (MEAN((phi_fit_pts))-MEAN((PHI_HX)))	
		if (L_HX(0) NE 0.) then begin
			ALL_L  = [ALL_L  , L_HX]
			ALL_P  = [ALL_P  , PHI_HX]
			ALL_D  = [ALL_D  , DPHI_HX]
			ALL_B  = [ALL_B  , 0*L_HX   + 4]
			ALL_ID = [ALL_ID , 0*L_HX   + 8]
		endif	


	ALL_L  = ALL_L[1:n_elements(ALL_L)-1]
	ALL_P  = ALL_P[1:n_elements(ALL_P)-1]
	ALL_D  = ALL_D[1:n_elements(ALL_D)-1]
	ALL_B  = ALL_B[1:n_elements(ALL_B)-1]
	ALL_ID =ALL_ID[1:n_elements(ALL_ID)-1]
	ok = where((FINITE(all_P) EQ 1) AND (all_P NE 0.) AND $
		(FINITE(all_P,/NAN) EQ 0) AND (all_L NE 0.) $
		AND (all_P GT -10.), n_data_pts)
	if (n_data_pts GT 0) then begin
		all_L = all_L[ok]
		all_P = all_P[ok]
		all_D = all_D[ok]
		all_B = all_B[ok]
		all_ID = all_ID[ok]
	endif
	L = all_L
	P = all_P
	D = all_D
	B = all_B
	ID= all_ID
end
