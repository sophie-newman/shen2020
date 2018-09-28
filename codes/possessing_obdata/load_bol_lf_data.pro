;; function to return the Nobs/Nmdl binned bolometric QLF w.r.t. the 
;;   best-fit at some redshift, using the full data set available
;;
pro load_BOL_LF_data, z, log_l_bol, log_phi, d_log_phi, band_id, sample_id, QUICK=QUICK, $
		CHI_RETURN=CHI_RETURN, PARAM_RETURN=PARAM_RETURN, DPARAM_RETURN=DPARAM_RETURN, $
		CPARAM_RETURN=CPARAM_RETURN, NO_RENORM=NO_RENORM
	COMMON BOL_LF_FITTING_PARAMS
	forward_function LF_shape_function, LF_fit_at_z, load_BOL_LF_FITTING_PARAMS, $
		convolved_LF, LF_shape_vs_z, lut_ir_NH, return_bol_LF_fitted, LF_params_at_z_fitted
	;; be sure to call dummy = load_BOL_LF_FITTING_PARAMS(1) before running

	load_ALL_LF_data, z, all_L, all_P, all_D, all_B, all_ID, /RENORM_KEY
	if keyword_set(NO_RENORM) then $
			load_ALL_LF_data, z, all_L, all_P, all_D, all_B, all_ID
		BAND_LIST = all_B
		REDSHIFT  = z
		all_z 	  = 0.0*all_L + REDSHIFT

	if (keyword_set(QUICK)) then begin
		P_fitted = LF_params_at_z_fitted(z)
		P_errors = 0.0*P_fitted
		P_covar  = fltarr(n_elements(P_fitted),n_elements(P_fitted))
		chi2 = 0.
		nu   = 1	
	endif else begin
	;; prep the fitter : 
	P_guess = [0.5,2.2,-5.0,12.0]		;; for double power-law
		if (z GE 3.0) then P_guess = [0.0, 1.3, -5.0, 12.0]
		if (z GE 5.0) then P_guess = [0.0, 1.3, -5.0, 11.0]
		P_errors = 0.0*P_guess
		p_fix  = {fixed:1, limited:[0,0], limits:[0.0,0.0]}
		p_gt0  = {fixed:0, limited:[1,0], limits:[0.0,0.0]}
		p_nocon= {fixed:0, limited:[0,0], limits:[0.0,0.0]}
		p_info = [p_nocon,p_nocon,p_nocon,p_nocon]
		if (z GE 3.0) then p_info = [p_fix,p_nocon,p_nocon,p_nocon]
		if (z GE 5.0) then p_info = [p_fix,p_nocon,p_nocon,p_fix]
		P_covar = fltarr(n_elements(P_guess),n_elements(P_guess))
	;; get the best-fit model bolometric QLF
	P_fitted = MPFITFUN('LF_fit_at_z',all_L,all_P,all_D,P_guess, $
				PERROR=P_errors,PARINFO=p_info,COVAR=P_covar,BESTNORM=chi2,DOF=nu,/QUIET)
	endelse

	PARAM = P_fitted
	D_PARAM = P_errors * (chi2/float(nu))
	PARAM_COVAR = P_covar * (chi2/float(nu))
	ALL_BESTNORM = chi2
	ALL_DOF = nu
	print, ' Bol LF at z = ',z,' given by best-fit parameters : '
		print, ' faint-end  slope : ',PARAM[0],' +/- ',D_PARAM[0]
		print, ' bright-end slope : ',PARAM[1],' +/- ',D_PARAM[1]
		print, ' log(phi_star)    : ',PARAM[2],' +/- ',D_PARAM[2]
		print, ' log(l_star)      : ',PARAM[3],' +/- ',D_PARAM[3]
		if (keyword_set(PARAM_RETURN))  then  PARAM_RETURN=PARAM
		if (keyword_set(DPARAM_RETURN)) then DPARAM_RETURN=D_PARAM
		if (keyword_set(CPARAM_RETURN)) then CPARAM_RETURN=PARAM_COVAR


	N_MODEL = 10^(LF_fit_at_z(all_L,P_fitted))
	N_OBS   = 10^(all_P)
		chi2 = TOTAL(((alog10(N_MODEL/N_OBS))/(all_D+0.015))^2)
		nu   = n_elements(N_OBS) - n_elements(P_guess)
		reduced_chi2 = chi2/nu
			if keyword_set(CHI_RETURN) then CHI_RETURN=reduced_chi2
		nustr = textoidl("\chi^{2}/\nu = ")+STRING(reduced_chi2,FORMAT='(F4.1)')
		print, '  ...reduced chi2 = '+STRING(reduced_chi2,FORMAT='(F5.2)'),'   ( ',ALL_BESTNORM,nu,' )'
	L_temp  = 7.0 + 0.1*findgen(111)
	N_full  = 10^(LF_shape_function(L_temp,P_fitted))

	;; correct points to bolometric luminosities
	l_bol_grid = 6. + 0.1*findgen(121)
	all_L_bol  = 0.0*all_L
	x_i = where(all_B EQ 0, n_i)
		if (n_i GT 0) then all_L_bol[x_i] = INTERPOL(l_bol_grid,fitted_band_lum(l_bol_grid,/BB),all_L[x_i])
	x_i = where(all_B EQ 1, n_i)
		if (n_i GT 0) then all_L_bol[x_i] = INTERPOL(l_bol_grid,fitted_band_lum(l_bol_grid,/SX),all_L[x_i])
	x_i = where(all_B EQ 2, n_i)
		if (n_i GT 0) then all_L_bol[x_i] = INTERPOL(l_bol_grid,fitted_band_lum(l_bol_grid,/HX),all_L[x_i])
	x_i = where(all_B EQ 3, n_i)
		if (n_i GT 0) then all_L_bol[x_i] = INTERPOL(l_bol_grid,fitted_band_lum(l_bol_grid,/IR),all_L[x_i])
	x_i = where(all_B EQ 4, n_i)
		if (n_i GT 0) then all_L_bol[x_i] = INTERPOL(l_bol_grid,fitted_band_lum(l_bol_grid,/HX),all_L[x_i])
	;; allow for Elvis corrections
	if (BOLCOR_KEY EQ 1) then $
		all_L_bol = (all_L + alog10(11.8))*(all_B EQ 0) + $
					(all_L + alog10(35.0*1.5))*(all_B EQ 1) + $
					(all_L + alog10(35.0))*(all_B EQ 2) + $
					(all_L + alog10(35.0*1.5))*(all_B EQ 3) + $
					(all_L + alog10(9.1))*(all_B EQ 4) 
	N_intp  = 10^(INTERPOL(alog10(N_full),L_temp,all_L_bol))
	N_renorm = alog10((N_OBS/N_MODEL) * N_intp)
		
		
	log_l_bol = all_L_bol
	log_phi   = N_renorm
	d_log_phi = all_D
	band_id   = all_B
	sample_id = all_ID
end
