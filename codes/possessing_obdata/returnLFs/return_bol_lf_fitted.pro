;; return the best-fit bolometric QLF at a given redshift 
;;   (from the paper with Gordon)
;;
function return_bol_LF_fitted, log_L, z, SCHECHTER=SCHECHTER, PLE=PLE, NLP=NLP, $
	RETURN_BESTFIT_PARAMS=RETURN_BESTFIT_PARAMS, $
	RETURN_BESTFIT_ERRORS=RETURN_BESTFIT_ERRORS, $
	RETURN_BESTFIT_COVARIANCE=RETURN_BESTFIT_COVARIANCE, $
	RETURN_BESTFIT_CHI2=RETURN_BESTFIT_CHI2, $
	RETURN_BESTFIT_DOF=RETURN_BESTFIT_DOF
		COMMON BOL_LF_FITTING_PARAMS
		forward_function LF_shape_vs_z
		REDSHIFT = z
	
	homedir = return_idl_routines_homedir(0)+'/luminosity_functions/'
	filena = homedir+'full.fits.params.dat'
		if (keyword_set(SCHECHTER)) then filena=homedir+'param.fits.schechter.dat'
		if (keyword_set(PLE)) then filena=homedir+'param.fits.ple.dat'
		if (keyword_set(NLP)) then filena=homedir+'param.fits.nlp.dat'
		  if (keyword_set(SCHECHTER)) then MODIFIED_SCHECHTER_FIT = 1
		
		OPENR,lunner,filena,/get_lun
			READF,lunner,n_params
			LF_FIT_PARAMS=fltarr(n_params)
			LF_FIT_ERRORS=fltarr(n_params)
			LF_FIT_COVARIANCE=fltarr(n_params,n_params)
			LF_FIT_BESTNORM=0.
			LF_FIT_DOF=0
			READF,lunner,LF_FIT_PARAMS
			READF,lunner,LF_FIT_ERRORS
			READF,lunner,LF_FIT_COVARIANCE
			READF,lunner,LF_FIT_BESTNORM
			READF,lunner,LF_FIT_DOF
		CLOSE,lunner	
		FREE_LUN,lunner
	
	if (keyword_set(RETURN_BESTFIT_PARAMS)) 	then RETURN_BESTFIT_PARAMS = LF_FIT_PARAMS
	if (keyword_set(RETURN_BESTFIT_ERRORS)) 	then RETURN_BESTFIT_ERRORS = LF_FIT_ERRORS
	if (keyword_set(RETURN_BESTFIT_COVARIANCE)) then RETURN_BESTFIT_COVARIANCE = LF_FIT_COVARIANCE
	if (keyword_set(RETURN_BESTFIT_CHI2)) 		then RETURN_BESTFIT_CHI2 = LF_FIT_CHI2
	if (keyword_set(RETURN_BESTFIT_DOF)) 		then RETURN_BESTFIT_DOF = LF_FIT_DOF

	f = LF_shape_vs_z(log_L,LF_FIT_PARAMS)
	return, f
end
