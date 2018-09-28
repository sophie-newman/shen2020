;;; 
;;; HMS 2005 LF data (supplants Miyaji LF data)
;;;
pro load_Hasinger_LF_data, L_SX, PHI_SX, DPHI_SX, z
	;; determine which redshift interval its in
	if ((z GE 0.0) AND (z LT 0.2)) then WHICH_BLOCK = 1
	if ((z GE 0.2) AND (z LT 0.4)) then WHICH_BLOCK = 2
	if ((z GE 0.4) AND (z LT 0.8)) then WHICH_BLOCK = 3
	if ((z GE 0.8) AND (z LT 1.6)) then WHICH_BLOCK = 4
	if ((z GE 1.6) AND (z LT 3.2)) then WHICH_BLOCK = 5
	if ((z GE 3.2) AND (z LE 4.8)) then WHICH_BLOCK = 6
	if (z GT 4.8) then WHICH_BLOCK = 7

	if (WHICH_BLOCK EQ 7) then begin
		L_SX = 0.
		PHI_SX = 0.
		DPHI_SX = 0.
	endif else begin
		OPENR, 1, return_idl_routines_homedir(0)+'/luminosity_functions/load_hasinger_lf_data.dat'
		zmin_HX = 0.0
		zmax_HX = 0.0
		N_LUM_BINS = 0
		head = lonarr(2)
		body = fltarr(8)
		for jj = 1, WHICH_BLOCK do begin
			READF, 1, head
			block = head[0]
			N_LUM_BINS = fix(head[1])
			L_SX_min = fltarr(N_LUM_BINS)
			L_SX_max = fltarr(N_LUM_BINS)
			PHI_SX   = fltarr(N_LUM_BINS)
			DPHIp_SX = fltarr(N_LUM_BINS)
			DPHIm_SX = fltarr(N_LUM_BINS)
			N_QSO_SX   = intarr(N_LUM_BINS)
		for i = 0, N_LUM_BINS-1 do begin
			READF,1,body
			zmin_HX     = body[0]
			zmax_HX     = body[1]
			L_SX_min[i] = body[2]
			L_SX_max[i] = body[3]
			PHI_SX[i]   = body[4]
			DPHIp_SX[i] = body[5] 
			DPHIm_SX[i] = body[6] 
			N_QSO_SX[i] = fix(body[7])
		endfor
		endfor	
		CLOSE, 1
		L_SX = 0.5 * (L_SX_min + L_SX_max)

		ok = where((L_SX LT 48.0)); AND (L_SX GT 42.))		;;; things go seriously wonky with the Ueda fits if these are included
		L_SX = L_SX[ok]
		PHI_SX = PHI_SX[ok]
		DPHIp_SX = DPHIp_SX[ok]
		DPHIm_SX = DPHIm_SX[ok]

		;; want to correct for the type-2 fraction, since this is only a Type-1 
		;;  luminosity function: do so following the data in Hasinger et al. 2004
		;;  as compiled in Simpson et al. 2005
		;;
		L0 = 10^(DOUBLE(42.37))/0.015 * 0.606	;; our spectrum HX->SX
		xi = 0.23
		f1 = 1. - ((1.+3.*((10^(DOUBLE(L_SX))/L0)^(1.-2.*xi)))^(-0.5))
		f1 = 1./(1.+ 10^(-DOUBLE(L_SX - 45.0) * 0.4)) * 3.0
		PHI_SX = PHI_SX + alog10(1./f1) 

		DPHI_SX = 0.5 * (DPHIp_SX + abs(DPHIm_SX))
		L_SX = (DOUBLE(L_SX)-alog10(3.9)-33.0)
		
		
	endelse
end
