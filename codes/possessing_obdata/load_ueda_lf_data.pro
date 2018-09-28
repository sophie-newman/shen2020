;;; 
;;; Ueda et al. 2003 data
;;;
pro load_Ueda_LF_data, L_HX, PHI_HX, DPHI_HX, z
	;; determine which redshift interval its in
	if ((z GE 0.0) AND (z LT 0.2)) then WHICH_BLOCK = 1
	if ((z GE 0.2) AND (z LT 0.4)) then WHICH_BLOCK = 2
	if ((z GE 0.4) AND (z LT 0.8)) then WHICH_BLOCK = 3
	if ((z GE 0.8) AND (z LT 1.6)) then WHICH_BLOCK = 4
	if ((z GE 1.6) AND (z LE 3.0)) then WHICH_BLOCK = 5
	if (z GT 3.0) then WHICH_BLOCK = 6
	
	if (WHICH_BLOCK EQ 6) then begin
		L_HX = 0.
		PHI_HX = 0.
		DPHI_HX = 0.
	endif else begin
		OPENR, 1, return_idl_routines_homedir(0)+'/luminosity_functions/load_ueda_lf_data.dat'
		zmin_SX = 0.0
		zmax_SX = 0.0
		N_LUM_BINS = 0
		head = fltarr(3)
		body = fltarr(5)
		for jj = 1, WHICH_BLOCK do begin
			READF, 1, head
			zmin_SX=head[0]
			zmax_SX=head[1]
			N_LUM_BINS=fix(head[2])
			L_HX = fltarr(N_LUM_BINS)
			PHI_HX = fltarr(N_LUM_BINS)
			DPHI_HX = fltarr(N_LUM_BINS)
			for i = 0, N_LUM_BINS-1 do begin
				READF,1,body
				L_HX[i] = body[2]
				PHI_HX[i]  = body[3]
				DPHI_HX[i] = body[4] 
			endfor
		endfor
		CLOSE, 1
		DPHIm_HX = -1.0*DPHI_HX
		DPHIp_HX = DPHI_HX
		L_solar = alog10(4.0) + 33.0
		L_HX = (L_HX - L_solar)
	endelse
end

