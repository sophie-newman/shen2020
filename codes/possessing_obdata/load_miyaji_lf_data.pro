pro load_Miyaji_LF_data, L_SX, PHI_SX, DPHI_SX, z
	;; determine which redshift interval its in
	if ((z GE 0.0) AND (z LT 0.2)) then WHICH_BLOCK = 1
	if ((z GE 0.2) AND (z LT 0.4)) then WHICH_BLOCK = 2
	if ((z GE 0.4) AND (z LT 0.8)) then WHICH_BLOCK = 3
	if ((z GE 0.8) AND (z LT 1.6)) then WHICH_BLOCK = 4
	if ((z GE 1.6) AND (z LT 2.3)) then WHICH_BLOCK = 5
	if ((z GE 2.3) AND (z LE 4.6)) then WHICH_BLOCK = 6
	if (z GT 4.6) then WHICH_BLOCK = 7

	if (WHICH_BLOCK EQ 7) then begin
		L_SX = 0.
		PHI_SX = 0.
		DPHI_SX = 0.
	endif else begin
		OPENR, 1, return_idl_routines_homedir(0)+'/luminosity_functions/load_miyaji_lf_data.dat'
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
		L_SX = 10^(DOUBLE(L_SX) - alog10(3.9) - 33.) * (0.7/0.5)^(-2.)
		DPHI_SX = 0.5 * (abs(DPHIp_SX) + abs(DPHIm_SX))
		ok = where(DPHI_SX GT 0.)
		L_SX = L_SX[ok]
		DPHI_SX = DPHI_SX[ok]
		PHI_SX = PHI_SX[ok] + 3.*alog10(7./5.)

		ok = where(alog10(L_SX) LT 47.0)		;;; things go seriously wonky with the Ueda fits if these are included
		L_SX = alog10(L_SX[ok])
		PHI_SX = PHI_SX[ok]
		DPHI_SX = DPHI_SX[ok]
	endelse
end
