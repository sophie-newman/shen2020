;;; COMBO-17 (Wolf et al.)
;;;   measuring M_1450
;;;
pro load_COMBO17_LF_data, L_BB, PHI_BB, DPHI_BB, z
	forward_function M_sun_BB

	if ((z LT 1.2) OR (z GT 4.8)) then begin
		L_BB = 0.
		PHI_BB = 0.
		DPHI_BB = 0.
	endif else begin
	if ((z GE 1.2) AND (z LT 1.8)) then WHICH_BLOCK = 0
	if ((z GE 1.8) AND (z LT 2.4)) then WHICH_BLOCK = 1
	if ((z GE 2.4) AND (z LT 3.0)) then WHICH_BLOCK = 2
	if ((z GE 3.0) AND (z LT 3.6)) then WHICH_BLOCK = 3
	if ((z GE 3.6) AND (z LT 4.2)) then WHICH_BLOCK = 4
	if ((z GE 4.2) AND (z LE 4.8)) then WHICH_BLOCK = 5

	N_Z_BINS = 6
	N_LUM_BINS = 6
	M_1450 = [-23.0,-24.0,-25.0,-26.0,-27.0,-28.0]
	M_B = M_1450 + 1.75 

	n_qso = intarr(N_LUM_BINS,N_Z_BINS)
	phi   = fltarr(N_LUM_BINS,N_Z_BINS)
	dphi  = fltarr(N_LUM_BINS,N_Z_BINS)

	n_qso[*,0] = [7,16,13,6,5,1]
	n_qso[*,1] = [9,22,17,10,5,1]
	n_qso[*,2] = [0,17,14,7,3,3]
	n_qso[*,3] = [0,5,5,3,3,0]
	n_qso[*,4] = [0,3,6,1,1,0]
	n_qso[*,5] = [0,0,0,3,1,0]
	
	phi[*,0] = [-5.78,-5.58,-5.69,-6.02,-6.10,-6.80]
	phi[*,1] = [-5.53,-5.36,-5.60,-5.85,-6.15,-6.85]
	phi[*,2] = [ 0.00,-5.42,-5.69,-6.01,-6.38,-6.38]
	phi[*,3] = [ 0.00,-5.85,-6.10,-6.36,-6.36, 0.00]
	phi[*,4] = [ 0.00,-6.02,-6.02,-6.81,-6.81, 0.00]
	phi[*,5] = [ 0.00, 0.00, 0.00,-6.33,-6.79, 0.00]

	dphi[*,0] = [ 0.10, 0.10, 0.15, 0.17, 0.50, 1.00]
	dphi[*,1] = [ 0.10, 0.10, 0.13, 0.17, 0.50, 0.50]
	dphi[*,2] = [ 0.00, 0.10, 0.10, 0.15, 0.25, 0.30]
	dphi[*,3] = [ 0.00, 0.20, 0.30, 0.30, 0.50, 0.00]
	dphi[*,4] = [ 0.00, 0.17, 0.17, 0.50, 0.50, 0.00]
	dphi[*,5] = [ 0.00, 0.00, 0.00, 0.35, 0.60, 0.00]

	P0  = phi[*,WHICH_BLOCK]
	N0  = n_qso[*,WHICH_BLOCK]
	D0  = dphi[*,WHICH_BLOCK]
	ok  = where(N0 GT 0)
	PHI_BB = P0[ok]
	M_BB   = M_B[ok]
	DPHI_BB= D0[ok]
	PHI_BB = PHI_BB[0:n_elements(PHI_BB)-1] + alog10(2.5) +3.*alog10(7./6.5) ;; converts to per magnitude
	M_BB   = M_BB[0:n_elements(M_BB)-1]+5.*alog10(7./6.5)
	DPHI_BB= DPHI_BB[0:n_elements(DPHI_BB)-1]
	L_BB   = (0.4*(M_sun_BB(0)-M_BB))
	
	endelse
end
