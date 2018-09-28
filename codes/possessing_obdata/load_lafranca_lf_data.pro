;;; 
;;; La Franca et al. 2005 data
;;;
pro load_LaFranca_LF_data, L_HX, PHI_HX, DPHI_HX, z
	;; determine which redshift interval its in
	if ((z GE 0.0) AND (z LT 0.5)) then WHICH_BLOCK = 1
	if ((z GE 0.5) AND (z LT 1.0)) then WHICH_BLOCK = 2
	if ((z GE 1.0) AND (z LT 1.5)) then WHICH_BLOCK = 3
	if ((z GE 1.5) AND (z LT 2.5)) then WHICH_BLOCK = 4
	if ((z GE 2.5) AND (z LE 3.5)) then WHICH_BLOCK = 5
	if (z GT 3.5) then WHICH_BLOCK = 6
	
	if (WHICH_BLOCK EQ 6) then begin
		L_HX = 0.
		PHI_HX = 0.
		DPHI_HX = 0.
	endif else begin
	if (WHICH_BLOCK EQ 1) then begin
		L_HX = [42.5, 43.5, 44.5, 45.5]
		P_HX = [-3.7, -4.6, -6.2, -9.05]
		D_HX = [ 0.2,  0.2,  0.2,  0.8]
	endif
	if (WHICH_BLOCK EQ 2) then begin
		L_HX = [42.5, 43.5, 44.5, 45.5]
		P_HX = [-3.25, -4., -5.45, -8.0]
		D_HX = [ 0.2,  0.2,  0.2,  0.33]
	endif
	if (WHICH_BLOCK EQ 3) then begin
		L_HX = [42.5, 43.5, 44.5, 45.5]
		P_HX = [-3.1, -4.05, -5., -7.0]
		D_HX = [ 0.3,  0.2,  0.2,  0.3]
	endif
	if (WHICH_BLOCK EQ 4) then begin
		L_HX = [43.5, 44.5, 45.5, 46.5]
		P_HX = [-4.1, -4.9, -6.9, -8.95]
		D_HX = [ 0.2,  0.2,  0.3,  0.5]
	endif
	if (WHICH_BLOCK EQ 5) then begin
		L_HX = [44.25, 45.0]
		P_HX = alog10([2.0 * 8.0d-6, 1.0 * 2.0d-6])
		D_HX = [ 0.2, 0.25]
	endif
	L_solar = alog10(4.0) + 33.0
	L_HX = (L_HX - L_solar)
	PHI_HX  = P_HX
	DPHI_HX = D_HX
	endelse
end
