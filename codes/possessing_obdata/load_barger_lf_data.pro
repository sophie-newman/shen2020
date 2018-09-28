;;; 
;;; Barger al. 2005 data
;;;
pro load_Barger_LF_data, L_HX, PHI_HX, DPHI_HX, z
	;; determine which redshift interval its in
	if ((z LT 0.1)) then WHICH_BLOCK = 7
	if ((z GE 0.1) AND (z LT 0.4)) then WHICH_BLOCK = 1
	if ((z GE 0.4) AND (z LT 0.8)) then WHICH_BLOCK = 2
	if ((z GE 0.8) AND (z LT 1.2)) then WHICH_BLOCK = 3
	if ((z GE 1.2) AND (z LT 1.5)) then WHICH_BLOCK = 6
	if ((z GE 1.5) AND (z LT 3.0)) then WHICH_BLOCK = 4
	if ((z GE 3.0) AND (z LT 5.0)) then WHICH_BLOCK = 5
	if ((z GE 5.0) AND (z LE 6.5)) then WHICH_BLOCK = 6
	if (z GT 6.5) then WHICH_BLOCK = 7
	
	if (WHICH_BLOCK EQ 7) then begin
		L_HX = 0.
		PHI_HX = 0.
		DPHI_HX = 0.
	endif else begin
	if (WHICH_BLOCK EQ 1) then begin
		L_HX = [42.25, 42.75, 43.25, 43.75]
		P_HX = [-4.15, -4.10, -4.52, -5.16]
		D_HX = [ 0.25,  0.20,  0.20,  0.30]
	endif
	if (WHICH_BLOCK EQ 2) then begin
		L_HX = [42.25, 42.75, 43.25, 43.75, 44.25, 44.75]
		P_HX = [-3.70, -4.15, -4.30, -4.55, -5.30, -6.40]
		D_HX = [ 0.20,  0.25,  0.35,  0.20,  0.25,  0.40]
	endif
	if (WHICH_BLOCK EQ 3) then begin
		L_HX = [42.25, 42.75, 43.25, 43.75, 44.25, 44.75, 45.25]
		P_HX = [-3.52, -4.05, -4.15, -4.40, -4.74, -5.89, -6.66]
		D_HX = [ 0.20,  0.25,  0.25,  0.25,  0.35,  0.35,  0.40]
	endif
	if (WHICH_BLOCK EQ 4) then begin
		L_HX = [42.25, 42.75, 43.25, 43.75, 44.25, 44.75, 45.25]
		P_HX = [-4.70, -4.75, -5.00, -5.10, -5.22, -5.58, -6.30]
		D_HX = [ 0.80,  0.70,  0.80,  0.25,  0.20,  0.35,  0.45]
	endif
	if (WHICH_BLOCK EQ 5) then begin
		L_HX = [42.75, 43.25, 43.75, 44.25, 44.75, 45.25]
		P_HX = [-5.30, -5.52, -5.30, -6.05, -6.52, -6.75]
		D_HX = [ 0.70,  0.45,  0.40,  0.30,  0.45,  0.60]
	endif
	
	;;; Barger et al. 2003 (highest-redshift) data
	if (WHICH_BLOCK EQ 6) then begin
		L_HX = [43.5]
		P_HX = [-5.95]
		D_HX = [ 0.60]
	endif
		
	L_solar = alog10(4.0) + 33.0
	L_HX = (L_HX - L_solar) + alog10(1.25) ;; convert the Barger 2-8kev to 2-10kev (gamma=1.4)
		;; (note if instead were to use a gamma=1.8 photon index would get 1.19, so very small difference
	PHI_HX  = P_HX
	DPHI_HX = D_HX
	endelse
end
