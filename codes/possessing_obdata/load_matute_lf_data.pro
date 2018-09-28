;;; Matute et al. 2006 -- 15 micron local and z=1.2
;;;
;;;   ADD THE TYPE 2's? -- weak contribution, although up to factor ~2 at low-L
;;;
pro load_Matute_LF_data, L_IR, PHI_IR, DPHI_IR, z
	L_solar = alog10(4.0) + 33.
	if (((z GT 0.2) AND (z LT 0.5)) OR (z GE 2.0)) then begin
		L_IR = 0.
		PHI_IR = 0.
		DPHI_IR = 0.
	endif else begin

	if (z LE 0.2) then begin
		L_15 = [42.6,  43.0,  43.4,  43.8,  44.2,  44.6,  45.0,  45.4,  45.8,  46.2]
		P_15 = [-4.10, -4.30, -4.90, -5.25, -5.74, -6.16, -7.12, -8.05, -8.70, -8.97] 
		D_15 = [ 0.25,  0.15,  0.25,  0.25,  0.25,  0.20,  0.40,  0.50,  0.55,  0.55]
	endif 
	if ((z GE 0.5) AND (z LE 2.0)) then begin
		L_15 = [43.8,  44.2,  44.6,  45.0,  45.4,  45.8,  46.2,  46.6]
		P_15 = [-4.72, -5.25, -5.30, -5.30, -6.13, -6.70, -7.50, -8.52]
		D_15 = [ 0.30,  0.30,  0.20,  0.25,  0.20,  0.20,  0.25,  0.55]	
	endif
	
	;; use below if converting to other bands
	;log_L_15_over_L_R = 0.23
	;L_R = L_15 - log_L_15_over_L_R - (alog10(3.9)+33.0)
	;L_HX   = L_R + alog10(17.5) - 2.*alog10(7./7.5)
	
	L_IR   = (DOUBLE(L_15-L_solar)) - 2.5*alog10(7./7.5) ;;+ 0.25
	PHI_IR = P_15 + alog10(2.5) + 3.*alog10(7./7.5)		 ;;- 0.2
	DPHI_IR= D_15
	endelse
end
