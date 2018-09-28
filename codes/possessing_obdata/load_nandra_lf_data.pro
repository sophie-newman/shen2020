;;; 
;;; Nandra et al. 2005 (LBG-selected) data
;;;
pro load_Nandra_LF_data, L_HX, PHI_HX, DPHI_HX, z
	;; determine which redshift interval its in
	WHICH_BLOCK = 2
	if ((z GT 2.75) AND (z LT 3.25)) then WHICH_BLOCK = 1
	
	if (WHICH_BLOCK EQ 2) then begin
		L_HX = 0.
		PHI_HX = 0.
		DPHI_HX = 0.
	endif else begin
		L_HX = [0.5*(43.0+44.5)]
		P_HX = alog10([4.2d-5]/1.5)
		D_HX = [0.15]
	L_solar = alog10(4.0) + 33.0
	L_HX = (L_HX - L_solar) 
	PHI_HX  = P_HX
	DPHI_HX = D_HX
	endelse
end
