;;; 
;;; Shinozaki et al. 2006 local data
;;;
pro load_Shinozaki_LF_data, L_HX, PHI_HX, DPHI_HX, z
	;; determine which redshift interval its in
	if (z GT 0.20) then begin
		L_HX = 0.
		PHI_HX = 0.
		DPHI_HX = 0.
	endif else begin
	 l = [ 42.5, 43.5, 44.25, 44.85, 45.5 - 0.25, 46.0 - 0.25]
	 p = alog10([1.2d-4, 7.8d-6, 7.0d-7, 1.8d-8, 2.0d-9, 2.9d-10])
	 d = [0.13, 0.13, 0.10, 0.30, 0.42, 0.42]

	L_solar = alog10(3.9) + 33.0
	L_HX = (L - L_solar) 
	PHI_HX  = p
	DPHI_HX = d

	endelse
end
