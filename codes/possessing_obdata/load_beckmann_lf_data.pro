;;; 
;;; Beckmann et al. 2006 local data
;;;
pro load_Beckmann_LF_data, L_HX, PHI_HX, DPHI_HX, z
	;; determine which redshift interval its in
	if (z GT 0.044) then begin
		L_HX = 0.
		PHI_HX = 0.
		DPHI_HX = 0.
	endif else begin
 ;; beckmann -- they measure 20-40 keV :: from our model spectrum, 
 ;;   L(20-40) = 1.15 * L(2-10)
	 l = [41.25, 41.75, 42.25, 42.75, 43.25, 43.75, 44.25, 44.75, 45.25] - alog10(1.1522883)
	 p = [-3.4,  -3.7,  -4.35, -4.7,   -5.25, -6.0,  -7.0, -8.2, -9.1]
	 d = [0.3,    0.18, 0.17, 0.14,     0.1,   0.11,  0.15, 0.25, 0.26]

	L_solar = alog10(3.9) + 33.0
	L_HX = (L - L_solar) 
	PHI_HX  = p
	DPHI_HX = d
 ;; still slightly low, but this sample is *really* local (z_bar = 0.022), so it's ok

	endelse
end
