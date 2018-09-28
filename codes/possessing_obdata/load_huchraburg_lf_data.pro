;;;
;;; Huchra & Burg local AGN luminosity function
;;;
pro load_HuchraBurg_LF_data, L_BB, PHI_BB, DPHI_BB, z
	forward_function M_sun_BB
	if (z GT 0.2) then begin
		L_BB = 0.
		PHI_BB = 0.
		DPHI_BB = 0.
	endif else begin

	M_BB   = [-17.75,-18.25,-18.75,-19.25,-19.75,-20.25,-20.75,-21.25,-21.75,-22.25]
	log_phi= [-3.83, -3.69, -4.08, -3.90, -4.32, -4.46, -4.65, -5.29, -6.45, -6.54 ]

	h = 0.7
	L_BB   = (0.4*(M_sun_BB(0)-M_BB)) - alog10(h*h)
	PHI_BB = log_phi + alog10(2.5) + alog10(h*h*h)
	DPHI_BB= 0.0*log_phi + 0.2
	
	endelse
end
