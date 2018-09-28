;;; KDC -- rest B-band
;;;
pro load_KDC_LF_data, L_BB, PHI_BB, DPHI_BB, z
	forward_function M_sun_BB, lum_correct_cosmology, phi_correct_cosmology
	if ((z LT 4.0) OR (z GT 4.5)) then begin
		L_BB = 0.
		PHI_BB = 0.
		DPHI_BB = 0.
	endif else begin
			
	M_B_o = [-27.0,-28.0]+ 5.*alog10(7./5.) - 2.5*alog10(lum_correct_cosmology(4.25))  ;; correct from h=0.5, Om=1, Ol=0
	phi_o  = [5.7e-9,5.5e-10] * 2.5 * (0.7/0.5)^3 * phi_correct_cosmology(4.25)
	dphi_o = [2.6e-9,3.6e-10] * 2.5 * (0.7/0.5)^3 * phi_correct_cosmology(4.25)

	L_BB   = (0.4*(M_sun_BB(0)-M_B_o))
	PHI_BB = alog10(PHI_o)
	DPHI_BB= dphi_o/phi_o/alog(10.)
	
	endelse
end
