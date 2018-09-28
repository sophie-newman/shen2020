;;;
;;; Hunt et al. (possibly being updated with new data)
;;;   Measuring M_1450 (converted to M_B following DR3 data)
;;;
pro load_Hunt_LF_data, L_BB, PHI_BB, DPHI_BB, z
	forward_function M_sun_BB, lum_correct_cosmology, phi_correct_cosmology
	if ((z LE 2.5) OR (z GE 3.5)) then begin
		L_BB = 0.
		PHI_BB = 0.
		DPHI_BB = 0.
	endif else begin

	M_1450 = [ -20.4,  -21.4,  -22.4,  -23.4,  -24.4,  -25.4] $
				+ 5.*alog10(7./5.) - 2.5*alog10(lum_correct_cosmology(3.0)) 
					;; correct from h=0.5, Om=1, Ol=0
	PHI    = [1.2d-6, 2.8d-6, 3.5d-7, 1.5d-6, 5.0d-7, 9.1d-7] $
				* 2.5	* (7./5.)^3 * phi_correct_cosmology(3.0)
	DPHI   = [   0.4,   0.25,    0.7,    0.2,    0.5,    0.3]
	M_BB   = M_1450 + 1.75
	M_BB   = M_1450 - 0.776
	L_BB   = (0.4*(M_sun_BB(0)-M_BB))
	PHI_BB = alog10(PHI)
	DPHI_BB= DPHI
	
	endelse
end
