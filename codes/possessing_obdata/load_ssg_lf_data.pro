;;; SSG -- broad range of redshift -- rest frame B-band
;;;
pro load_SSG_LF_data, L_BB, PHI_BB, DPHI_BB, z
	forward_function M_sun_BB, phi_correct_cosmology, lum_correct_cosmology
	if ((z LT 2.75) OR (z GT 4.75)) then begin
		L_BB = 0.
		PHI_BB = 0.
		DPHI_BB = 0.
	endif else begin
			
	pc = alog10(2.5 * phi_correct_cosmology(3.75)*((7./5.)^3))
	M_B_o = [-25.8,-26.3,-26.8,-27.5]+ 5.*alog10(7./5.) $
		- 2.5*alog10(lum_correct_cosmology(3.75)) ;; correct from h=0.5, Om=1, Ol=0
	phi_o = [-6.90, -7.30, -7.450, -8.20]  + pc
	phi_p = [-6.80, -7.20, -7.375, -8.05]  + pc
	phi_m = [-7.05, -7.40, -7.550, -8.40]  + pc
	dphi_o= 0.5*(phi_p-phi_m)

	L_BB   = (0.4*(M_sun_BB(0)-M_B_o))
	PHI_BB = phi_o
	DPHI_BB= dphi_o
	
	endelse
end
