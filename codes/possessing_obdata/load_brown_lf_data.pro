;;;
;;;  Brown et al. 2005; 1<z<5 combined data set - but binned only given at z=2
;;;		8 micron data
;;;
;;;  potential issue -- the Brown et al. sample is only Type 1's. which even at the brightest 
;;;     luminosities will be affected sign. by obsc. since lower-lum heavily obscured 
;;;     sources can be up-scattered up here (see e.g. the matute et al. LFs)
;;;
pro load_Brown_LF_data, L_IR, PHI_IR, DPHI_IR, z
	forward_function M_sun_BB
	L_solar = 3.9d33
	if ((z LT 1.5) OR (z GT 2.5)) then begin
		L_IR = 0.
		PHI_IR = 0.
		DPHI_IR = 0.
	endif else begin
		M_8  = [-27.2, -27.6, -28.0, -28.4, -28.8, -29.2, -29.6, -30.0, -30.4, -31.2]
		P_8  = [4.3d-6,2.4d-6,3.0d-6,1.2d-6,6.1d-7,2.6d-7,2.4d-7,9.2d-8,4.7d-8,2.3d-8]
		D_8  = [ 0.21,  0.18,  0.10,  0.15,  0.15,  0.18,  0.18,  0.35,  0.40,  0.45]

	;; use below if converting to another band
	;;M_B = M_8 + 3.08
	;;L_B = 10^(DOUBLE(0.4*(M_sun_BB(0) - M_B)))		
	;;L_HX   = alog10(L_B) + alog10(10.0)		

	L_8 = 10^(0.4*(-28.5 - M_8)) * 4.11d45	;; conversion given therein
		L_8 = L_8 / L_solar
		L_15_over_L_8 = 10^(+0.060001373)	;; from Hatziminaoglou model spectrum
			;L_15_over_L_8 = 10^(+0.0263996)	;; from Hatziminaoglou model spectrum
		L_15 = L_15_over_L_8 * L_8
	L_IR = alog10(L_15)
	PHI_IR = alog10(P_8) + alog10(2.5) ;;+ alog10(1./(1.-0.55))
							;; correction for mean Type-2 fraction (NOTE: not 
							;;   obscured in the IR, but not optical qso so 
							;;   no redshift data in this limited sample
	DPHI_IR= D_8
	endelse
end
