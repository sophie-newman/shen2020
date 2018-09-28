;;; Hao et al. local SDSS -- choose their OIII measure, as its least 
;;;   contaminated by star formation -- optical technically, but really an isotropic measure
;;;   that seems to be correlated with L(2-10kev), and even better correlated with harder 
;;;   X-rays (Mulchaey 1994, Heckman et al. 2005 ApJ 634, 161)
;;;
;; remember, should ignore obscuration when doing the narrow line fits, 
;;   as the Heckman study demonstrates
;;
pro load_Hao_LF_data, L_HX, PHI_HX, DPHI_HX, z
	forward_function M_sun_BB

	if (z GT 0.15) then begin
		L_HX = 0.
		PHI_HX = 0.
		DPHI_HX = 0.
	endif else begin
	
	LO3 = [2.0d5,  3.0d5,  4.4d5,  6.5d5,  9.5d5,  1.4d6,  2.1d6,  3.1d6,  4.8d6,  7.0d6, $
		   1.0d7,  1.5d7,  2.0d7,  3.1d7,  5.0d7,  7.0d7,  1.0d8,  1.6d8,  2.2d8,  3.5d8,  5.0d8]
	PO31= [1.1d-2, 5.0d-3, 2.5d-3, 1.4d-3, 1.3d-3, 5.1d-4, 3.0d-4, 2.6d-4, 2.8d-4, 1.9d-4, $
	       1.3d-4, 8.0d-5, 5.5d-5, 2.9d-5, 1.4d-5, 4.5d-6, 7.0d-6, 1.4d-6, 5.0d-8, 1.4d-6, 2.0d-8]
	PO32= [1.1d-2, 5.0d-3, 2.5d-3, 1.4d-3, $
	       0.5*1.3d-3, 0.5*5.1d-4, 0.5*3.0d-4, 0.5*2.6d-4, $
	       0.3*2.8d-4, 0.3*1.9d-4, 0.3*1.3d-4, 0.3*8.0d-5, 0.3*5.5d-5, $
	       0.5*2.9d-5, 0.5*1.4d-5, $
	       4.5d-6, 2.0d-6, 9.0d-8, 1.8d-6, 0.1d-8, 0.1d-8]	;; Type 2 contribution to each bin
	PO3 = alog10(PO31 + PO32)
	D00 = alog10(3./2.)
	DO3 = [   D00,    D00,    D00,    D00,    D00,    D00,    D00,    D00,    D00,    D00, $
	          D00,    D00,    D00,    D00,    0.2,    0.3,    0.3,   0.48,    0.4,   0.35,   0.8]
	h = 0.7	;; they adopt H0 = 100
	PO3 = PO3 + 3.*alog10(h)			
	LO3 = LO3/h/h

	;; alternatively, go with the *all H-alpha* LF, which is substantially higher, 
	;;   but may be significantly contaminated
	if (2 EQ 0) then begin
		L_star = 3.16d8 /h/h
		alpha  = -1.94
		phi_star = 9.17d-6 *h*h*h
		PO3 = alog10(alog(10.)*phi_star*((LO3/L_star)^(alpha+1.0))*exp(-(LO3/L_star)))	

		;; all O3
		L_star = 2.13d10 /h/h
		alpha  = -2.15
		phi_star = 6.19d-9 *h*h*h
		PO3 = alog10(alog(10.)*phi_star*((LO3/L_star)^(alpha+1.0))*exp(-(LO3/L_star)))	

		;; type-2 O3
		L_star = 4.70d8 /h/h
		alpha  = -2.14
		phi_star = 4.13d-7 *h*h*h
;		PO3 = alog10(alog(10.)*phi_star*((LO3/L_star)^(alpha+1.0))*exp(-(LO3/L_star)))	

		;; type-1 O3
		L_star = 5.17d6 /h/h
		alpha  = -1.59
		phi_star = 1.36d-4 *h*h*h
;		PO3 = alog10(10^(PO3) + alog(10.)*phi_star*((LO3/L_star)^(alpha+1.0))*exp(-(LO3/L_star)))	

	endif


	;; Mulchaey calibration considerably older, for 2-10; 
	;;   Heckman's more recent study has more (still only ~50), shows a problem -- 
	;;   2-10 significantly affected by obscuration (compton-thick) at the lowest L, 
	;;   and so very different ratio inferred from e.g. type 1 & type 2 ( & not clear 
	;;   that the type 1's are completely unaffected). also forces splitting the 
	;;   sample -> worse stats. But 3-20 seems clean - no type 1/2 dep't. Mean there 
	;;   is a 2.15 dex offset, then from our model spectrum we can correct to 
	;;   2-10 = (2/3) * (3-20); since we're interested in the intrinsic spectrum that 
	;;   should be ok?

	L_O3_over_L_2to10keV = 0.01 ;4 	;; Mulchaey calibration
	;L_O3_over_L_2to10keV = 10^(-2.00015) / (0.661) ;; very similar -- at this point, just splitting hairs
	;L_O3_over_L_2to10keV = 10^(-1.59)  ;; direct for unobscured type-1s, less certain (could be ~0.2 dex larger)
		;; Heckman calibration with our 2-10/3-20 intrinsic spectrum
	L_HX   = alog10( LO3 / L_O3_over_L_2to10keV )
	PHI_HX = PO3	;; assuming a linear proportionality for now
	DPHI_HX= DO3

	;; alternatively, adopt the correlations calibrated in Hao et al
	;;   -- ugh, altogether gives wierd results -- better off going with the direct calibration
	if (2 EQ 0) then begin
		M_i = (alog10(LO3) + 1.38)/(-0.470)
		M_B = M_i + 0.55
		MB0 = M_sun_BB(0)
		L_B = (0.4*(MB0-M_B))
		lbg = 6.+0.1*findgen(121)
		lbb = fitted_band_lum(lbg,/BB)
		lhx = fitted_band_lum(lbg,/HX)
		L_HX = INTERPOL(lhx,lbb,L_B)
		L_HX = 10^(DOUBLE(L_HX))
		dL_over_dLhx = fitted_band_lum(lbg,/HX,/JACOBIAN)
		dL_over_dLbb = fitted_band_lum(lbg,/BB,/JACOBIAN)
		dLalpha_over_dLbb = (0.470/0.4)	;; from M_i(L_Halpha)
		PHI_HX = PHI_HX + alog10(dLalpha_over_dLbb*dL_over_dLhx/dL_over_dLbb)
	endif

	endelse
end
