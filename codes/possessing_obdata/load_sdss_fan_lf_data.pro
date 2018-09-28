;;; Fan et al. -- effectively supplanted by DR3 LF at z<5.0
;;;		Measures M_1450, but with a different convention for it than others
;;;
pro load_SDSS_FAN_LF_data, L_BB, PHI_BB, DPHI_BB, z
	forward_function M_sun_BB
	
	if (((z LT 3.6) OR (z GT 6.2)) OR (((z GT 5.0) AND (z LT 5.75)))) then begin
		L_BB = 0.
		PHI_BB = 0.
		DPHI_BB = 0.
	endif else begin
	Mcorr = +5.*alog10(7./6.5)
	Pcorr = (7./6.5)^3				;; cancelled to ~1% by the change in Omega_Lambda (*small*)
	Mcorr = 0.
	Pcorr = 1.
	if ((z GE 3.6) AND (z LE 3.9)) then begin
			M_1450_o = [-25.83,-26.30,-26.80] + Mcorr	; z = 3.6-3.9
			phi_o    = [6.0e-8,2.5e-8,7.7e-9] * Pcorr
			phi_p    = [9.2e-8,3.5e-8,1.1e-8] * Pcorr
			phi_m    = [2.4e-8,1.5e-8,4.4e-9] * Pcorr
	endif
	if ((z GT 3.9) AND (z LT 4.4)) then begin
			M_1450_o = [-26.00,-26.45,-27.18] + Mcorr	; z = 3.9-4.4
			phi_o    = [2.5e-8,1.1e-8,6.0e-9] * Pcorr
			phi_p    = [4.0e-8,1.7e-8,9.0e-9] * Pcorr
			phi_m    = [1.0e-8,8.0e-9,3.2e-9] * Pcorr
	endif
	if ((z GE 4.4) AND (z LE 5.0)) then begin
			M_1450_o = [-26.45,-27.08] + Mcorr		; z = 4.4-5.0
			phi_o    = [5.2e-9,2.2e-9] * Pcorr
			phi_p    = [8.0e-9,3.5e-9] * Pcorr
			phi_m    = [2.6e-9,9.0e-10] * Pcorr
	endif
	if ((z GE 5.75) AND (z LE 6.2)) then begin
			M_1450_o = [-26.7,-27.7]			; z = 5.several - 6.highest
			phistar  = 3.3d-9
			dphistar = 3.8d-9
			beta     = -3.2
			dbeta    = 0.8
			phi_o    = phistar*10^(-0.4*(M_1450_o+26.0)*(beta+1.0))
			phi_p    = 0.0*phi_o + (3.3+3.8)*1.0d-9
			phi_m    = 0.0*phi_o + (3.3-1.6)*1.0d-9


			;; updated with Fan 2004 data -- ''anchor for 03-04 data at about 
			;;   M_1450 = -27.12
			M_1450_o = [ -26.11, -27.12, -27.88 ] 
				;; min M_1450 of a z~6 qso, anchor point, and max M_1450
			phi_26 = 3.3d-9
			phip_26 = 3.8d-9
			phim_26 = 1.6d-9
			slp_26 = -3.2
			slp_p  = -2.5
			slp_m  = -4.0
				
			phi_anchor = 0.34d9
			phimin_anchor = 0.32d9	;; letting norm -> +1 sigma, slope -> -1 sigma
			phimax_anchor = 0.36d9	;; letting norm -> -1 sigma, slope -> +1 sigma
			
			p_0 = [3.3 , 0.34, 0.073]	;; including amplified Lyalpha emitted
			p_p = [7.1 , 0.36, 0.127]
			p_m = [1.7 , 0.32, 0.039]
			
			p_0 = [2.18 , 0.34, 0.096]	;; excluding amplified Lyalpha emitted
			p_p = [5.26 , 0.36, 0.148]		;; also more consistent with '03 results
			p_m = [1.00 , 0.32, 0.058]
			
	;; mockup to get fit leverage ::
			M_1450_o = [ -26.11, -26.75, -27.12, -27.88 ] 
			p_0 = [2.18, 0.53, 0.34, 0.096]		;; excluding amplified Lyalpha emitted
			p_p = [5.26, 0.65, 0.36, 0.148]		;; also more consistent with '03 results
			p_m = [1.00, 0.44, 0.32, 0.058]

			M_1450_o = [ -26.11, -26.75, -27.12, -27.88 ] 
			p_0 = [2.18, 0.53, 0.34, 0.096]			;; excluding amplified Lyalpha emitted
			p_p = [5.26, 0.65, 0.36, 0.148]			;; also more consistent with '03 results
			p_m = [1.00, 0.44, 0.32, 0.058]
			
			;;p_0 = p_0 + [0.089,-0.056,-0.058,-0.014]
			
			phi_o = p_0 * 1.0d-9
			phi_p = p_p * 1.0d-9 * 1.2	;; ends up giving a min ~0.1 dex scatter -- approp.
			phi_m = p_m * 1.0d-9 * 0.8	;;   for the 19-object sample - GOOD		
			
	endif
	
	
	M_B_o = M_1450_o - 0.83	
		;; note this is a slightly different convention for M_1450 than other places, 
		;;   but is directly from the paper (& used therein to define M_1450)
	dphi  = 0.5*(alog10(phi_p) - alog10(phi_m))
	L_BB  = (0.4*(M_sun_BB(0)-M_B_o))
	PHI_BB= alog10(2.5*phi_o)
	DPHI_BB=dphi	
	endelse
end
