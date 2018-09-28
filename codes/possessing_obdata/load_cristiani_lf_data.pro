;;; Cristiani 2004 faint at z=3.5-5.2
;;;
pro load_Cristiani_LF_data, L_BB, PHI_BB, DPHI_BB, z
	forward_function M_sun_BB
	
	if ((z LE 4.0) OR (z GT 5.2)) then begin
		L_BB = 0.
		PHI_BB = 0.
		DPHI_BB = 0.
	endif else begin
		M_1450_o = [-22.3]			; z = 4.0 - 5.2
		phi_o    = [10^(-5.9)]
		phi_p    = [10^(-5.55)]
		phi_m    = [10^(-6.4)]
	M_B_o = M_1450_o - 0.83	
		;; note this is a different convention for M_1450 than other places, 
		;;   but is directly from the paper (& used therein to define M_1450)
	dphi  = 0.5*(alog10(phi_p) - alog10(phi_m))
	L_BB  = (0.4*(M_sun_BB(0)-M_B_o))
	PHI_BB= alog10(2.5*phi_o)
	DPHI_BB=dphi	
	endelse
end
