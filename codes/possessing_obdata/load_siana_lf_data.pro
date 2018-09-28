;;;
;;; Siana et al. 2006 (IR + optical selection (SWIRE)) at z~3
;;;   Measuring M_1450 -- convert to B-band following DR3 corrections
;;;
pro load_Siana_LF_data, L_BB, PHI_BB, DPHI_BB, z
	forward_function M_sun_BB
	if ((z LE 2.9) OR (z GE 3.4)) then begin
		L_BB = 0.
		PHI_BB = 0.
		DPHI_BB = 0.
	endif else begin

	M_1450 = [-23.75, -24.25, -24.75, -25.25, -25.75, -26.25]
	PHI    = [1.0d-6, 8.0d-7, 5.5d-7, 3.0d-7, 1.2d-7, 1.8d-7] * 2.5
	DPHI   = [ 0.097,  0.056,  0.105,  0.125,  0.176,  0.143]	;; way too optimistic
	DPHI   = [ 0.097,  0.125,  0.155,  0.176,  0.269,  0.234]	;; use lower errors - much better log-errors
	M_BB   = M_1450 - 0.776
	L_BB   = (0.4*(M_sun_BB(0)-M_BB))
	PHI_BB = alog10(PHI)
	DPHI_BB= DPHI
	
	endelse
end
