;; Function to return the analytical Miyaji et al. (2000) luminosity function 
;;   for a list of soft XR luminosities L0_list (in SOLAR luminosities) at redshift z
;;   (for an Omega_M = 0.3, Omega_Lambda = 0.7 cosmology)
;;
function return_Miyaji_LF_fitted, L0_list, z
	OMEGA_MATTER = 0.3	;; Cosmology - not important
	OMEGA_LAMBDA = 0.7  
	h_50 = 7.0/5.0		;; h in 50 km/s/Mpc
	;N_QSO = ~400?		;; Number used in sample ??
	P_KS  = 0.9			;; K-S probability of this fit

	GAMMA_1 = 0.66
	GAMMA_2 = 2.19
	
	L_X_star_0 = 7.14286e9		;; (converts to solar luminosities)
	
	A_0 = 4.41784e-6			;; Normalization - in AGNs/Mpc^3/log(L)
	
	p1 = 5.3		;; Account for redshift evolution
	p2 = 0.0
	;L_a = 10^(44.4)
	L_a = 6.27972e10
	z_cut = 1.58
	alpha = 2.6

	if (z GT z_cut) then begin
		e_z = 0.0*L0_list + (1.0 + z_cut)^p1 * ((1.0 + z)/(1.0+z_cut))^p2
	endif else begin
		e_z = 0.0*L0_list + (1.0 + z)^p1
		zl  = where(10^(L0_list) LT L_a, count)
		if (count GT 0) then begin 
			q  = p1 - alpha*alog10(L_a/10^(L0_list[zl]))
			ql = where(q LT 0.0, count_q)
			if (count_q GT 0) then q[ql] = 0.0
			e_z[zl] = (1.0 + z)^q
		endif
	endelse
	
	x = 10^(L0_list) / L_X_star_0
	PHI = A_0 * e_z / (x^GAMMA_1 + x^GAMMA_2)

	return, alog10(PHI)
end
