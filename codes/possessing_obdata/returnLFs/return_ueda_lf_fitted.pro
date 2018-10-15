function return_Ueda_LF_fitted, L0_list, z
	OMEGA_MATTER = 0.3	;; Cosmology - not important
	OMEGA_LAMBDA = 0.7  
	N_QSO = 141			;; Number used in sample ??
	P_KS  = 0.9			;; K-S probability of this fit

	GAMMA_1 = 0.86
	GAMMA_2 = 2.23
	
	;L_X_star_0 = 10^(43.94)
	L_X_star_0 = 2.1774e10		;; (converts to solar luminosities)
	
	A_0 = 5.04e-6	;; Normalization - in AGNs/Mpc^3/log(L)
	
	p1 = 4.23		;; Account for redshift evolution
	p2 = -1.5
	;L_a = 10^(44.6)
	L_a = 1.0e11
	z_cut_star = 1.9
	alpha = 0.335

	z_cut = 0.0*L0_list + z_cut_star
	zcut_ldde = where(10^(L0_list) LT L_a, count)
	if (count GT 0) then z_cut[zcut_ldde] = z_cut_star * (10^(L0_list[zcut_ldde])/L_a)^alpha
	e_z = 0.0*L0_list + (1.0 + z)^p1
	zl = where(z_cut LE z, count)
	if (count GT 0) then e_z[zl] = ((1.+z_cut[zl])^(p1)) * ((1.+z)/(1.+z_cut[zl]))^p2
	
	x = 10^(L0_list) / L_X_star_0
	PHI = A_0 * e_z / (x^GAMMA_1 + x^GAMMA_2)

	return, alog10(PHI)
end
