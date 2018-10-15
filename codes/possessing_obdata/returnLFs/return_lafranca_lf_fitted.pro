function return_LaFranca_LF_fitted, L0_list, z	
	L_X_star_0 = 10^(DOUBLE(44.25 - alog10(3.9) - 33.0))		;; (converts to solar luminosities)
	GAMMA_1 = 1.01
	GAMMA_2 = 2.38
	
	A_0 = 1.21d-6	;; Normalization - in AGNs/Mpc^3/log(L)	
	p1 = 4.62		;; Account for redshift evolution
	p2 = -1.15
	L_a = 10^(DOUBLE(45.74 - alog10(3.9) - 33.0))
	z_cut_star = 2.49
	alpha = 0.20

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
