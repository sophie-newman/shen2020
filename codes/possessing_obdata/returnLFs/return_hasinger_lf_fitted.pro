;; Function to return the analytical Hasinger et al. (2005) luminosity function 
;;   for a list of soft XR luminosities L0_list (in SOLAR luminosities) at redshift z
;;   (for an Omega_M = 0.3, Omega_Lambda = 0.7 cosmology)
;;
function return_Hasinger_LF_fitted, L0_list, z
	L0_x 	= 10^(33.0d0) * 3.9 * 10^(L0_list)
	GAMMA_1 = 0.87
	GAMMA_2 = 2.57	
	L_X_star_0 = 10^(43.94d0)	;; (converts to solar luminosities)
	
	A_44  = 2.62d-7						;; Normalization - in AGNs/Mpc^3/log(L)
	
	p1_44 = 4.7
	p2_44 = -1.5
	beta_1 = 0.7
	beta_2 = 0.6
	p1_L  = p1_44 + beta_1 * (alog10(L0_x)-44.0)
	p2_L  = p2_44 + beta_2 * (alog10(L0_x)-44.0)

	zc_44 = 1.42
	L_c   = 10^(44.67d0)
	alpha = 0.21
	zc_0  = zc_44 / (((1.0d44)/(L_c))^alpha)

	PHI_44   = 1./((1.0d44/L_X_star_0)^(GAMMA_1) + (1.0d44/L_X_star_0)^(GAMMA_2))
	ed_44_zc = (1.+zc_44)^p1_44
	
	A0_44 = A_44 / PHI_44 ;;*ed_44_zc

	PHI_0 = 1./((L0_x/L_X_star_0)^(GAMMA_1) + (L0_x/L_X_star_0)^(GAMMA_2))
	PHI_0 = PHI_0 * A0_44
	z_c_1 = zc_0*((L0_x/L_c)^(alpha))
	z_c_2 = 0.0*L0_x + zc_0
	z_c   = z_c_1*(L0_x LE L_c) + z_c_2*(L0_x GT L_c)
	e_d_1 = 0.0*L0_x + (1.+z)^p1_L
	e_d_2 = ((1.+z_c)^(p1_L)) * (((1.+z)/(1.+z_c))^(p2_L))
	e_d   = e_d_1*(z LE z_c) + e_d_2*(z GT z_c)
	PHI   = PHI_0 * e_d

	;; want to correct for the type-2 fraction, since this is only a Type-1 
	;;  luminosity function: do so following the data in Hasinger et al. 2004
	;;  as compiled in Simpson et al. 2005
	;;
		L0 = 10^(DOUBLE(42.37))/0.015/3.9d33 * 0.606	;; our spectrum HX->SX
		;L0 = 10^(DOUBLE(42.37))/0.015/3.9d33 * 1.606	;; our spectrum HX->SX
		xi = 0.23
		f1 = 1. - ((1.+3.*((10^(DOUBLE(L0_x))/L0)^(1.-2.*xi)))^(-0.5))
		PHI = PHI/f1 * 1.3

	return, alog10(PHI)
end
