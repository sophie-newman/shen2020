function return_Matute_LF_fitted, L0_list, z	
	;; given L in L_solar at 15 microns
	h_corr = 7./7.5

	;;; TYPE-1 OBJECTS
	;;log_L_15_over_L_R = 0.23
	;;L = (L0_list) - alog10(17.0) + log_L_15_over_L_R 
	;;L = (L0_list) - alog10(23.0) + log_L_15_over_L_R 
	L = (L0_list)
	phi_star = 10^(-6.26) * 2.5 
	L15_0    = 10^(44.78 - (alog10(3.9)+33.0)) * h_corr^2
		;;if (z LT 0.2) then L15_0 = L15_0 * 10^(+0.25)
	alpha_1  = 0.69
	alpha_2  = 5.92
	alpha_3  = 0.64
	beta     = 2.76
	k_L      = 2.87
	z_cut    = 2.00
	alpha    = alpha_1 * EXP(-alpha_2*z) + alpha_3	
	L15_ast = L15_0 * ((1. + z)^k_L)
	if (z GT z_cut) then L15_ast  = L15_0 * ((1. + z_cut)^k_L)
	x = (10^L)/L15_ast
	phi = phi_star / (x^alpha + x^beta) * h_corr^3


	;;; TYPE-2 OBJECTS
	;;log_L_15_over_L_R = 0.23
	;;L_R  = (L0_list) - alog10(17.0)
	;;L = (L_R - 5.02) / (1.0 - 0.47)
	;;L = (L0_list) - alog10(17.0) + log_L_15_over_L_R 
	L = (L0_list)
	phi_star = 10^(-4.69) * 2.5 
	L15_0    = 10^(44.11 - (alog10(3.9)+33.0)) * h_corr^2
	alpha_1  = 0.90
	alpha_2  = 6.00
	alpha_3  = 0.00
	beta     = 2.66
	k_L      = 2.04
	z_cut    = 2.00
	alpha    = alpha_1 * EXP(-alpha_2*z) + alpha_3	
	L15_ast = L15_0 * ((1. + z)^k_L)
	if (z GT z_cut) then L15_ast  = L15_0 * ((1. + z_cut)^k_L)
	x = (10^L)/L15_ast
	if (z GT 0.25) then phi = phi + phi_star / (x^alpha + x^beta) * h_corr^3

	
	;;if (z LT 0.2) then begin phi = phi * 10^(-0.2)
	return, alog10(phi) 
end
