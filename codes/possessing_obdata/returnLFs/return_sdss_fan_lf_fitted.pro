;; Function to return the analytical Fan et al. (2001) luminosity function 
;;   for a list of B-band luminosities L0_list (in SOLAR luminosities) at redshift z
;;   (for an Omega_M = 0.3, Omega_Lambd = 0.7 cosmology)
;;
function return_SDSS_FAN_LF_fitted, L0_list, z
	forward_function M_sun_BB
	
	PHI_STAR = 7.2d-8
	alpha    = 0.75
	beta     = -2.58

	M_B    = M_sun_BB(0) - 2.5*(L0_list)
	M_1450 = M_B + 0.83		
		;; note this is a different convention for M_1450 than other places, 
		;;   but is directly from the paper (& used therein to define M_1450)
	PHI    = 2.5*PHI_STAR*10^(-0.4*(M_1450 + 26.0 - alpha*(z-3))*(beta+1.0))

	return, alog10(PHI)
end
