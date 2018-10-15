;; Function to return the analytical Richards et al. (2005) luminosity function 
;;   for a list of B-band luminosities L0_list (in SOLAR luminosities) at redshift z
;;   (for an Omega_M = 0.3, Omega_Lambd = 0.7 cosmology)
;;
function return_2SLAQ_LF_fitted, L0_list, z
	forward_function M_sun_BB
	OMEGA_MATTER = 0.3	;; Cosmology - not important
	OMEGA_LAMBDA = 0.7  
	N_QSO = 6180		;; Number used in sample
	P_KS  = 0.075		;; K-S probability of this fit

	ALPHA = 3.28
	BETA  = 1.78

	M_B_star_0 = -22.68	 + 0.05 ;; Break magnitude at z=0, with b/g conversion
	L_B_star_0 = 10^((M_sun_BB(0)-M_B_star_0)/2.5)
	k1 =  1.37			;; Determines the polynomial evolution of the break
	k2 = -0.32			;;   magnitude, as below
	
	PHI_STAR = 5.96e-7	;; Normalization - in AGNs/Mpc^3/mag
	;; Convert this to AGNs/Mpc^3/log(L) for luminosity units, 
	;;   simply using dmag = -2.5 d(log(L)), but then get cancelling
	;;   factor of 2.5 from other side of the equation
	PHI_STAR_L = PHI_STAR * 2.5
	
	;; Account for redshift evolution
	;; M_B_star = M_B_star_0 - 2.5 * (k1*z + k2*z*z), so have 
	L_B_star = L_B_star_0 * 10^(k1*z + k2*z*z)

	x = 10^(L0_list) / L_B_star
	PHI = PHI_STAR_L / (x^(ALPHA-1.) + x^(BETA-1.))
	return, alog10(PHI)
end

