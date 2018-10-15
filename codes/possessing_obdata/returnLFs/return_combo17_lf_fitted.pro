;; Function to return the analytical Wolf et al. (2003) luminosity function 
;;   for a list of B-band luminosities L0_list (in SOLAR luminosities) at redshift z
;;   (for an Omega_M = 0.3, Omega_Lambd = 0.7 cosmology)
;;
function return_COMBO17_LF_fitted, L0_list, z
	forward_function M_sun_BB
	
	A0 = -5.620
	A1 = 0.1845
	A2 = -0.02652
	M0_145 = -25.0
		;; convert to a B-band M0
		M0_B = M0_145 + 1.75 +5.*alog10(7./6.5)
	B1 = 1.3455
	B2 = -80.845
	B3 = 127.32
	C1 = 0.3599
	C2 = -15.574
	
	zeta = alog10((1.+z)/(1.+2.))
	Mz_145 = M0_145 + B1*zeta + B2*zeta*zeta + B3*zeta*zeta*zeta

	Mz_B   = M0_B   + B1*zeta + B2*zeta*zeta + B3*zeta*zeta*zeta
	M_B_t  = M_sun_BB(0) - 2.5*(L0_list)
	mu     = M_B_t - Mz_B
	logphi = A0 + A1*mu + A2*mu*mu		;; PLE


	mu     = M_B_t - Mz_B
	logphi = A0 + A1*mu + A2*mu*mu + C1*zeta + C2*zeta*zeta	;; PDE

	PHI    = logphi + alog10(2.5) + 3.*alog10(7./6.5) 
	return, PHI
end
