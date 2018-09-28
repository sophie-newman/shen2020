pro load_Silverman_HX_LF_data, L_HX, PHI_HX, DPHI_HX, z
	;; determine which redshift interval its in
	WHICH_BLOCK = 8
	if ((z GE 0.2) AND (z LT 0.5)) then WHICH_BLOCK = 1
	if ((z GE 1.5) AND (z LT 2.0)) then WHICH_BLOCK = 4
	if ((z GE 3.0) AND (z LT 4.0)) then WHICH_BLOCK = 6
	
	if (WHICH_BLOCK EQ 8) then begin
		L_HX = [0.]
		PHI_HX = [0.]
		DPHI_HX = [0.]
	endif else begin
	if (WHICH_BLOCK EQ 1) then begin
		L_HX = [42.25,42.75,43.25,43.75,44.25,44.75]
		P_HX = [2.3d-4,1.7d-4,4.3d-5,2.3d-5,7.0d-6,2.0d-7]
		D_HX = [0.115,  0.081, 0.145, 0.157, 0.146, 0.222]
	endif
	if (WHICH_BLOCK EQ 2) then begin
		L_HX = [  42.5,  43.5,  44.5]
		P_HX = [2.2d-4,6.0d-5,2.5d-6]
		D_HX = [ 0.125, 0.125, 0.215]
	endif
	if (WHICH_BLOCK EQ 3) then begin
		L_HX = [  43.5,  44.5]
		P_HX = [5.0d-5,7.5d-6]
		D_HX = [ 0.130, 0.204]
	endif
	if (WHICH_BLOCK EQ 4) then begin
		L_HX = [43.25,43.75,44.25,44.75,45.25]
		P_HX = [2.1d-5,4.0d-5,2.8d-5,7.0d-6,1.2d-6]
		D_HX = [ 0.226, 0.138, 0.101, 0.109, 0.221]
	endif
	if (WHICH_BLOCK EQ 5) then begin
		L_HX = [  43.5,  44.5]
		P_HX = [6.8d-6,2.8d-6]
		D_HX = [ 0.167, 0.196]
	endif
	if (WHICH_BLOCK EQ 6) then begin
		L_HX = [43.75,44.25,44.75,45.25]
		P_HX = [9.0d-6,3.0d-6,1.6d-6,8.0d-8]
		D_HX = [ 0.222, 0.204, 0.155, 0.574]
	endif
	if (WHICH_BLOCK EQ 7) then begin
		L_HX = [  43.5,  0.5*(44.5+45.5)]	;; the latter is soft X-ray selected
		P_HX = [1.2d-5,2.0d-7]
		D_HX = [ 0.368, 0.398]
	endif
	L_solar = alog10(4.0) + 33.0
	L_HX = (L_HX - L_solar) + alog10(1.19) ;; convert the Barger 2-8kev to 2-10kev (gamma=1.8)
	PHI_HX  = alog10(P_HX)
	DPHI_HX = D_HX + 0.1
	endelse
end
