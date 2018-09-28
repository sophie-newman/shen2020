;;;
;;; Jiang et al. 2006 (low-L SDSS quasars)
;;;
pro load_Jiang_LF_data, L_BB, PHI_BB, DPHI_BB, z
	forward_function M_sun_BB

	if ((z LT 0.5) OR (z GT 3.6)) then begin
		L_BB = 0.
		PHI_BB = 0.
		DPHI_BB = 0.
	endif else begin
	if ((z GE 0.5) AND (z LT 1.0)) then WHICH_BLOCK = 0
	if ((z GE 1.0) AND (z LT 1.5)) then WHICH_BLOCK = 1
	if ((z GE 1.5) AND (z LT 2.0)) then WHICH_BLOCK = 2
	if ((z GE 2.0) AND (z LT 2.5)) then WHICH_BLOCK = 3
	if ((z GE 2.5) AND (z LT 3.0)) then WHICH_BLOCK = 4
	if ((z GE 3.0) AND (z LE 3.6)) then WHICH_BLOCK = 5

	if (WHICH_BLOCK EQ 0) then begin
	  M_g = [-21.25, -21.75, -22.50, -23.25, -24.25]
	  phi = [2.3d-6, 2.3d-6, 1.4d-6, 1.0d-6, 4.3d-7]
	  dphi= [0.3000, 0.3000, 0.3000, 0.3000, 0.3500]
	endif
	if (WHICH_BLOCK EQ 1) then begin
	  M_g = [-22.00, -22.50, -23.00, -23.50, -24.50]
	  phi = [3.5d-6, 3.0d-6, 2.3d-6, 2.2d-6, 7.8d-7]
	  dphi= [0.2500, 0.2400, 0.2400, 0.2300, 0.2400]
	endif
	if (WHICH_BLOCK EQ 2) then begin
	  M_g = [-22.60, -23.10, -23.60, -24.20, -24.75, -25.20]
	  phi = [4.0d-6, 3.0d-6, 2.9d-6, 1.7d-6, 2.25d-6, 1.05d-6]
	  dphi= [0.2000, 0.1800, 0.1800, 0.1800, 0.18000, 0.20000]
	endif
	if (WHICH_BLOCK EQ 3) then begin
	  M_g = [-23.00, -23.80, -24.15, -24.50, -25.40]
	  phi = [3.3d-6, 1.8d-6, 2.0d-6, 1.0d-6, 6.0d-7]
	  dphi= [0.4300, 0.3200, 0.3000, 0.3000, 0.4000]
	endif
	if (WHICH_BLOCK EQ 4) then begin
	  M_g = [-23.35, -23.70, -24.60, -26.00]
	  phi = [3.0d-6, 1.4d-6, 8.1d-7, 3.0d-7]
	  dphi= [0.4250, 0.3680, 0.3580, 0.4840]
	endif
	if (WHICH_BLOCK EQ 5) then begin
	  M_g = [-23.90, -24.60, -25.20]
	  phi = [1.0d-6, 7.9d-7, 4.3d-7]
	  dphi= [0.3520, 0.3500, 0.4820]
	endif
	phi = alog10(phi*2.5) ;; convert to per log(L)
	dphi= dphi/2. ;; numbers above for full +/- 1sigma range
	M_B = M_g + 0.05 ;; ;; bJ = g + 0.05 from Fukugita et al. (R06 DR3 ppr)

	L_BB   = (0.4*(M_sun_BB(0)-M_B))
	PHI_BB = PHI
	DPHI_BB= DPHI
	
	endelse
end

;;
;; here are the fit params including this data in the fitting procedure :: 
;;   although the fit is technically changed, looking at it, it's clear that 
;;   it's more or less *exactly* the same (like 0.01 dex change at most) at 
;;   all luminosities of interest for our fitting/comparisons, so 
;;   we know adding this data won't change our results
;; {P0=-4.80530; P1=13.0073; P2=0.629385; P3=-11.7977; P4=-14.3867; 
;;  P5=0.406049; P6=-0.649790; P7=2.17013;  P8=1.46220;  P9=-0.800668; 
;;  P10=0.;P11=0.;P12=0.;P13=0.;P14=0.;}
;;


pro tester, z_key
 z = 0.
 if (z_key EQ 0) then z = 0.75
 if (z_key EQ 1) then z = 1.25
 if (z_key EQ 2) then z = 1.75
 if (z_key EQ 3) then z = 2.25
 if (z_key EQ 4) then z = 2.75
 if (z_key EQ 5) then z = 3.30
 
 l = 8. + 0.1*findgen(81)
 f = alog10(qlf_calculator(l,z,/BB,/FULL))
 
 m0 = 4.77
 plot,m0-2.5*l,f,xstyle=1,xrange=[-20.,-26.1],yrange=[-6.5,-4.5],ystyle=1

 f = alog10(qlf_calculator(l,z,/BB,/PLE))
 oplot,m0-2.5*l,f,color=250

 load_jiang_lf_data,x,y,dy,z
 plotsym,0,/fill
 oploterror,m0-2.5*x,y,dy,psym=8
 
 f=load_cosmological_parameters(SET_OMEGA_MATTER=0.3,SET_OMEGA_LAMBDA=0.7,SET_HUBBLE=0.7)
 M_g = -19.50 - 1.08*7.50*frac_lookback_time(z) + 0.05
 L0  = 10^(-0.4*(M_g-m0))
 phi = 1.84d-6 * 2.5 / ((10^(l)/L0)^(-(-1.25+1.)) + (10^(l)/L0)^(-(-3.25+1.)))
 if (z GT 2.0) then $
 phi = 1.02d-6 * 2.5 / ((10^(l)/L0)^(-(-1.55+1.)) + (10^(l)/L0)^(-(-3.25+1.))) * 10^(-0.45*(z-2.0))
 oplot,m0-2.5*l,alog10(phi),linestyle=5,thick=4.
 
 phi_int = INTERPOL(alog10(phi),l,x)
 phi_renorm = TOTAL((phi_int-y)/dy/dy)/TOTAL(1./dy/dy)
 !p.color=210
 oploterror,m0-2.5*x,y+MEAN(phi_int)-MEAN(y),dy,psym=8
 ;oploterror,m0-2.5*x,y+phi_renorm,dy,psym=8

 !p.color=250
 load_2slaq_lf_data,x,y,dy,z
 plotsym,0,/fill
 oploterror,m0-2.5*x,y,dy,psym=8

 !p.color=80
 load_combo17_lf_data,x,y,dy,z
 plotsym,0,/fill
 oploterror,m0-2.5*x,y,dy,psym=8
 
 !p.color=150
 load_sdss_dr3_lf_data,x,y,dy,z
 plotsym,0,/fill
 oploterror,m0-2.5*x,y,dy,psym=8

 !p.color=0
end


