;;; 
;;; Sazonov & Revnivtsev al. 2004 local data
;;;
pro load_SazRev_LF_data, L_HX, PHI_HX, DPHI_HX, z
	;; determine which redshift interval its in
	if (z GT 0.15) then begin
		L_HX = 0.
		PHI_HX = 0.
		DPHI_HX = 0.
	endif else begin
		L_HX = [41.25,  41.75,  42.25,  42.75,  43.25,  43.75,  44.25,  44.75,  45.25,  45.75]
		P_HX = [1.0d-3, 2.9d-4, 1.0d-4, 3.5d-5, 7.3d-6, 1.5d-6, 1.7d-7, 2.2d-8, 4.8d-10,1.4d-10]
		D_HX = [ 0.4,    0.35,   0.30,   0.20,   0.20,   0.15,   0.20,   0.25,   0.35,   0.4]

	L_solar = alog10(3.9) + 33.0
	L_HX = (L_HX - L_solar) + alog10( 0.66676 * (0.75/0.7)^2 )
	;L_HX = (L_HX - L_solar) + alog10( (0.75/0.7)^2 ) ;+ alog10(1./1.25)
		;; convert the 2-30kev to 2-10kev (using our spectrum model in the X-rays)
		;;  and from H0 = 75 to our 70
	PHI_HX  = alog10(P_HX *1.43* (0.75/0.7)^3)	;; includes the completeness correction they estimate
	DPHI_HX = D_HX

 	;; Sazonov et al. 2006 (17-60 keV, 0.0 <~ z <~ 0.1-1.0 (ideal 0.15))
	 l = [ 40.25,  41.75,  42.25,  42.75,  43.25,  43.75,  44.25,  44.75,  45.25] - 2.*alog10(0.7/0.75) - alog10(2.0)
	L_HX = [L_HX,l-L_solar]
	 p = alog10([3.2d-2, 7.5d-4, 3.0d-4, 8.0d-5, 3.5d-5, 4.3d-6, 5.0d-7, 8.0d-8, 1.3d-9] * (0.7/0.75)^3)
	PHI_HX = [PHI_HX,p]
	 d = [0.4,    0.3,    0.24,   0.16,   0.09,   0.15,   0.17,   0.22,   0.39]
	DPHI_HX = [DPHI_HX,d]

	endelse
end

pro tester, z_key
 z = 0.0
 
 l = 8. + 0.1*findgen(81)
 f = alog10(qlf_calculator(l,z,/HX,/FULL))
 
 plot,l+alog10(3.9)+33.+alog10(0.5),f,xstyle=1,xrange=[39.9,47.],ystyle=1,yrange=[-10.5,-1.]

 load_SazRev_lf_data,x,y,dy,z
 plotsym,0,/fill
 oploterror,x+alog10(4.0)+33.,y,dy,psym=8

 !p.color=250
 load_ueda_lf_data,x,y,dy,z
 oploterror,x+alog10(3.9)+33.,y,dy,psym=8

 !p.color=210
 load_hao_lf_data,x,y,dy,z
 oploterror,x+alog10(3.9)+33.,y,dy,psym=8


 ;; shinozaki
 !p.color=150
 load_shinozaki_lf_data,x,y,dy,z
 oploterror,x+alog10(3.9)+33.,y,dy,psym=8


 ;; beckmann
 !p.color=50
 load_beckmann_lf_data,x,y,dy,z
 oploterror,x+alog10(3.9)+33.,y,dy,psym=8


 !p.color=0 
 l = 40. + 0.1*findgen(101)

 ;; shinozaki sample
 p = 3.68 / (10^((l-44.02)*0.93) + 10^((l-44.02)*(2.51))) * 1.0d-6
 oplot,l,alog10(p),thick=2.,color=80

 ;; U03 sample
 p = 2.64 / (10^((l-44.11)*0.93) + 10^((l-44.11)*(2.23))) * 1.0d-6
; oplot,l,alog10(p),thick=2.,color=250


 ;; sazrev update 
 l0 = 43.40 - 2.*alog10(0.7/0.75) - alog10(2.0)
 p0 = 3.55 * 1.0d-5 * (0.7/0.75)^3
 p = p0 / (10^((l-l0)*0.76) + 10^((l-l0)*(2.28))) 
 ;oplot,l,alog10(p),thick=2.,color=150

 ;; beckmann
 l0 = 43.38 - alog10(1.1522883)
 p0 = 0.7 * 1.0d-5 
 p = p0 / (10^((l-l0)*0.80) + 10^((l-l0)*(2.11))) 
 oplot,l,alog10(p),thick=2.,color=50


 !p.color=0
end
