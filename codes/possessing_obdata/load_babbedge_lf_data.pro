pro load_babbedge_LF_data ;;, L_IR, PHI_IR, DPHI_IR, z
	h = 0.7
	DLOGL = alog10(8.)/5.
	c = 3.0d8
	mm= 1.0d-6
	nu_norm = alog10(c/(15.0*mm))
	f_norm  = richards_return_spectrum(nu_norm,/BLUE)
		!P.COLOR=0
		plot,findgen(2),findgen(2),/nodata,xstyle=1,ystyle=1,xrange=[9.5,12.1],yrange=[-7.0,-4.5]
		plotsym,0,/fill

	;; 3.6 micron
	L = 9.55506 + DLOGL*findgen(7)
	P = [2.3d-5, 2.0d-5, 1.3d-5, 7.0d-6, 4.5d-6, 3.2d-6, 1.9d-6]
	dp= [0.18,     0.17,   0.20,   0.43,   0.45,   0.43,   0.52]
		L = L - 2.*alog10(h) + alog10(f_norm/richards_return_spectrum(alog10(c/(3.6*mm)),/BLUE))
		P = alog10(P) + 3.*alog10(h)
		!P.COLOR=80
			print, L
			print, P
			oploterror,L + 1.05,P,dp,PSYM=8		
	;; 4.5 micron
	L = 9.55506 + DLOGL/2. + DLOGL*findgen(5)
	P = [2.2d-5, 1.5d-5, 9.0d-6, 6.0d-6, 3.1d-6]
	dp= [0.12,     0.33,   0.43,   0.43,   0.55]
		L = L - 2.*alog10(h) + alog10(f_norm/richards_return_spectrum(alog10(c/(4.5*mm)),/BLUE))
		P = alog10(P) + 3.*alog10(h)
		!P.COLOR=110
			oploterror,L + 1.05,P,dp,PSYM=8		
	;; 5.8 micron
	L = 9.55506 + 3.*DLOGL + DLOGL*findgen(10)
	P = [9.0d-6, 7.8d-6, 9.0d-6, 7.8d-6, 7.8d-6, 7.8d-6, 6.5d-6,  5.5d-6, 4.0d-6, 2.9d-6]
	dp= [0.46,     0.24,   0.35,   0.40,   0.38,   0.18,   0.33,    0.43,   0.43,   0.55]
		L = L - 2.*alog10(h) + alog10(f_norm/richards_return_spectrum(alog10(c/(5.8*mm)),/BLUE))
		P = alog10(P) + 3.*alog10(h)
		!P.COLOR=150
			oploterror,L,P,dp,PSYM=8		
	;; 8.0 micron
	L = 9.55506 + 6.*DLOGL + DLOGL*findgen(5)
	P = [1.0d-5, 9.0d-6, 8.0d-6, 6.8d-6, 3.9d-6]
	dp= [0.30,     0.30,   0.60,   0.30,   0.62]
		L = L - 2.*alog10(h) + alog10(f_norm/richards_return_spectrum(alog10(c/(8.0*mm)),/BLUE))
		P = alog10(P) + 3.*alog10(h)
		!P.COLOR=210
			oploterror,L,P,dp,PSYM=8		
	;; 24.0 micron
	L = 9.55506 + 4.*DLOGL + DLOGL*findgen(6)
	P = [8.0d-6, 8.5d-6, 8.0d-6, 6.5d-6, 5.5d-6, 4.0d-6]
	dp= [0.60,     0.42,   0.30,   0.40,   0.40,   0.65]
		L = L - 2.*alog10(h) + alog10(f_norm/richards_return_spectrum(alog10(c/(24.*mm)),/BLUE))
		P = alog10(P) + 3.*alog10(h)
		!P.COLOR=250
			oploterror,L,P,dp,PSYM=8		


		!P.COLOR=0
	lambda_ang = 10^(1. + 0.05 * findgen((7.-2.)/0.05+1.)) ;; angstroms
		lambda = lambda_ang * 1.0d-10
	nu = c/lambda
	f = richards_return_spectrum(alog10(nu),/IRLUM)
		;; f = nu f_nu = lambda * f_lambda
		f_lam = f/lambda
		fnorm = INTERPOL(f_lam,lambda_ang,[3.1d4])
		f_lam = f_lam * (1.0d-4 / fnorm(0))
	plot,lambda_ang,f_lam,/xlog,/ylog,xstyle=1,xrange=[2.0d2,1.0d5],ystyle=1,yrange=[1.0d-4,1.0]
	plot,lambda_ang,f_lam,/xlog,/ylog,xstyle=1,xrange=[1.0d3,5.0d5],ystyle=1,yrange=[1.0d-4,1.0]
		oplot,lambda_ang,f_lam*50.,color=110
		micron = 1.0d4
		l0 = 3.6 * micron
		z_l = [0.5, 1.5, 2.5, 3.5]
		oplot,[l0,l0],[1.0d-10,1.0d10]
			for i=0, n_elements(z_l)-1 do begin
				corr = 1./(1.+z_l[i])
				oplot,[l0*corr,l0*corr],[1.0d-10,1.0d10],LINESTYLE=i+1
			endfor
		oplot,lambda_ang,0.4*(lambda_ang/2.0d2)^(-1.7),COLOR=150
		oplot,lambda_ang,15.*0.4*(lambda_ang/2.0d2)^(-1.7),COLOR=150
		oplot,lambda_ang,3.*15.*0.4*(lambda_ang/2.0d2)^(-1.5),COLOR=150

		fm = marconi_return_spectrum(nu,1.0d12)
		fm = 10^(fm-12.)
			fm_lam = fm/lambda
			fmnorm = INTERPOL(fm_lam,lambda_ang,[3.1d4])
			fm_lam = fm_lam * (1.0d-4 / fmnorm(0))
			oplot,lambda_ang,fm_lam,color=250


if (2 EQ 0) then begin
	plot,nu,f,/xlog,/ylog
		oplot,nu,fm,color=250
		nu0 = c/(3.6d-6)
		oplot,[nu0,nu0],[1.0d-10,1.0d10]
		nu0 = c/(4.5d-6)
		oplot,[nu0,nu0],[1.0d-10,1.0d10]
		nu0 = c/(5.8d-6)
		oplot,[nu0,nu0],[1.0d-10,1.0d10]
		nu0 = c/(8.0d-6)
		oplot,[nu0,nu0],[1.0d-10,1.0d10]
		nu0 = c/(24.0d-6)
		oplot,[nu0,nu0],[1.0d-10,1.0d10]

		nu0 = c/(3.6d-6) * 2.5
		oplot,[nu0,nu0],[1.0d-10,1.0d10],linestyle=2
		nu0 = c/(4.5d-6) * 2.5
		oplot,[nu0,nu0],[1.0d-10,1.0d10],linestyle=2
		nu0 = c/(5.8d-6) * 2.5
		oplot,[nu0,nu0],[1.0d-10,1.0d10],linestyle=2
		nu0 = c/(8.0d-6) * 2.5
		oplot,[nu0,nu0],[1.0d-10,1.0d10],linestyle=2
		nu0 = c/(24.0d-6) * 2.5
		oplot,[nu0,nu0],[1.0d-10,1.0d10],linestyle=2
	oplot,nu,0.120*(nu/(c/(3.6d-6)))^(2.0),linestyle=1
	oplot,nu,0.120*(nu/(c/(3.6d-6)))^(0.1),linestyle=1
	oplot,nu,0.320*(nu/(c/(3.6d-6)))^(0.7),linestyle=5,color=150
	oplot,nu,0.150*(nu/(c/(3.6d-6)))^(-0.7),linestyle=5,color=110
endif
;;;
;;; end story -- their photo-z's are *way* too inaccurate ( & their accuracy changes as a function 
;;;   of redshift, which could create severe bias), they ignore the hot dust when doing the 
;;;   photo-z's then fit it for the mid-IR excess -- their AGN template (at least for the photo-zs) 
;;;   looks nothing like a real AGN -- and so at the end of the day, their 
;;;   calculated 3.6 micron and 4.5 micron LFs are seriously discrepant from the longer wavelength LFs
;;;

		!P.COLOR=0
	
;	L_IR   = (DOUBLE(L_15-L_solar)) - 5.*alog10(7./7.5)
;	PHI_IR = P_15 + alog10(2.5) + 3.*alog10(7./7.5)
;	DPHI_IR= D_15
end
