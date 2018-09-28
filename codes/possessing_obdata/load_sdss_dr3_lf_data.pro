;;;
;;; Richards et al. SDSS DR3
;;;   Measuring M_i(z=2) (outputs as M_B)
;;;
pro load_SDSS_DR3_LF_data, L_BB, PHI_BB, DPHI_BB, z
	forward_function M_sun_BB
	;; determine which redshift interval its in (bins follow Croom 2004)
	if ((z LT 0.4) OR (Z GT 5.0)) then begin
		L_BB = 0.
		PHI_BB = 0.
		DPHI_BB = 0.
	endif else begin	
	if ((z GE 0.40) AND (z LT 0.68)) then filena = 'z.0.49'
	if ((z GE 0.68) AND (z LT 1.06)) then filena = 'z.0.87'
	if ((z GE 1.06) AND (z LT 1.44)) then filena = 'z.1.25'
	if ((z GE 1.44) AND (z LT 1.82)) then filena = 'z.1.63'
	if ((z GE 1.82) AND (z LT 2.21)) then filena = 'z.2.01'
	if ((z GE 2.21) AND (z LE 2.60)) then filena = 'z.2.40'
	if ((z GE 2.60) AND (z LE 3.03)) then filena = 'z.2.80'
	if ((z GE 3.03) AND (z LE 3.50)) then filena = 'z.3.25'
	if ((z GE 3.50) AND (z LE 4.00)) then filena = 'z.3.75'
	if ((z GE 4.00) AND (z LE 4.50)) then filena = 'z.4.25'
	if ((z GE 4.50) AND (z LE 5.00)) then filena = 'z.4.75'

	filena = return_idl_routines_homedir(0)+'/luminosity_functions/load_sdss_dr3_lf_data.dat/'+filena
	OPENR, lunner, filena, /get_lun

	N_M_BINS = 0
	READF,lunner,N_M_BINS
	M_i  = fltarr(N_M_BINS)
	phi  = fltarr(N_M_BINS)
	dphi = fltarr(N_M_BINS)
	fill = fltarr(N_M_BINS)
	zbar = fltarr(N_M_BINS)
	Nqso = fltarr(N_M_BINS)
	Ncor = fltarr(N_M_BINS)
	for i = 0, N_M_BINS-1 do begin
		dumm = fltarr(7)
		READF,lunner,dumm
		M_i[i]  = dumm[0]
		phi[i]  = dumm[1]
		dphi[i] = ((dumm[2]*1.0d-9)/(10^(dumm[1])))/alog(10.)
	endfor
	M_B = M_i + 0.66	;; Correction from the paper from Mi(z=2) to M_Bj
	M_B = M_i + 0.71	;; Correction from the paper from Mi(z=2) to M_Bj
		
	L_BB = (0.4*(M_sun_BB(0)-M_B))
	PHI_BB = phi + alog10(2.5)
	DPHI_BB = dphi

	CLOSE,lunner
	FREE_LUN,lunner
	endelse

end
