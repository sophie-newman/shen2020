;;; 2SLAQ (Richards et al. -- supplants Boyle et al. & Croom et al.)
;;;     rest-frame B-band
;;;
pro load_2SLAQ_LF_data, L_BB, PHI_BB, DPHI_BB, z
	forward_function M_sun_BB
	;; determine which redshift interval its in (bins follow Croom 2004)
	if (z LT 0.4) then WHICH_BLOCK = 0
	if ((z GE 0.40) AND (z LT 0.68)) then WHICH_BLOCK = 1
	if ((z GE 0.68) AND (z LT 0.97)) then WHICH_BLOCK = 2
	if ((z GE 0.97) AND (z LT 1.25)) then WHICH_BLOCK = 3
	if ((z GE 1.25) AND (z LT 1.53)) then WHICH_BLOCK = 4
	if ((z GE 1.53) AND (z LT 1.81)) then WHICH_BLOCK = 5
	if ((z GE 1.81) AND (z LE 2.10)) then WHICH_BLOCK = 6
	if (z GT 2.1) then WHICH_BLOCK = 7

	if ((WHICH_BLOCK EQ 0) OR (WHICH_BLOCK EQ 7)) then begin
		L_BB = 0.
		PHI_BB = 0.
		DPHI_BB = 0.
	endif else begin

	OPENR, 1, return_idl_routines_homedir(0)+'/luminosity_functions/load_2slaq_lf_data.dat'

	N_LUM_BINS = 34
	N_Z_BINS   = 6
	z_list = fltarr(N_LUM_BINS,N_Z_BINS)
	M_B    = z_list
	phi    = z_list
	dphi   = z_list
	nqso   = z_list
	bin_flags = intarr(N_LUM_BINS,N_Z_BINS)

	for i = 0, N_Z_BINS-1 do begin
	for j = 0, N_LUM_BINS-1 do begin
		block = fltarr(6)
		READF, 1, block
		z_list[j,i] = block[0]
		M_B[j,i] = block[1]
		phi[j,i] = block[2]
		dphi[j,i] = block[3]
		nqso[j,i] = block[4]
		bin_flags[j,i] = block[5]
	endfor
	endfor
	CLOSE, 1
	M_B = M_B[*,WHICH_BLOCK-1] + 0.05	;; bJ = g + 0.05 from Fukugita et al. (R06 DR3 ppr)
	phi = phi[*,WHICH_BLOCK-1]
	dphi = dphi[*,WHICH_BLOCK-1]
	bin_flags = bin_flags[*,WHICH_BLOCK-1]
	ok = where(bin_flags GE 0)
	
	L_BB = (0.4*(M_sun_BB(0)-M_B[ok])) 
	PHI_BB = alog10(phi[ok]) + alog10(2.5)
	DPHI_BB = alog10(1.+dphi[ok]/phi[ok])
	endelse
end

