from data import *
import numpy as np 
from lf_shape import *
import scipy.interpolate as inter
from convolve import *
from scipy.optimize import curve_fit

# fit the luminosity function based on datasets at a given redshift
#L_BB_tmp = bolometric_correction(L_bol_grid,-1)

dset_ids={
	"DR3":0,
	"2SLAQ":0,
	"HUNT":0
}

zmins={
	"DR3": np.array([0.40, 0.68, 1.06, 1.44, 1.82, 2.21, 2.60, 3.03, 3.50, 4.00, 4.50]),
	"2SLAQ":np.array([0.40, 0.68, 0.97, 1.25, 1.53, 1.81]),
	"HUNT":np.array([2.0]),
	"SIANA":np.array([2.9]),
	"CRISTIANI":np.array([4.0]),
	"KDC":np.array([4.0]),
	"FAN":np.array([3.6, 3.9, 4.4, 5.8]),
	"SSG":np.array([2.75]),
	"COMBO17":np.array([1.2, 1.8, 2.4, 3.0, 3.6, 4.2]),
	#soft xray
	"HMS":np.array([0.0, 0.2, 0.4, 0.8, 1.6, 3.2]),
	"MIYAJI":np.array([0.0, 0.2, 0.4, 0.8, 1.6, 2.3]),
	#hard xray
	"UEDA":np.array([0.0, 0.2, 0.4, 0.8, 1.6]),
	"LAFRANCA":np.array([0.0, 0.5, 1.0, 1.5, 2.5]),
	"SILVERMAN_SX":np.array([0.5, 1.0, 2.0, 4.0]),
	"SILVERMAN_HX":np.array([0.2, 1.5, 3.0]),
	"BARGER":np.array([0.1, 0.4, 0.8, 1.5, 3.0]),
	"NANDRA":np.array([2.75]),
	"SAZREV":np.array([0.0]),
	#INFRARED
	"BROWN":np.array([1.5]),
	"MATUTE":np.array([0.0, 0.5]),
	#narrow lines
	"HAO":np.array([0.00])
}
zmaxs={
	"DR3": np.array([0.68, 1.06, 1.44, 1.82, 2.21, 2.60, 3.03, 3.50, 4.00, 4.50, 5.00]),
	"2SLAQ":np.array([0.68, 0.97, 1.25, 1.53, 1.81, 2.10]),
	"HUNT":np.array([4.0]),
	"SIANA":np.array([3.4]),
	"CRISTIANI":np.array([5.2]),
	"KDC":np.array([4.5]),
	"FAN":np.array([3.9, 4.4, 5.0, 6.2]),
	"SSG":np.array([4.75]),
	"COMBO17":np.array([1.8, 2.4, 3.0, 3.6, 4.2, 4.8]),
	"HMS":np.array([0.2, 0.4, 0.8, 1.6, 3.2, 4.8]),
	"MIYAJI":np.array([0.2, 0.4, 0.8, 1.6, 2.3, 4.6]),
	"UEDA":np.array([0.2, 0.4, 0.8, 1.6, 3.0]),
	"LAFRANCA":np.array([0.5, 1.0, 1.5, 2.5, 3.5]),
	"SILVERMAN_SX":np.array([1.0, 1.5, 3.0, 5.5]),
	"SILVERMAN_HX":np.array([0.5, 2.0, 4.0]),
	"BARGER":np.array([0.4, 0.8, 1.2, 3.0, 5.0]),
	"NANDRA":np.array([3.25]),
	"SAZREV":np.array([0.1]),
	"BROWN":np.array([2.5]),
	"MATUTE":np.array([0.2, 2.0]),
	"HAO":np.array([0.10])
}
load_LF_data={
	"DR3":load_sdss_dr3_lf_data,
	"2SLAQ":load_2slaq_lf_data,
	"HUNT":load_hunt_lf_data,
	"SIANA":load_siana_lf_data,
	"CRISTIANI":load_cristiani_lf_data,
	"KDC":load_kdc_lf_data,
	"FAN":load_fan_lf_data,
	"SSG":load_ssg_lf_data,
	"COMBO17":load_combo17_lf_data,
	"HMS":load_hms_lf_data,
	"MIYAJI":load_miyaji_lf_data,
	"UEDA":load_ueda_lf_data,
	"LAFRANCA":load_lafranca_lf_data,
	"SILVERMAN_SX":load_silverman_sx_lf_data,
	"SILVERMAN_HX":load_silverman_hx_lf_data,
	"BARGER":load_barger_lf_data,
	"NANDRA":load_nandra_lf_data,
	"SAZREV":load_sazrev_lf_data,
	"BROWN":load_brown_lf_data,
	"MATUTE":load_matute_lf_data,
	"HAO":load_hao_lf_data
}

parameters_init = np.array([0.10 , 	0.608, 	-6.58, 	8.427,  1.02, 
					    -10.23, -3.00,	0.000, 	1.403,  -1.07,  
					    0.0, 	0.00,   0.0, 	1.5, 	-3.19])		
parameters_info = np.array(["nocon", "nocon", "nocon", "nocon", "nocon", 
		           "nocon", "nocon", "fix"	, "nocon", "nocon", 
				   "fix",  	"fix",	 "fix",   "nocon", "nocon"])

def get_fit_data(alldata,parameters,zmin,zmax,dset_name,dset_id):
	z_lis = 0.5*(zmin + zmax)
	for iz in range(len(z_lis)):
		redshift = z_lis[iz]
		L_BB, PHI_BB, DPHI_BB = load_LF_data[dset_name](redshift)
		#phi_fit_tmp = np.log10( return_LF[dset_name](L_BB_tmp, redshift))
		#phi_fit_pts = inter.interp1d(L_BB_tmp, phi_fit_tmp)(L_BB)

		#if (KEYS["RENORM_IN_FITTING"]==True):
		#	PHI_BB = PHI_BB + (MEAN((phi_fit_pts))-MEAN((PHI_BB)))	
		if (len(L_BB) > 1):
			L_B, PHI_B = convolve(LF_at_z(L_bol_grid,parameters,redshift,"Fiducial"), -1) 

			phi_i = np.power(10.,interp1d(LB, np.log10(PHI_B))(L_BB))

			alldata["P_PRED"] = np.append(alldata["P_PRED"] , phi_i)
			alldata["L_OBS"]  = np.append(alldata["L_OBS"]  , L_BB)
			alldata["P_OBS"]  = np.append(alldata["P_OBS"]  , PHI_BB)
			alldata["D_OBS"]  = np.append(alldata["D_OBS"]  , DPHI_BB)
			alldata["Z_TOT"]  = np.append(alldata["Z_TOT"]  , np.ones(len(phi_i)) * redshift)
			alldata["B"]      = np.append(alldata["B"]      , 0*phi_i)
			alldata["ID"]     = np.append(alldata["ID"]     , np.ones(len(phi_i)) * dset_id)

			if KEYS["VERBOSE"]==True:
				if (iz == 0): chi2tott, doftott  = 0., 0.			
				chi2tott += np.sum(((phi_i-PHI_BB)/DPHI_BB)**2)
				doftott  += len(PHI_BB)
				if (iz == (len(z_lis)-1)):  print dset_name+' ->',chi2tott,doftott

def get_lf_pred(parameters): 
	alldata_tem={"P_PRED":0,"L_OBS":0,"P_OBS":0,"D_OBS":0,"Z_TOT":0,"B":0,"ID":0}
	for key in dset_ids.keys():
		get_fit_data(alldata_tem,parameters,zmins[key],zmaxs[key],key,dset_ids[key])
	bad = np.invert(np.isfinite(alldata_tem["P_PRED"]))
	if (np.count_nonzero(bad) > 0): alldata_tem["P_PRED"][bad] = -40.0

	print np.sum(((alldata_tem["P_PRED"]-alldata_tem["P_OBS"])/alldata_tem["D_OBS"])**2)
	print len(alldata_tem["L_OBS"])

	return alldata_tem["P_PRED"]

parameters_init=np.array([   ])

alldata={"P_PRED":0,"L_OBS":0,"P_OBS":0,"D_OBS":0,"Z_TOT":0,"B":0,"ID":0}
for key in dset_ids.keys():
	get_fit_data(alldata,parameters_init,zmins[key],zmaxs[key],key,dset_ids[key])

paremeters_fit, pcov = curve_fit(get_lf_pred,alldata["L_OBS"],alldata["P_OBS"],p0=parameters_init,sigma=alldata["D_OBS"],maxfev=100000)

print parameters

