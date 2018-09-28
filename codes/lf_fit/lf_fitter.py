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
	"2SLAQ":0,
	"HUNT":0
}
zmaxs={
	"DR3": np.array([0.68, 1.06, 1.44, 1.82, 2.21, 2.60, 3.03, 3.50, 4.00, 4.50, 5.00]),
	"2SLAQ":0,
	"HUNT":0
}
load_LF_data={
	"DR3":0,
	"2SLAQ":0,
	"HUNT":0
}
return_LF={
	"DR3":0,
	"2SLAQ":0
}

def get_fit_data(alldata,parameters,zmin,zmax,dset_name,dset_id):
	z_lis = 0.5*(z_min + z_max)
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

paremeters_init=np.array([   ])

alldata={"P_PRED":0,"L_OBS":0,"P_OBS":0,"D_OBS":0,"Z_TOT":0,"B":0,"ID":0}
for key in dset_ids.keys():
	get_fit_data(alldata,parameters_init,zmins[key],zmaxs[key],key,dset_ids[key])

paremeters_fit,_ = curve_fit(get_lf_pred,alldata["L_OBS"],alldata["P_OBS"],p0=parameters_init,sigma=alldata["D_OBS"])

print parameters

