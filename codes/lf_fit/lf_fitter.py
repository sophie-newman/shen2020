from data import *
import numpy as np 
from lf_shape import *
import scipy.interpolate as inter
from convolve import *
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import least_squares
# fit the luminosity function based on datasets at a given redshift
from lf_fitter_data import *
from ctypes import *
import ctypes

bestfit2        = np.array([0.39856372, 2.19426155, -4.73290265, 12.97241616, 0.43717880, -11.63470663, -11.72357430, -0.72820353, 1.36234503, -0.79697647])
bestfit         = np.array([0.38181256, 2.16955741, -4.70557406, 12.94851374, 0.43771614, -11.42561263, -11.34952214, -0.75528960, 1.32130027, -0.77768681])
parameters_init = np.array([0.41698725, 2.17443860, -4.82506430, 13.03575300, 0.63150872, -11.76356000, -14.24983300, -0.62298947, 1.45993930, -0.79280099])
parameters_info = np.array(["gamma1_0", "gamma2_0", "logphis"  , "logLs_0"  , "k1"      , "k2"        , "k3"        , "k_gamma1" , "k_gamma2_1", "k_gamma2_2"])
parameters_bound= (np.array([0,0,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf]),np.array([np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]))

#load the shared object file
c_extenstion = CDLL(homepath+'codes/c_lib/convolve.so')
convolve_c = c_extenstion.convolve
convolve_c.restype = ctypes.POINTER(ctypes.c_double * N_bol_grid)

def get_fit_data(alldata,parameters,zmin,zmax,dset_name,dset_id):
	z_lis = 0.5*(zmin + zmax)
	
	alldata_tem={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([])}

	for iz in range(len(z_lis)):
		redshift = z_lis[iz]
		L_data, PHI_data, DPHI_data = load_LF_data[dset_name](redshift)
		
		#L_tmp=bolometric_correction(L_bol_grid,dset_id)
		#if return_LF[dset_name]!=None:
		#	phi_fit_tmp = return_LF[dset_name](L_tmp, redshift)
		#	phi_fit_pts = np.interp(L_data ,L_tmp, phi_fit_tmp)
		#	PHI_data = PHI_data + (np.mean((phi_fit_pts))-np.mean((PHI_data)))	
		if (len(L_data) > 0):
			if dset_id==-1.1: 
				L_model = bolometric_correction(L_bol_grid,-1)
				nu_c = c_double(-1)
			else:
				L_model = bolometric_correction(L_bol_grid,dset_id)
				nu_c = c_double(dset_id)
			input_c= np.power(10.,LF_at_z(L_bol_grid,parameters,redshift,"Fiducial")).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
			res = convolve_c(input_c,nu_c)
			res = [i for i in res.contents]
			PHI_model = np.array(res,dtype=np.float64)
			#L_model, PHI_model = convolve(np.power(10.,LF_at_z(L_bol_grid,parameters,redshift,"Fiducial")), dset_id) 
			
			phi_i = np.interp(L_data, L_model, np.log10(PHI_model))
			alldata["P_PRED"] = np.append(alldata["P_PRED"] , phi_i)
			alldata["L_OBS"]  = np.append(alldata["L_OBS"]  , L_BB)
			alldata["P_OBS"]  = np.append(alldata["P_OBS"]  , PHI_BB)
			alldata["D_OBS"]  = np.append(alldata["D_OBS"]  , DPHI_BB)# + 0.01)
			alldata["Z_TOT"]  = np.append(alldata["Z_TOT"]  , np.ones(len(L_BB)) * redshift)
			alldata["ID"]     = np.append(alldata["ID"]     , np.ones(len(L_BB)) * dset_id)
			
			alldata_tem["P_PRED"] = np.append(alldata_tem["P_PRED"] , phi_i)
			alldata_tem["L_OBS"]  = np.append(alldata_tem["L_OBS"]  , L_BB)
			alldata_tem["P_OBS"]  = np.append(alldata_tem["P_OBS"]  , PHI_BB)
			alldata_tem["D_OBS"]  = np.append(alldata_tem["D_OBS"]  , DPHI_BB)# + 0.01)

			#print "NAME:",dset_name,"; redshift",redshift,";  chisq:", np.sum(((phi_i-PHI_BB)/DPHI_BB)**2)," / ",len(L_BB)
	print "NAME:",dset_name,";  CHISQ:", np.sum(((alldata_tem["P_PRED"]-alldata_tem["P_OBS"])/alldata_tem["D_OBS"])**2)," / ",len(alldata_tem["L_OBS"])

def chisq(parameters):
	alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z_TOT":np.array([]),"B":np.array([]),"ID":np.array([])}
	for key in dset_ids.keys():
		get_fit_data(alldata,parameters,zmins[key],zmaxs[key],key,dset_ids[key])

	bad = np.invert(np.isfinite(alldata["P_PRED"]))
	if (np.count_nonzero(bad) > 0): alldata["P_PRED"][bad] = -40.0

	return np.sum(((alldata["P_PRED"]-alldata["P_OBS"])/alldata["D_OBS"])**2)#,len(alldata["L_OBS"])-10

out = minimize(chisq,x0=parameters_init,method='Nelder-Mead')#bounds=parameters_bound)
print out



