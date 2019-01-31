from data import *
import numpy as np 
from lf_shape import *
import scipy.interpolate as inter
from convolve import *
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import least_squares
import lmfit
# fit the luminosity function based on datasets at a given redshift
from lf_fitter_data import *
from ctypes import *
import ctypes

#parameters_init = np.array([0.41698725, 2.17443860, -4.82506430, 13.03575300, 0.63150872, -11.76356000, -14.24983300, -0.62298947, 1.45993930, -0.79280099])
#parameters_info = np.array(["gamma1_0", "gamma2_0", "logphis"  , "logLs_0"  , "k1"      , "k2"        , "k3"        , "k_gamma1" , "k_gamma2_1", "k_gamma2_2"])
parameters_init = np.array([1.705501200816446294e+00 ,-1.078700539443752326e+00,1.417091812716838317e-01 ,-4.913503567223131800e-03,
2.342208072486966941e+00 , 2.067114074110976141e+00,-1.000403268800853951e+00, 1.399154527410910287e+00,
-3.860125192409941342e+00, -3.531881721705402155e-01,
1.099643031838996166e+01, 7.302490517201445819e+00, -3.252882882822146215e-01, 1.239403614340442328e+00])
parameters_info = np.array(["gamma1", "gamma2", "logphis"  , "logLs"])

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
		
		if dset_id==-1.1: L_tmp=bolometric_correction(L_bol_grid,-1)
        	elif dset_id==-1:
                	L_tmp=bolometric_correction(L_bol_grid,dset_id)
                	L_tmp = (M_sun_Bband_AB -2.5*L_tmp) + 0.706
                	L_tmp = np.sort(L_tmp)
        	else: L_tmp=bolometric_correction(L_bol_grid,dset_id)
        
		if return_LF[dset_name]!=None:
                	phi_fit_tmp = return_LF[dset_name](L_tmp, redshift)
                	phi_fit_pts = np.interp(L_data ,L_tmp, phi_fit_tmp)
                	PHI_data = PHI_data + (np.mean((phi_fit_pts))-np.mean((PHI_data)))

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
			
			if dset_id==-1:
				L_Bband = (M_sun_Bband_AB-(L_data - 0.706))/2.5
				phi_i = np.interp(L_Bband, L_model, np.log10(PHI_model))
				phi_i = phi_i - np.log10(2.5)
			else:
				phi_i = np.interp(L_data, L_model, np.log10(PHI_model))

			alldata["P_PRED"] = np.append(alldata["P_PRED"] , phi_i)
			alldata["L_OBS"]  = np.append(alldata["L_OBS"]  , L_data)
			alldata["P_OBS"]  = np.append(alldata["P_OBS"]  , PHI_data)
			alldata["D_OBS"]  = np.append(alldata["D_OBS"]  , DPHI_data + 0.01)
			alldata["Z_TOT"]  = np.append(alldata["Z_TOT"]  , np.ones(len(L_data)) * redshift)
			alldata["ID"]     = np.append(alldata["ID"]     , np.ones(len(L_data)) * dset_id)
			
			alldata_tem["P_PRED"] = np.append(alldata_tem["P_PRED"] , phi_i)
			alldata_tem["L_OBS"]  = np.append(alldata_tem["L_OBS"]  , L_data)
			alldata_tem["P_OBS"]  = np.append(alldata_tem["P_OBS"]  , PHI_data)
			alldata_tem["D_OBS"]  = np.append(alldata_tem["D_OBS"]  , DPHI_data + 0.01)

			#print "NAME:",dset_name,"; redshift",redshift,";  chisq:", np.sum(((phi_i-PHI_BB)/DPHI_BB)**2)," / ",len(L_BB)

	#print "NAME:",dset_name,";  CHISQ:", np.sum(((alldata_tem["P_PRED"]-alldata_tem["P_OBS"])/alldata_tem["D_OBS"])**2)," / ",len(alldata_tem["L_OBS"])

def residual(pars):
	parvals = pars.valuesdict()
	parameters=np.array([parvals['p0'],parvals['p1'],parvals['p2'],parvals['p3'],parvals['p4'],parvals['p5'],parvals['p6'],parvals['p7'],parvals['p8'],parvals['p9'],parvals['p10'],parvals['p11'],parvals['p12'],parvals['p13']])

	alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z_TOT":np.array([]),"B":np.array([]),"ID":np.array([])}
	for key in dset_ids.keys():
		get_fit_data(alldata,parameters,zmins[key],zmaxs[key],key,dset_ids[key])

	bad = np.invert(np.isfinite(alldata["P_PRED"]))
	if (np.count_nonzero(bad) > 0): alldata["P_PRED"][bad] = -40.0

	chitot = np.sum(((alldata["P_PRED"]-alldata["P_OBS"])/alldata["D_OBS"])**2)
	print chitot, len(alldata["L_OBS"])
	return (alldata["P_PRED"]-alldata["P_OBS"])/alldata["D_OBS"]


params = lmfit.Parameters()
# add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
params.add_many(('p0' , parameters_init[0], True, None, None, None, None),
                ('p1' , parameters_init[1], True, None, None, None, None),
                ('p2' , parameters_init[2], True, None, None, None, None),
                ('p3' , parameters_init[3], True, None, None, None, None),
		('p4' , parameters_init[4], True, None, None, None, None),
		('p5' , parameters_init[5], True, None, None, None, None),
		('p6' , parameters_init[6], True, None, None, None, None),
		('p7' , parameters_init[7], True, None, None, None, None),
		('p8' , parameters_init[8], True, None, None, None, None),
		('p9' , parameters_init[9], True, None, None, None, None),
		('p10' , parameters_init[10], True, None, None, None, None),
		('p11' , parameters_init[11], True, None, None, None, None),
		('p12' , parameters_init[12], True, None, None, None, None),
		('p13' , parameters_init[13], True, None, None, None, None))

fitter = lmfit.Minimizer(residual, params, scale_covar=True,nan_policy='raise',calc_covar=True)
#result=fitter.minimize(method='emcee',burn=300, steps=1000,nwalkers=100,workers=8)
result=fitter.minimize(method='leastsq')
print "bestfit:"
result.params.pretty_print()
print "all messages:"
cov=result.covar
print cov

