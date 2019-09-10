from data import *
import numpy as np 
from lf_shape import *
import scipy.interpolate as inter
from convolve import *
from scipy.optimize import minimize
from scipy.optimize import leastsq
import lmfit
import corner
# fit the luminosity function based on datasets at a given redshift
from lf_fitter_data_limited import *
from ctypes import *
import ctypes
import sys
redshift=float(sys.argv[1])
dtg = return_dtg(redshift)
FIX = int(sys.argv[2])
FIX = bool(FIX)

parameters_init = np.array([0.41698725, 2.17443860, -4.82506430, 13.03575300, 0.4])
parameters_info = np.array(["gamma1_0", "gamma2_0", "logphis"  , "logLs_0", "dtg"])
parameters_bound=(np.array([-np.inf,-np.inf,-np.inf,-np.inf,0]),np.array([0,0,np.inf,np.inf,np.inf]))

#load the shared object file
c_extenstion = CDLL(homepath+'codes/c_lib/convolve.so')
convolve_c = c_extenstion.convolve
convolve_c.restype = ctypes.POINTER(ctypes.c_double * N_bol_grid)

def get_fit_data(alldata,parameters,zmin,zmax,dset_name,dset_id):

	alldata_tem={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([])}
	if load_LF_data[dset_name](redshift)!=False:
		L_data, PHI_data, DPHI_data = load_LF_data[dset_name](redshift)
	else:
		return False

	if dset_id!=-5: return False

	'''
	if dset_id==-5: 
		L_1450 = bolometric_correction(L_bol_grid,dset_id) + L_solar
                M_1450 = -2.5*( L_1450 - np.log10(Fab*con.c.value/1450e-10) )
                L_tmp  = np.sort(M_1450)
	else: L_tmp=bolometric_correction(L_bol_grid,dset_id)

	if dset_id!=-5:
		if return_LF[dset_name]!=None:
			phi_fit_tmp = return_LF[dset_name](L_tmp, redshift)
			phi_fit_pts = np.interp(L_data ,L_tmp, phi_fit_tmp)
			PHI_data = PHI_data + (np.mean((phi_fit_pts))-np.mean((PHI_data)))	
	'''
	if (len(L_data) > 0):
			L_model = bolometric_correction(L_bol_grid,dset_id)
			nu_c = c_double(dset_id)
			redshift_c = c_double(redshift)
			dtg_c = c_double(dtg)
			input_c= np.power(10.,LF(L_bol_grid,parameters)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
			res = convolve_c(input_c,nu_c,redshift_c,dtg_c)
			res = [i for i in res.contents]
			PHI_model = np.array(res,dtype=np.float64)
			#L_model, PHI_model = convolve(np.power(10.,LF_at_z(L_bol_grid,parameters,redshift,"Fiducial")), dset_id) 
			
			if dset_id==-5:
				L_1450 = (-0.4*L_data) + np.log10(Fab*(con.c.value/1450e-10)) - L_solar
                                phi_i = np.interp(L_1450, L_model, np.log10(PHI_model))
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

	#if dset_id==-4: print "NAME:",dset_name,";  CHISQ:", np.sum(((alldata_tem["P_PRED"]-alldata_tem["P_OBS"])/alldata_tem["D_OBS"])**2)," / ",len(alldata_tem["L_OBS"])

def chisq(pars):
	parvals = pars.valuesdict()
    	gamma1 = parvals['gamma1']
    	gamma2 = parvals['gamma2']
    	logphis= parvals['logphis']
	Lbreak = parvals['Lbreak']
	parameters=np.array([gamma1,gamma2,logphis,Lbreak])

        alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z_TOT":np.array([]),"B":np.array([]),"ID":np.array([])}
        for key in dset_ids.keys():
                get_fit_data(alldata,parameters,zmins[key],zmaxs[key],key,dset_ids[key])

        bad = np.invert(np.isfinite(alldata["P_PRED"]))
        if (np.count_nonzero(bad) > 0): alldata["P_PRED"][bad] = -40.0

        chitot = np.sum(((alldata["P_PRED"]-alldata["P_OBS"])/alldata["D_OBS"])**2)
        print chitot, len(alldata["L_OBS"])
        return chitot#, len(alldata["L_OBS"])

def residual(pars):
	parvals = pars.valuesdict()
        gamma1 = parvals['gamma1']
        gamma2 = parvals['gamma2']
        logphis= parvals['logphis']
        Lbreak = parvals['Lbreak']
        parameters=np.array([gamma1,gamma2,logphis,Lbreak])

	alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z_TOT":np.array([]),"B":np.array([]),"ID":np.array([])}
	for key in dset_ids.keys():
		get_fit_data(alldata,parameters,zmins[key],zmaxs[key],key,dset_ids[key])

	bad = np.invert(np.isfinite(alldata["P_PRED"]))
	if (np.count_nonzero(bad) > 0): alldata["P_PRED"][bad] = -40.0

	chitot = np.sum(((alldata["P_PRED"]-alldata["P_OBS"])/alldata["D_OBS"])**2)
	#print chitot, len(alldata["L_OBS"])
	return (alldata["P_PRED"]-alldata["P_OBS"])/alldata["D_OBS"]

params = lmfit.Parameters()
# add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
if FIX==False:
	params.add_many(('gamma1' , parameters_init[0], True, None, None, None, None),
                	('gamma2' , parameters_init[1], True, None, None, None, None),
                        ('logphis', parameters_init[2], True, None, None, None, None),
                        ('Lbreak' , parameters_init[3], True, None, None, None, None))
else:
	if (redshift>=0.4) and (redshift<=2.8):
                params.add_many(('gamma1' , parameters_init[0], True, None, None, None, None),
                                ('gamma2' , parameters_init[1], True, None, None, None, None),
                                ('logphis', parameters_init[2], True, None, None, None, None),
                                ('Lbreak' , parameters_init[3], True, None, None, None, None))
	elif (redshift<5.8):
		logphis_fixed=-3.94340268-0.21766248*(1+redshift)
		params.add_many(('gamma1' , parameters_init[0], True, None, None, None, None),
                	        ('gamma2' , parameters_init[1], True, None, None, None, None),
                	        ('logphis', logphis_fixed     ,False, None, None, None, None),
                	        ('Lbreak' , parameters_init[3], True, None, None, None, None))
	else: # at z>5.8, use single power law to do te fit
		logphis_fixed=-3.94340268-0.21766248*(1+redshift)
		params.add_many(('gamma1' , parameters_init[0], True, None, None, None, None),
                	        ('gamma2' , parameters_init[1], True, None, None, "gamma1", None),
                	        ('logphis', logphis_fixed     ,False, None, None, None, None),
                	        ('Lbreak' , parameters_init[3], True, None, None, None, None))

fitter = lmfit.Minimizer(residual, params, scale_covar=True,nan_policy='raise',calc_covar=True)
result=fitter.minimize(method='leastsq')
#print "bestfit:"
#result.params.pretty_print()

print redshift, result.params['gamma1'].value, result.params['gamma1'].stderr, result.params['gamma2'].value, result.params['gamma2'].stderr, result.params['logphis'].value, result.params['logphis'].stderr, result.params['Lbreak'].value, result.params['Lbreak'].stderr, result.nfree, result.redchi
