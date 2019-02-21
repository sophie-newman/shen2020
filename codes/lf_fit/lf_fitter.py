from data import *
import numpy as np 
from lf_shape import *
import scipy.interpolate as inter
from convolve import *

import emcee
from emcee.utils import MPIPool
import corner
import matplotlib.pyplot as plt
# fit the luminosity function based on datasets at a given redshift
from lf_fitter_data import *
from ctypes import *
import ctypes

import time

#parameters_init = np.array([0.41698725, 2.17443860, -4.82506430, 13.03575300, 0.63150872, -11.76356000, -14.24983300, -0.62298947, 1.45993930, -0.79280099])
#parameters_info = np.array(["gamma1_0", "gamma2_0", "logphis"  , "logLs_0"  , "k1"      , "k2"        , "k3"        , "k_gamma1" , "k_gamma2_1", "k_gamma2_2"])
'''
parameters_init = np.array([4.5638212952980256e-01,-5.3535826442735933e-02,8.291438766390475321e-03,
 2.23952286876408957e+00,9.6563649154324882e-01,-1.46845995121186257e+00,7.653560660287425099e-01,
 -3.855561029913782800e+00, -3.548048158301563837e-01,
 9.92811918927010062e+00,-1.92756399797551145e-03,-1.2185864295518154e+00,2.072801872747117857e-01])
'''
parameters_init = np.array([0.05,0.20,-0.004,  2.4,1.0,-1.7,0.8,  -3.85,-0.6,   12.0,0.5,-0.6,0.16])

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
			alldata["Z"]  = np.append(alldata["Z"]  , np.ones(len(L_data)) * redshift)
			alldata["WEIGHT"]  = np.append(alldata["WEIGHT"]  , np.ones(len(L_data)) * np.sqrt(1+redshift))
			alldata["ID"] = np.append(alldata["ID"] , np.ones(len(L_data)) * dset_id)
			
			alldata_tem["P_PRED"] = np.append(alldata_tem["P_PRED"] , phi_i)
			alldata_tem["L_OBS"]  = np.append(alldata_tem["L_OBS"]  , L_data)
			alldata_tem["P_OBS"]  = np.append(alldata_tem["P_OBS"]  , PHI_data)
			alldata_tem["D_OBS"]  = np.append(alldata_tem["D_OBS"]  , DPHI_data + 0.01)

			#print "NAME:",dset_name,"; redshift",redshift,";  chisq:", np.sum(((phi_i-PHI_BB)/DPHI_BB)**2)," / ",len(L_BB)

	#print "NAME:",dset_name,";  CHISQ:", np.sum(((alldata_tem["P_PRED"]-alldata_tem["P_OBS"])/alldata_tem["D_OBS"])**2)," / ",len(alldata_tem["L_OBS"])

def residual(pars):
	parameters = pars

	alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z":np.array([]),"WEIGHT":np.array([]),"ID":np.array([])}
	for key in dset_ids.keys():
		get_fit_data(alldata,parameters,zmins[key],zmaxs[key],key,dset_ids[key])

	bad = np.invert(np.isfinite(alldata["P_PRED"]))
	if (np.count_nonzero(bad) > 0): alldata["P_PRED"][bad] = -40.0

	chitot = np.sum( ((alldata["P_PRED"]-alldata["P_OBS"])/alldata["D_OBS"])**2)
	#print chitot, len(alldata["L_OBS"])
	return (alldata["P_PRED"]-alldata["P_OBS"])/alldata["D_OBS"], alldata["D_OBS"]

def lnlike(pars):
        res, sigma = residual(pars)
        return -0.5*np.sum(res**2 + np.log(2*np.pi*sigma**2))

def lnprior(pars):
        return 0.0

def lnprob(pars):
        lp = lnprior(pars)
	llike = lnlike(pars)
        if not np.isfinite(lp):  return -np.inf
	if not np.isfinite(llike):  return -np.inf
        return lp + llike

start_time = time.time()

ndim, nwalkers = 13, 100
pos = np.array([np.random.randn(ndim) for i in range(nwalkers)])
for i in range(pos.shape[0]):
	pos[i,:] = pos[i,:] * 0.2 * parameters_init + parameters_init + pos[i,:] * 0.1

f = open("output/chain.dat", "w")
f.close()

with MPIPool() as pool:
	if not pool.is_master():
        	pool.wait()
        	sys.exit(0)

	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool)
	nsteps = 4000
	#sampler.run_mcmc(pos, nsteps)
	print "begin sampling"
	for i, result in enumerate(sampler.sample(pos, iterations=nsteps)):
        	#if pool.is_master(): 
		if (i+1) % 10 == 0:
        		print "{0:5.1%}".format(float(i) / nsteps), "  t:", time.time()-start_time 
		
		position = result[0]
		f = open("output/chain.dat", "a")
		for k in range(position.shape[0]):
			f.write("{0:3d} {1:7d} {2:s}\n".format(k, i, np.array2string(position[k]).strip('[').strip(']').replace('\n',' ') ))
		f.close()
		
burn = 200
samples = sampler.chain[:, burn:, :].reshape((-1, ndim))

np.save("output/chain_main.npy",samples)

