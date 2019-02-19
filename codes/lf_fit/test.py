import os
os.environ["OMP_NUM_THREADS"] = "1"

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
import sys
redshift=float(sys.argv[1])

parameters_init = np.array([0.41698725, 2.17443860, -4.82506430, 13.03575300])
parameters_info = np.array(["gamma1_0", "gamma2_0", "logphis"  , "logLs_0"])
parameters_bound=(np.array([-np.inf,-np.inf,-np.inf,-np.inf]),np.array([0,0,np.inf,np.inf]))

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
			input_c= np.power(10.,LF(L_bol_grid,parameters)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
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

	#if dset_id==-1: print "NAME:",dset_name,";  CHISQ:", np.sum(((alldata_tem["P_PRED"]-alldata_tem["P_OBS"])/alldata_tem["D_OBS"])**2)," / ",len(alldata_tem["L_OBS"])

def residual(pars):
	parameters = pars

	alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z_TOT":np.array([]),"B":np.array([]),"ID":np.array([])}
	for key in dset_ids.keys():
		get_fit_data(alldata,parameters,zmins[key],zmaxs[key],key,dset_ids[key])

	bad = np.invert(np.isfinite(alldata["P_PRED"]))
	if (np.count_nonzero(bad) > 0): alldata["P_PRED"][bad] = -40.0

	chitot = np.sum(((alldata["P_PRED"]-alldata["P_OBS"])/alldata["D_OBS"])**2)
	#print chitot, len(alldata["L_OBS"])
	return (alldata["P_PRED"]-alldata["P_OBS"])/alldata["D_OBS"], alldata["D_OBS"]

def lnlike(pars):
	res, sigma = residual(pars)
	return -0.5*np.sum(res**2 + np.log(2*np.pi*sigma**2))

def lnprior(pars):
	gamma1, gamma2, logphis, Lbreak = pars
        return 0.0

def lnprob(pars):
	lp = lnprior(pars)
	if not np.isfinite(lp):  return -np.inf
	return lp + lnlike(pars)

ndim, nwalkers = 4, 8
pos = np.array([0.5 * np.random.randn(ndim) for i in range(nwalkers)])
pos = pos + parameters_init

with MPIPool() as pool:
	if not pool.is_master():
		pool.wait()
		sys.exit(0)

	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool)
	nsteps = 40
	#width = 80
	print "begin sampling"
	sampler.run_mcmc(pos, nsteps, progress=True)
	#for i, result in enumerate(sampler.sample(pos, iterations=nsteps)):
	#	if (i+1) % 20 == 0:
        #                print("{0:5.1%}".format(float(i) / nsteps))

print "end sampling"

burn = 20
samples = sampler.chain[:, burn:, :].reshape((-1, ndim))

np.save("dummy",samples)

print "chain saved"

fig = corner.corner(samples, labels=["g1", "g2", "Phi" , "L"])
fig.savefig("triangle.png")
#np.savetxt("chain.dat",result.flatchain)
