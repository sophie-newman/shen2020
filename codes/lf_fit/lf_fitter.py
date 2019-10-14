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

import resource
import time
import os

zref = 2.

def weightfunc(z):
	if z<3.: return 1.
	elif z<4: return np.power( (1.+z)/(1.+3.), 1.)
	else: return np.power( (1.+4.)/(1.+3.), 1.) * np.power( (1.+z)/(1.+4.), 1.5)

#parameters_init = np.array([0.79847831 ,-0.25102297,0.02265792 ,1.846936805,0.316915805,-2.80025713,0.4003701  ,-3.5469101 ,-0.39856649,11.32661766,0.27481868 ,-0.82451577,0.25281821])
parameters_init = np.array([0.79847831 ,-0.25102297,0.02265792 , 1.846936805, -2.80025713 , 0.4003701  ,-3.5469101 ,-0.39856649, 11.32661766 ,-0.82451577,0.25281821])
#parameters_vary = np.array([0.4, 0.2, 0.2, 0.9, 0.3, 1.4, 0.3, 1.7, 0.3, 3., 0.2, 0.4, 0.3])
parameters_vary = np.array([0.4, 0.2, 0.2, 0.9, 1.4, 0.3, 1.7, 0.3, 3., 0.4, 0.3])
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
		
        	if dset_id==-5:
                	L_1450 = bolometric_correction(L_bol_grid,dset_id) + L_solar
                	M_1450 = -2.5*( L_1450 - np.log10(Fab*con.c.value/1450e-10) )
                	L_tmp  = np.sort(M_1450)
        	else: L_tmp=bolometric_correction(L_bol_grid,dset_id)
        
		if return_LF[dset_name]!=None:
                	phi_fit_tmp = return_LF[dset_name](L_tmp, redshift)
                	phi_fit_pts = np.interp(L_data ,L_tmp, phi_fit_tmp)
                	PHI_data = PHI_data + (np.mean((phi_fit_pts))-np.mean((PHI_data)))

		if (len(L_data) > 0):
			L_model = bolometric_correction(L_bol_grid,dset_id)
			nu_c = c_double(dset_id)

			bolLF_model = LF_at_z(L_bol_grid,parameters,redshift,"Fiducial")
			if len(bolLF_model[np.invert(np.isfinite(bolLF_model))])!=0:
                                alldata_tem = None
                                return False
			input_c= np.power(10., bolLF_model).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
			redshift_c = c_double(redshift)
			dtg_c = c_double(return_dtg(redshift))
			res = convolve_c(input_c,nu_c,redshift_c,dtg_c)
			res = [i for i in res.contents]
			PHI_model = np.array(res,dtype=np.float64)
			#L_model, PHI_model = convolve(np.power(10.,LF_at_z(L_bol_grid,parameters,redshift,"Fiducial")), dset_id) 
			if len(PHI_model[np.invert(np.isfinite(PHI_model))])!=0: 
				alldata_tem = None
				return False
			
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
			alldata["Z"]  = np.append(alldata["Z"]  , np.ones(len(L_data)) * redshift)
			alldata["WEIGHT"]  = np.append(alldata["WEIGHT"]  , np.ones(len(L_data)) * weightfunc(redshift))
			alldata["ID"] = np.append(alldata["ID"] , np.ones(len(L_data)) * dset_id)
			
			alldata_tem["P_PRED"] = np.append(alldata_tem["P_PRED"] , phi_i)
			alldata_tem["L_OBS"]  = np.append(alldata_tem["L_OBS"]  , L_data)
			alldata_tem["P_OBS"]  = np.append(alldata_tem["P_OBS"]  , PHI_data)
			alldata_tem["D_OBS"]  = np.append(alldata_tem["D_OBS"]  , DPHI_data + 0.01)

	alldata_tem = None
	#print "NAME:",dset_name,"; redshift",redshift,";  chisq:", np.sum(((phi_i-PHI_BB)/DPHI_BB)**2)," / ",len(L_BB)

	#print "NAME:",dset_name,";  CHISQ:", np.sum(((alldata_tem["P_PRED"]-alldata_tem["P_OBS"])/alldata_tem["D_OBS"])**2)," / ",len(alldata_tem["L_OBS"])

def residual(pars):
	parameters = pars

	alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z":np.array([]),"WEIGHT":np.array([]),"ID":np.array([])}
	for key in dset_ids.keys():
		get_fit_data(alldata,parameters,zmins[key],zmaxs[key],key,dset_ids[key])

	if len(alldata["P_PRED"])==0: return np.inf, np.inf, False

	bad = np.invert(np.isfinite(alldata["P_PRED"]))
	if (np.count_nonzero(bad) > 0): alldata["P_PRED"][bad] = -40.0

	chitot = np.sum( ((alldata["P_PRED"]-alldata["P_OBS"])/alldata["D_OBS"])**2)
	print chitot, len(alldata["L_OBS"])
	residuals, sigmas = alldata["WEIGHT"]*(alldata["P_PRED"]-alldata["P_OBS"])/alldata["D_OBS"],alldata["D_OBS"]
	alldata=None	

	return residuals, sigmas, True

def lnlike(pars):
        res, sigma, success = residual(pars)
	if success==True:
        	return -0.5*np.sum(res**2 + np.log(2*np.pi*sigma**2))
	else: return -np.inf

def lnprior(pars):
	P=pars
	zdummy = 4.
	xsi = 1.+ zdummy
	gamma1=P[0]*T0(xsi)+P[1]*T1(xsi)+P[2]*T2(xsi)
	gamma2=doublepower(zdummy,[P[3],zref,P[4],P[5]])
	Phis  =P[6]*T0(xsi)+P[7]*T1(xsi)
	Lbreak=doublepower(zdummy,[P[8],zref,P[9],P[10]])
	if (np.isfinite(gamma1)) and (np.isfinite(gamma2)) and (np.isfinite(Phis)) and (np.isfinite(Lbreak)): 
		if (gamma1>-5) and (gamma1<5) and (gamma2>-5) and (gamma2<5) and (Phis>-15) and (Phis<5) and (Lbreak>5) and (Lbreak<20):
			return 0.0
		else: return -np.inf
	else: return -np.inf

def lnprob(pars):
        lp = lnprior(pars)
	llike = lnlike(pars)
        if not np.isfinite(lp):  return -np.inf
	if not np.isfinite(llike):  return -np.inf
	return lp + llike

start_time = time.time()
np.random.seed(4267)

ndim, nwalkers = 11, 100
pos = np.array([np.random.randn(ndim) for i in range(nwalkers)])
for i in range(pos.shape[0]):
	pos[i,:] = pos[i,:] * parameters_vary + parameters_init

if os.path.isfile("output/chain_end.dat")==True:
	pos = np.genfromtxt("output/chain_end.dat")[:,2:]

if os.path.isfile("output/chain.dat")==False:
	f = open("output/chain.dat", "w")
	f.close()

with MPIPool() as pool:
	if not pool.is_master():
        	pool.wait()
        	sys.exit(0)

	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool)
	nsteps = 1800
	#sampler.run_mcmc(pos, nsteps)
	print "begin sampling"
	mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
	print "initial memory consumption: ",  mem
	for i, result in enumerate(sampler.sample(pos, iterations=nsteps, storechain=False)):
        	#if pool.is_master(): 
		if (i+1) % 10 == 0:
			mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        		print "{0:5.1%}".format(float(i) / nsteps), "  t:", time.time()-start_time 
			print "memory: ", mem		

		position = result[0]
		f = open("output/chain.dat", "a")
		for k in range(position.shape[0]):
			f.write("{0:3d} {1:7d} {2:s}\n".format(k, i, np.array2string(position[k]).strip('[').strip(']').replace('\n',' ') ))
		f.close()

	# save final position for future use
	f = open("output/chain_end.dat", "w")
	for k in range(position.shape[0]):
                        f.write("{0:3d} {1:7d} {2:s}\n".format(k, i, np.array2string(position[k]).strip('[').strip(']').replace('\n',' ') ))
	f.close()

