from data import *
import numpy as np 
from lf_shape import *
import scipy.interpolate as inter
from convolve import *
from scipy.optimize import minimize
from scipy.optimize import leastsq
# fit the luminosity function based on datasets at a given redshift
from lf_fitter_data import *
import matplotlib.pyplot as plt

def weightfunc(z):
        if z<3.: return 1.
        elif z<4: return np.power( (1.+z)/(1.+3.), 1.)
        else: return np.power( (1.+4.)/(1.+3.), 1.) * np.power( (1.+z)/(1.+4.), 1.5)

def get_data(alldata,zmin,zmax,dset_name,dset_id):
        z_lis = 0.5*(zmin + zmax)
        for iz in range(len(z_lis)):
                redshift = z_lis[iz]
                L_data, PHI_data, DPHI_data = load_LF_data[dset_name](redshift)

                alldata["L_OBS"]  = np.append(alldata["L_OBS"]  , L_data)
                alldata["P_OBS"]  = np.append(alldata["P_OBS"]  , PHI_data)
                alldata["Z"]  = np.append(alldata["Z"]  , np.ones(len(L_data)) * redshift)
                alldata["WEIGHT"]  = np.append(alldata["WEIGHT"]  , np.ones(len(L_data)) * weightfunc(redshift))
		return len(L_data)	

def loop_over_all_dset():
	alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z":np.array([]),"WEIGHT":np.array([]),"ID":np.array([])}
	for key in dset_ids.keys():
		#if dset_ids[key]==-5:
			get_data(alldata,zmins[key],zmaxs[key],key,dset_ids[key])
			print key

	return alldata["Z"], alldata["WEIGHT"]

bins = np.linspace(0,7,11)
z,w=loop_over_all_dset()
plt.hist(z,bins=bins,weights=w,color='b')
plt.hist(z,bins=bins,color='r',alpha=0.3)
plt.show()



