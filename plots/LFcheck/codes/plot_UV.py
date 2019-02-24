from data import *
import numpy as np 
from lf_shape import *
import scipy.interpolate as inter
from convolve import *
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import least_squares
# fit the luminosity function based on datasets at a given redshift
from lf_data_compilation import *
from ctypes import *
import ctypes
import sys

redshift=float(sys.argv[1])

parameters_init = np.array([0.41698725, 2.17443860, -4.82506430, 13.03575300, 0.63150872, -11.76356000, -14.24983300, -0.62298947, 1.45993930, -0.79280099])
#parameters_info = np.array(["gamma1_0", "gamma2_0", "logphis"  , "logLs_0"  , "k1"      , "k2"        , "k3"        , "k_gamma1" , "k_gamma2_1", "k_gamma2_2"])
#parameters_bound= (np.array([0,0,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf]),np.array([np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]))

fit_res=np.genfromtxt("../../fitresult/fit_at_z.dat",names=True)
id=fit_res["z"]==redshift
parameters=np.array([ fit_res["gamma1"][id],fit_res["gamma2"][id],fit_res["phi_s"][id],fit_res["L_s"][id]])

#load the shared object file
c_extenstion = CDLL(homepath+'codes/c_lib/convolve.so')
convolve_c = c_extenstion.convolve
convolve_c.restype = ctypes.POINTER(ctypes.c_double * N_bol_grid)

def get_fit_data(alldata,parameters,zmin,zmax,dset_name,dset_id,newdata=False):
	alldata_tem={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([])}
	
	if dset_id!=-1: return False
	if newdata==False: 
		if load_LF_data[dset_name](redshift)!=False:
			L_data, PHI_data, DPHI_data = load_LF_data[dset_name](redshift)
		else: return False
	else: 
		if new_load_LF_data[dset_name](redshift)!=False:
			L_data, PHI_data, DPHI_data = new_load_LF_data[dset_name](redshift)
		else: return False

	L_tmp=bolometric_correction(L_bol_grid,dset_id)
	if newdata==False:
		if (return_LF[dset_name]!=None):
			phi_fit_tmp = return_LF[dset_name](L_tmp, redshift)
			phi_fit_pts = np.interp(L_data ,L_tmp, phi_fit_tmp)
			PHI_data = PHI_data + (np.mean((phi_fit_pts))-np.mean((PHI_data)))	
	else:
		M_1450_tmp = (M_sun_Bband_AB -2.5*L_tmp) + 0.706
		M_1450_tmp = np.sort(M_1450_tmp)
		phi_fit_tmp = return_kk18_lf_fitted(M_1450_tmp, redshift)
		phi_fit_pts = np.interp(L_data ,M_1450_tmp, phi_fit_tmp)
		PHI_data = PHI_data + (np.mean((phi_fit_pts))-np.mean((PHI_data)))

	if len(L_data)>0:
			alldata["L_OBS"]  = np.append(alldata["L_OBS"]  , L_data)
			alldata["P_OBS"]  = np.append(alldata["P_OBS"]  , PHI_data)
			alldata["D_OBS"]  = np.append(alldata["D_OBS"]  , DPHI_data)# + 0.01)
	print "NAME:",dset_name

def get_data(parameters=parameters_init,newdata=False):
        alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z_TOT":np.array([]),"B":np.array([]),"ID":np.array([])}
	if newdata==False:
        	for key in dset_ids.keys():
                	get_fit_data(alldata,parameters,zmins[key],zmaxs[key],key,dset_ids[key])

        	return alldata["L_OBS"],alldata["P_OBS"],alldata["D_OBS"],alldata["P_PRED"]
	if newdata==True:
		for key in new_dset_ids.keys():
                        get_fit_data(alldata,parameters,new_zmins[key],new_zmaxs[key],key,new_dset_ids[key],newdata=True)

		return alldata["L_OBS"],alldata["P_OBS"],alldata["D_OBS"],alldata["P_PRED"]

import matplotlib.pyplot as plt 
import matplotlib

matplotlib.style.use('classic')
matplotlib.rc('xtick.major', size=15, width=3)
matplotlib.rc('xtick.minor', size=7.5, width=3)
matplotlib.rc('ytick.major', size=15, width=3)
matplotlib.rc('ytick.minor', size=7.5, width=3)
matplotlib.rc('lines',linewidth=4)
matplotlib.rc('axes', linewidth=4)

fig=plt.figure(figsize = (15,10))
ax = fig.add_axes([0.13,0.12,0.79,0.83])

L_B = bolometric_correction(L_bol_grid,-1)
nu_c = c_double(-1)
input_c= np.power(10.,LF(L_bol_grid,parameters)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
res = convolve_c(input_c,nu_c)
res = [i for i in res.contents]
PHI_B = np.array(res,dtype=np.float64)
x = (M_sun_Bband_AB -2.5*L_B) + 0.706
y = np.log10(PHI_B) - np.log10(2.5)
ax.plot(x,y,'--',dashes=(25,15),c='black',label=r'$\rm new$ $\rm fit$')

L_B = bolometric_correction(L_bol_grid,-1)
nu_c = c_double(-1)
input_c= np.power(10.,LF_at_z_H07(L_bol_grid,parameters_init,redshift,"Fiducial")).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
res = convolve_c(input_c,nu_c)
res = [i for i in res.contents]
PHI_B = np.array(res,dtype=np.float64)
x = (M_sun_Bband_AB -2.5*L_B) + 0.706
y = np.log10(PHI_B) - np.log10(2.5)
ax.plot(x,y,'--',dashes=(25,15),c='gray',label=r'$\rm old$ $\rm fit$')

x,y,dy,yfit=get_data(newdata=True)
ax.errorbar(x,y,yerr=dy,capsize=10,linestyle='',c='crimson',marker='o',markeredgewidth=0, ms=10,alpha=0.6,label=r'$\rm new$ $\rm data$')

x,y,dy,yfit=get_data()
x = (M_sun_Bband_AB -2.5*x) + 0.706
y = y - np.log10(2.5)
ax.errorbar(x,y,yerr=dy,capsize=10,linestyle='',c='royalblue',marker='o',markeredgewidth=0, ms=10,alpha=0.6,label=r'$\rm old$ $\rm data$')

prop = matplotlib.font_manager.FontProperties(size=30.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
#ax.set_xlabel(r'$\log{(L_{\rm B}/{\rm L}_{\odot})}$',fontsize=40,labelpad=2.5)
ax.set_xlabel(r'$M_{\rm 1450}$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\phi[{\rm mag}^{-1}{\rm Mpc}^{-3}])}$',fontsize=40,labelpad=5)
ax.text(0.88, 0.92, r'${\rm z\sim'+str(redshift)+'}$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=40)

ax.set_xlim(-19.5,-32)
ax.set_ylim(-10.2,-4.3)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
#plt.savefig("../figs/UV_"+str(redshift)+".pdf",fmt='pdf')
plt.show()

