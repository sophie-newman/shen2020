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

bestfit2=np.array([0.39856372,2.19426155,-4.73290265,12.97241616,0.4371788,-11.63470663,-11.7235743,-0.72820353,1.36234503,-0.79697647])
bestfit=np.array([0.38181256,2.16955741,-4.70557406,12.94851374,0.43771614,-11.42561263,-11.34952214,-0.7552896,1.32130027,-0.77768681])
parameters_init = np.array([0.41698725,2.1744386,-4.8250643,13.035753,0.63150872,-11.763560,-14.249833,-0.62298947,1.4599393,-0.79280099])
parameters_info = np.array(["gamma1_0", "gamma2_0", "logphis", "logLs_0", "k1", "k2", "k3" ,"k_gamma1" , "k_gamma2_1", "k_gamma2_2"])
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
		L_BB, PHI_BB, DPHI_BB = load_LF_data[dset_name](redshift)
		
		L_BB_tmp=bolometric_correction(L_bol_grid,dset_id)
		if return_LF[dset_name]!=None:
			phi_fit_tmp = return_LF[dset_name](L_BB_tmp, redshift)
			phi_fit_pts = np.interp(L_BB ,L_BB_tmp, phi_fit_tmp)
			PHI_BB = PHI_BB + (np.mean((phi_fit_pts))-np.mean((PHI_BB)))	
		if (len(L_BB) > 0):
			L_B = bolometric_correction(L_bol_grid,dset_id)
			nu_c = c_double(dset_id)
			input_c= np.power(10.,LF_at_z(L_bol_grid,parameters,redshift,"Fiducial")).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
			res = convolve_c(input_c,nu_c)
			res = [i for i in res.contents]
			PHI_B = np.array(res,dtype=np.float64)
			#L_B, PHI_B = convolve(np.power(10.,LF_at_z(L_bol_grid,parameters,redshift,"Fiducial")), dset_id) 
			
			phi_i = np.interp(L_BB, L_B, np.log10(PHI_B))

			alldata["P_PRED"] = np.append(alldata["P_PRED"] , phi_i)
			alldata["L_OBS"]  = np.append(alldata["L_OBS"]  , L_BB)
			alldata["P_OBS"]  = np.append(alldata["P_OBS"]  , PHI_BB)
			alldata["D_OBS"]  = np.append(alldata["D_OBS"]  , DPHI_BB)# + 0.01)
			alldata["Z_TOT"]  = np.append(alldata["Z_TOT"]  , np.ones(len(L_BB)) * redshift)
			#alldata["B"]      = np.append(alldata["B"]      , 0*L_BB)
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

	#print "call once"
	#id= (alldata["Z_TOT"]>=2.7) & (alldata["Z_TOT"]<=3.3) & (alldata["ID"]==-4)
	#return alldata["L_OBS"][id],alldata["P_OBS"][id],alldata["D_OBS"][id],alldata["P_PRED"][id]
	return np.sum(((alldata["P_PRED"]-alldata["P_OBS"])/alldata["D_OBS"])**2)#,len(alldata["L_OBS"])-10

t0=time.time()
print chisq(parameters_init)
print time.time()-t0
#exit()

print chisq(bestfit)
print chisq(bestfit2)
'''
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
x,y,dy,yfit=chisq(parameters_init)
ax.errorbar(x,y,yerr=dy,capsize=10,linestyle='',c='crimson',marker='o',markeredgewidth=0, ms=12,alpha=1,label=r'$\rm{Observations}$')
ax.plot(x,yfit,'.',c='royalblue',marker='o',ms=12,label=r'$\rm Best-fit$')
prop = matplotlib.font_manager.FontProperties(size=30.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\log{(L_{\rm HX}/{\rm L}_{\odot})}$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\Phi[{\rm dex}^{-1}{\rm Mpc}^{-3}])}$',fontsize=40,labelpad=5)
ax.text(0.85, 0.92, r'${z=2.7-3.3}$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=40)

#ax.set_xlim(1.5,10.5)
#ax.set_ylim(-2.9,-1.3)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.show()
'''
#out = minimize(chisq,x0=parameters_init,method='Nelder-Mead')#bounds=parameters_bound)
#print out



