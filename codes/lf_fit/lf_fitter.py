from data import *
import numpy as np 
from lf_shape import *
import scipy.interpolate as inter
from convolve import *
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import least_squares
# fit the luminosity function based on datasets at a given redshift
#L_BB_tmp = bolometric_correction(L_bol_grid,-1)

from lf_fitter_data import *

#parameters_init = np.array([0.10 , 	0.608, 	-6.58, 	8.427,  1.02, 
#					    -10.23, -3.00,	0.000, 	1.403,  -1.07,  
#					    0.0, 	0.00,   0.0, 	1.5, 	-3.19])		
#parameters_info = np.array(["gamma1_0", "gamma2_0", "logphis_0", "logLs_0", "k1", 
#		           "k2", "k_gamma1", "k_gamma1_2(0)" , "k_gamma2_1", "k_gamma2_2", 
#				   "nothing",  	"k3",	 "k4(0)",   "z_phis(0)", "k_phis(0)"])
parameters_init = np.array([0.417, 2.174, -4.825, 13.036, 0.632, -11.76, -14.25, -0.623, 1.460, -0.793])
parameters_info = np.array(["gamma1_0", "gamma2_0", "logphis", "logLs_0", "k1", "k2", "k3" ,"k_gamma1" , "k_gamma2_1", "k_gamma2_2"])
parameters_bound= (np.array([0,0,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf]),np.array([np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]))


def get_fit_data(alldata,parameters,zmin,zmax,dset_name,dset_id):
	z_lis = 0.5*(zmin + zmax)
	for iz in range(len(z_lis)):
		redshift = z_lis[iz]
		L_BB, PHI_BB, DPHI_BB = load_LF_data[dset_name](redshift)
		
		#L_BB_tmp=bolometric_correction(L_bol_grid,dset_id)
		#if return_LF[dset_name]!=None:
		#	phi_fit_tmp = return_LF[dset_name](L_BB_tmp, redshift)
		#	phi_fit_pts = np.interp(L_BB ,L_BB_tmp, phi_fit_tmp)
		#	PHI_BB = PHI_BB + (np.mean((phi_fit_pts))-np.mean((PHI_BB)))	
		if (len(L_BB) > 1):
			L_B, PHI_B = convolve(np.power(10.,LF_at_z(L_bol_grid,parameters,redshift,"Fiducial")), dset_id) 
			#L_B, PHI_B = convolve( 1e-3*L_bol_grid, dset_id)
			phi_i = np.interp(L_BB, L_B, np.log10(PHI_B))

			alldata["P_PRED"] = np.append(alldata["P_PRED"] , phi_i)
			alldata["L_OBS"]  = np.append(alldata["L_OBS"]  , L_BB)
			alldata["P_OBS"]  = np.append(alldata["P_OBS"]  , PHI_BB)
			alldata["D_OBS"]  = np.append(alldata["D_OBS"]  , DPHI_BB)
			alldata["Z_TOT"]  = np.append(alldata["Z_TOT"]  , np.ones(len(L_BB)) * redshift)
			#alldata["B"]      = np.append(alldata["B"]      , 0*L_BB)
			alldata["ID"]     = np.append(alldata["ID"]     , np.ones(len(L_BB)) * dset_id)

def chisq(parameters):
	alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z_TOT":np.array([]),"B":np.array([]),"ID":np.array([])}
	for key in dset_ids.keys():
		get_fit_data(alldata,parameters,zmins[key],zmaxs[key],key,dset_ids[key])

	bad = np.invert(np.isfinite(alldata["P_PRED"]))
	if (np.count_nonzero(bad) > 0): alldata["P_PRED"][bad] = -40.0

	id= (alldata["Z_TOT"]>=2.7) & (alldata["Z_TOT"]<=3.3) & (alldata["ID"]==-4)
	return alldata["L_OBS"][id],alldata["P_OBS"][id],alldata["D_OBS"][id],alldata["P_PRED"][id]
	#return np.sum(((alldata["P_PRED"]-alldata["P_OBS"])/alldata["D_OBS"])**2),len(alldata["L_OBS"])

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

#out = minimize(chisq,x0=parameters_init)#bounds=parameters_bound)
#out = least_squares(chisq,x0=parameters_init,bounds=parameters_bound)
#print out



