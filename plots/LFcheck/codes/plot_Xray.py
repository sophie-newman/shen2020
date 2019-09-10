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
import sys

redshift=float(sys.argv[1])

fit_res=np.genfromtxt("../../../codes/lf_fit/output/fit_at_z_fix.dat",names=True)
id=fit_res["z"]==redshift
parameters=np.array([ fit_res["gamma1"][id],fit_res["gamma2"][id],fit_res["phi_s"][id],fit_res["L_s"][id]])

#load the shared object file
c_extenstion = CDLL(homepath+'codes/c_lib/convolve.so')
convolve_c = c_extenstion.convolve
convolve_c.restype = ctypes.POINTER(ctypes.c_double * N_bol_grid)

def get_fit_data(alldata,zmin,zmax,dset_name,dset_id):
	alldata_tem={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([])}
	
	if dset_id!=-4: return False
	
	if load_LF_data[dset_name](redshift)!=False:
		L_data, PHI_data, DPHI_data = load_LF_data[dset_name](redshift)
	else: return False
	
	L_tmp=bolometric_correction(L_bol_grid,dset_id)
	if (return_LF[dset_name]!=None):
		phi_fit_tmp = return_LF[dset_name](L_tmp, redshift)
		phi_fit_pts = np.interp(L_data ,L_tmp, phi_fit_tmp)
		PHI_data = PHI_data + (np.mean((phi_fit_pts))-np.mean((PHI_data)))	
	
	if len(L_data)>0:
			alldata["L_OBS"]  = np.append(alldata["L_OBS"]  , L_data)
			alldata["P_OBS"]  = np.append(alldata["P_OBS"]  , PHI_data)
			alldata["D_OBS"]  = np.append(alldata["D_OBS"]  , DPHI_data)# + 0.01)
	print "NAME:",dset_name

def get_data():
        alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z_TOT":np.array([]),"B":np.array([]),"ID":np.array([])}
        for key in dset_ids.keys():
                get_fit_data(alldata,zmins[key],zmaxs[key],key,dset_ids[key])

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

L_HX = bolometric_correction(L_bol_grid,-4)
nu_c = c_double(-4)
redshift_c = c_double(redshift)
dtg_c = c_double(return_dtg(redshift))
input_c= np.power(10.,LF(L_bol_grid,parameters)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
res = convolve_c(input_c,nu_c,redshift_c,dtg_c)
res = [i for i in res.contents]
PHI_HX = np.array(res,dtype=np.float64)
x = L_HX + L_solar
y = np.log10(PHI_HX)
ax.plot(x,y,'--',dashes=(25,15),c='black',label=r'$\rm new$ $\rm fit$')

x,y,dy,yfit=get_data()
x = x + L_solar
y = y

ax.errorbar(x,y,yerr=dy,capsize=6,linestyle='',lw=2,c='crimson',mec='crimson',marker='o', ms=10,capthick=2,label=r'$\rm data$')

prop = matplotlib.font_manager.FontProperties(size=30.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\log{(L_{\rm HX}[{\rm erg}\,{\rm s}^{-1}])}$',fontsize=40,labelpad=2.5)
#ax.set_xlabel(r'$M_{\rm 1450}$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\phi[{\rm dex}^{-1}{\rm Mpc}^{-3}])}$',fontsize=40,labelpad=5)
ax.text(0.88, 0.92, r'${\rm z\sim'+str(redshift)+'}$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=40)

ax.set_xlim(41.7, 47.7)
ax.set_ylim(-10.2,-2.3)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
#plt.savefig("../figs/HX_"+str(redshift)+".pdf",fmt='pdf')
plt.show()

