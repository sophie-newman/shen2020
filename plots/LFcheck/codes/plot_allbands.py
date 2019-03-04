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

#bestfit        = np.array([0.39856372, 2.19426155, -4.73290265, 12.97241616, 0.43717880, -11.63470663, -11.72357430, -0.72820353, 1.36234503, -0.79697647])
#bestfit         = np.array([0.38181256, 2.16955741, -4.70557406, 12.94851374, 0.43771614, -11.42561263, -11.34952214, -0.75528960, 1.32130027, -0.77768681])
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
	#if len(PHI_data)==0:
	#	print dset_name

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

def get_data(parameters,dataid):
        alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z_TOT":np.array([]),"B":np.array([]),"ID":np.array([])}
        for key in dset_ids.keys():
                get_fit_data(alldata,parameters,zmins[key],zmaxs[key],key,dset_ids[key])

        bad = np.invert(np.isfinite(alldata["P_PRED"]))
        if (np.count_nonzero(bad) > 0): alldata["P_PRED"][bad] = -40.0

	if (dataid!=-1) and  (dataid!=-1.1):
		select = alldata["ID"]==dataid
		logLbol = 0.0*alldata["L_OBS"][select]
		for i,Lband in enumerate(alldata["L_OBS"][select]):
			logLbol[i] = bolometric_correction_inverse(Lband,dataid)
	elif (dataid==-1.1):
		select = alldata["ID"]==dataid
                logLbol = 0.0*alldata["L_OBS"][select]
		for i,Lband in enumerate(alldata["L_OBS"][select]):
                        logLbol[i] = bolometric_correction_inverse(Lband,-1)
	elif (dataid==-1):
		select = alldata["ID"]==dataid
		L_Bband = (M_sun_Bband_AB-(alldata["L_OBS"][select] - 0.706))/2.5
		logLbol = 0.0*L_Bband
		for i,Lband in enumerate(L_Bband):
			logLbol[i] = bolometric_correction_inverse(Lband,dataid)

	x = L_bol_grid
	y = LF(L_bol_grid,parameters)
	phi_fit_pts = np.interp(logLbol , x, y)

        return logLbol+L_solar, phi_fit_pts+(alldata["P_OBS"][select]-alldata["P_PRED"][select]),alldata["D_OBS"][select]


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

x = L_bol_grid + L_solar 
y = LF(L_bol_grid,parameters)
ax.plot(x,y,'--',dashes=(25,15),c='black',label=r'$\rm new$ $\rm fit$')

x = L_bol_grid + L_solar 
y = LF_at_z_H07(L_bol_grid,parameters_init,redshift,"Fiducial")
ax.plot(x,y,'--',dashes=(25,15),c='gray',label=r'$\rm old$ $\rm fit$')

x,y,dy=get_data(parameters,dataid=-1)
ax.errorbar(x,y,yerr=dy,linestyle='none',c='crimson',mec='crimson',marker='o',ms=10,capsize=6,capthick=2,label=r'$\rm UV$ $\rm 1450\AA$')
x,y,dy=get_data(parameters,dataid=-1.1)
ax.errorbar(x,y,yerr=dy,linestyle='none',c='crimson',mec='crimson',marker='o',ms=10,capsize=6,capthick=2)

x,y,dy=get_data(parameters,dataid=-4)
ax.errorbar(x,y,yerr=dy,linestyle='none',c='royalblue',mec='royalblue',marker='o',ms=10,capsize=6,capthick=2,label=r'$\rm Hard$ $\rm Xray$')

ax.axvline(parameters[3]+L_solar,color='cyan',lw=2)
ax.axhline(parameters[2]-np.log10(2.),color='cyan',lw=2)

prop = matplotlib.font_manager.FontProperties(size=30.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\log{(L_{\rm bol}[{\rm erg}\,{\rm s}^{-1}])}$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\phi[{\rm dex}^{-1}{\rm Mpc}^{-3}])}$',fontsize=40,labelpad=5)
ax.text(0.88, 0.92, r'${\rm z\sim'+str(redshift)+'}$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=40)

ax.set_xlim(42.5,51.3)
ax.set_ylim(-11.2,-2.3)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/bol_"+str(redshift)+".pdf",fmt='pdf')
#plt.show()

