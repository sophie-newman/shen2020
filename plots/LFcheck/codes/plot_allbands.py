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
#dtg=float(sys.argv[2])
dtg=return_dtg(redshift)

# parameters of the H07 model
parameters_init = np.array([0.41698725, 2.17443860, -4.82506430, 13.03575300, 0.63150872, -11.76356000, -14.24983300, -0.62298947, 1.45993930, -0.79280099])

# our best-fits
fit_res=np.genfromtxt("../../../codes/lf_fit/output/fit_at_z_nofix.dat",names=True)
id=fit_res["z"]==redshift
parameters_free_local=np.array([ fit_res["gamma1"][id],fit_res["gamma2"][id],fit_res["phi_s"][id],fit_res["L_s"][id]])

fit_res=np.genfromtxt("../../../codes/lf_fit/output/fit_at_z_fix.dat",names=True)
id=fit_res["z"]==redshift
parameters_fix_local=np.array([ fit_res["gamma1"][id],fit_res["gamma2"][id],fit_res["phi_s"][id],fit_res["L_s"][id]])

fit_evolve=np.genfromtxt("../../Fit_parameters/codes/zevolution_fit.dat",names=['gamma1','gamma2','phis','Lbreak'])
parameters_global_1 = pars_at_z(fit_evolve,redshift)

source = np.genfromtxt("../../Fit_parameters/codes/zevolution_fit_global.dat",names=True)
zref = 2.
p=source['value'][ source['paraid']==0 ]
gamma1 = polynomial(redshift,p,2)
p=source['value'][ source['paraid']==1 ]
gamma2 = doublepower(redshift,(p[0],zref, p[1], p[2]))
p=source['value'][ source['paraid']==2 ]
logphi = polynomial(redshift,p,1) 
p=source['value'][ source['paraid']==3 ]
Lbreak = doublepower(redshift,(p[0],zref, p[1], p[2]))
parameters_global_2 = np.array([gamma1,gamma2,logphi,Lbreak])

############################
# parameters of the local QLF
p=source['value'][ source['paraid']==0 ]
gamma1 = polynomial(0.1,p,2)
p=source['value'][ source['paraid']==1 ]
gamma2 = doublepower(0.1,(p[0],zref, p[1], p[2]))
p=source['value'][ source['paraid']==2 ]
logphi = polynomial(0.1,p,1)
p=source['value'][ source['paraid']==3 ]
Lbreak = doublepower(0.1,(p[0],zref, p[1], p[2]))
parameters_global_z0 = np.array([gamma1,gamma2,logphi,Lbreak])

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
	
        if dset_id==-5:
                L_1450 = bolometric_correction(L_bol_grid,dset_id) + L_solar
                M_1450 = -2.5*( L_1450 - np.log10(Fab*con.c.value/1450e-10) ) 
                L_tmp  = np.sort(M_1450)
        else: L_tmp=bolometric_correction(L_bol_grid,dset_id)

        if return_LF[dset_name]!=None:
		#if (redshift!=0.4) or (dset_id!=-2):
		#if (dset_id!=-2):# and (dset_id==-5):
                	phi_fit_tmp = return_LF[dset_name](L_tmp, redshift)
                	phi_fit_pts = np.interp(L_data ,L_tmp, phi_fit_tmp)
                	PHI_data = PHI_data + (np.mean((phi_fit_pts))-np.mean((PHI_data)))		
	#if len(PHI_data)==0:
	#	print dset_name
	#if dset_id == -5:
	#	print dset_name, L_data

        if (len(L_data) > 0):
                        L_model = bolometric_correction(L_bol_grid,dset_id)
                        nu_c = c_double(dset_id)
			#if dset_id==-5:
			#	L_model = bolometric_correction(L_bol_grid,-1)
                        #	nu_c = c_double(-1)
			redshift_c = c_double(redshift)
			dtg_c = c_double(dtg)
                        input_c= np.power(10.,LF(L_bol_grid,parameters)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                        res = convolve_c(input_c,nu_c,redshift_c,dtg_c)
                        res = [i for i in res.contents]
                        PHI_model = np.array(res,dtype=np.float64)
                        #L_model, PHI_model = convolve(np.power(10.,LF_at_z(L_bol_grid,parameters,redshift,"Fiducial")), dset_id)
			
			if dset_id==-5:
				#L_Bband = (M_sun_Bband_AB-(L_data - 0.706))/2.5
                                #phi_i = np.interp(L_Bband, L_model, np.log10(PHI_model))

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

def get_data(parameters,dataid):
        alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z_TOT":np.array([]),"B":np.array([]),"ID":np.array([])}
        for key in dset_ids.keys():
                get_fit_data(alldata,parameters,zmins[key],zmaxs[key],key,dset_ids[key])

        bad = np.invert(np.isfinite(alldata["P_PRED"]))
        if (np.count_nonzero(bad) > 0): alldata["P_PRED"][bad] = -40.0

	if (dataid!=-5):
		select = alldata["ID"]==dataid
		logLbol = 0.0*alldata["L_OBS"][select]
		for i,Lband in enumerate(alldata["L_OBS"][select]):
			logLbol[i] = bolometric_correction_inverse(Lband,dataid)
	elif (dataid==-5):
		select = alldata["ID"]==dataid
		L_1450 = -0.4*alldata["L_OBS"][select] + np.log10(Fab*(con.c.value/1450e-10)) - L_solar
		logLbol = 0.0*L_1450
		for i,Lband in enumerate(L_1450):
			logLbol[i] = bolometric_correction_inverse(Lband,dataid)

	x = L_bol_grid
	y = LF(L_bol_grid,parameters)
	phi_fit_pts = np.interp(logLbol, x, y)
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
#y = LF(L_bol_grid,parameters_free_local)
#ax.plot(x,y,'--',dashes=(25,15),c='navy',alpha=0.7,label=r'$\rm Local$ $\rm fit$ ($\rm free$)')
y = LF(L_bol_grid,parameters_fix_local)
ax.plot(x,y,'--',dashes=(25,15),c='chocolate',alpha=0.7,label=r'$\rm Local$ $\rm fit$ ($\phi_{\ast}(z)$ $\rm fixed$)')
#y = LF(L_bol_grid,parameters_global_1)
#ax.plot(x,y,':',c='tan',label=r'$\rm Fit$ $\rm on$ $\rm local$ $\rm fits$')
y = LF(L_bol_grid,parameters_global_2)
ax.plot(x,y,'--',dashes=(25,15),c='darkorchid',label=r'$\rm Global$ $\rm fit$')

y = LF(L_bol_grid,parameters_global_z0)
ax.plot(x,y,'--',dashes=(25,15),lw=2,c='cyan',label=r'$\rm Global$ $\rm fit$ ($\rm z\sim0$)')

x = L_bol_grid + L_solar 
y = LF_at_z_H07(L_bol_grid,parameters_init,redshift,"Fiducial")
ax.plot(x,y,'--',dashes=(25,15),c='crimson',label=r'$\rm Hopkins+$ $\rm 2007$')

xcollect = np.array([])

x,y,dy=get_data(parameters_fix_local,dataid=-5)
ax.errorbar(x,y,yerr=dy,linestyle='none',c='seagreen',mec='seagreen',marker='o',ms=10,capsize=6,capthick=2,lw=2,label=r'$\rm UV$ $\rm 1450\AA$')
xcollect = np.append(xcollect,x)

x,y,dy=get_data(parameters_fix_local,dataid=-1)
ax.errorbar(x,y,yerr=dy,linestyle='none',c='pink',mec='pink',marker='o',ms=10,capsize=6,capthick=2,lw=2,label=r'$\rm B$ $\rm Band$')
xcollect = np.append(xcollect,x)

x,y,dy=get_data(parameters_fix_local,dataid=-4)
ax.errorbar(x,y,yerr=dy,linestyle='none',c='royalblue',mec='royalblue',marker='o',ms=10,capsize=6,capthick=2,lw=2,label=r'$\rm Hard$ $\rm X-ray$')
xcollect = np.append(xcollect,x)

x,y,dy=get_data(parameters_fix_local,dataid=-3)
ax.errorbar(x,y,yerr=dy,linestyle='none',c='gray',mec='gray',marker='o',ms=10,capsize=6,capthick=2,lw=2,label=r'$\rm Soft$ $\rm X-ray$')
xcollect = np.append(xcollect,x)

x,y,dy=get_data(parameters_fix_local,dataid=-2)
ax.errorbar(x,y,yerr=dy,linestyle='none',c='olive',mec='olive',marker='o',ms=10,capsize=6,capthick=2,lw=2,label=r'$\rm Mid$ $\rm IR$')
xcollect = np.append(xcollect,x)

xcollect = np.sort(xcollect-L_solar)
print redshift, xcollect[2], xcollect[4], len(xcollect)

ax.axvline(parameters_fix_local[3]+L_solar,color='gold',lw=3,alpha=1)
ax.axhline(parameters_fix_local[2]-np.log10(2.),color='gold',lw=3,alpha=1)

prop = matplotlib.font_manager.FontProperties(size=22.0)
if not (redshift in [2,3,4,5,6]):
	ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\log{(L_{\rm bol}[{\rm erg}\,{\rm s}^{-1}])}$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\phi[{\rm dex}^{-1}{\rm cMpc}^{-3}])}$',fontsize=40,labelpad=5)

ax.text(0.88, 0.92, r'${\rm z\sim'+str(redshift)+'}$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=40)

#ax.set_xlim(42.5,51.3)
ax.set_xlim(42.5,50.3)
ax.set_ylim(-11.2,-2.3)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/bol_"+str(redshift)+".pdf",fmt='pdf')
#plt.show()

