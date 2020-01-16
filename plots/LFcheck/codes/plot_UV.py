from data import *
import numpy as np 
from lf_shape import *
import scipy.interpolate as inter
from convolve import *
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import least_squares
# fit the luminosity function based on datasets at a given redshift
from lf_fitter_data_limited import *
from ctypes import *
import ctypes
import sys

from new_load_giallongo15_lf_shape import *

redshift=float(sys.argv[1])
dtg = return_dtg(redshift)

#load a special fit done for the Giallongo15 data
fit_res=np.genfromtxt("../../../codes/lf_fit/output/special_fit.dat",names=True)
id=fit_res["z"]==redshift
parameters_special=np.array([ fit_res["gamma1"][id],fit_res["gamma2"][id],fit_res["phi_s"][id],fit_res["L_s"][id]])

#load our global best fit at this redshift
source = np.genfromtxt("../../Fit_parameters/codes/zevolution_fit_global.dat",names=True)
zref = 2.
p=source['value'][ source['paraid']==0 ]
gamma1 = polynomial(redshift,p,2)
p=source['value'][ source['paraid']==1 ]
gamma2 = doublepower(redshift,(p[0],zref,p[1],p[2]))
p=source['value'][ source['paraid']==2 ]
logphi = polynomial(redshift,p,1)
p=source['value'][ source['paraid']==3 ]
Lbreak = doublepower(redshift,(p[0],zref,p[1],p[2]))
parameters_global = np.array([gamma1,gamma2,logphi,Lbreak])

#load the shared object file
c_extenstion = CDLL(homepath+'codes/c_lib/convolve.so')
convolve_c = c_extenstion.convolve
convolve_c.restype = ctypes.POINTER(ctypes.c_double * N_bol_grid)

#functions to get binned estimations compiled
def get_fit_data(alldata,zmin,zmax,dset_name,dset_id,rescale=False):
	alldata_tem={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([])}
	
	if dset_id!=-5: return False
	if load_LF_data[dset_name](redshift)!=False:
		L_data, PHI_data, DPHI_data = load_LF_data[dset_name](redshift)
	else: return False

        if dset_id==-5:
                L_1450 = bolometric_correction(L_bol_grid,dset_id) + L_solar
                M_1450 = -2.5*( L_1450 - np.log10(Fab*con.c.value/1450e-10) )
                L_tmp  = np.sort(M_1450)
        else: L_tmp=bolometric_correction(L_bol_grid,dset_id)

	if rescale==True:	
        	if return_LF[dset_name]!=None:
			phi_fit_tmp = return_LF[dset_name](L_tmp, redshift)
			phi_fit_pts = np.interp(L_data ,L_tmp, phi_fit_tmp)
			PHI_data = PHI_data + (np.mean((phi_fit_pts))-np.mean((PHI_data)))
	
	if len(L_data)>0:
			alldata["L_OBS"]  = np.append(alldata["L_OBS"]  , L_data)
			alldata["P_OBS"]  = np.append(alldata["P_OBS"]  , PHI_data)
			alldata["D_OBS"]  = np.append(alldata["D_OBS"]  , DPHI_data)# + 0.01)
	print "NAME:",dset_name

def get_data(rescale=False):
        alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z_TOT":np.array([]),"B":np.array([]),"ID":np.array([])}
	
        for key in dset_ids.keys():
                get_fit_data(alldata,zmins[key],zmaxs[key],key,dset_ids[key],rescale=rescale)

        return alldata["L_OBS"],alldata["P_OBS"],alldata["D_OBS"],alldata["P_PRED"]

def get_data_g15(rescale=False):
        alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z_TOT":np.array([]),"B":np.array([]),"ID":np.array([])}

        for key in dset_ids_special.keys():
                get_fit_data(alldata,zmins[key],zmaxs[key],key,dset_ids_special[key],rescale=rescale)

        return alldata["L_OBS"],alldata["P_OBS"],alldata["D_OBS"],alldata["P_PRED"]

#####
def get_fit_data_Xray(alldata,parameters,zmin,zmax,dset_name,dset_id):
        alldata_tem={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([])}
        if load_LF_data[dset_name](redshift)!=False:
                L_data, PHI_data, DPHI_data = load_LF_data[dset_name](redshift)
        else:
                return False

        L_tmp=bolometric_correction(L_bol_grid,dset_id)

        if return_LF[dset_name]!=None:
                        phi_fit_tmp = return_LF[dset_name](L_tmp, redshift)
                        phi_fit_pts = np.interp(L_data ,L_tmp, phi_fit_tmp)
                        PHI_data = PHI_data + (np.mean((phi_fit_pts))-np.mean((PHI_data)))

	if (len(L_data) > 0):
                        L_model = bolometric_correction(L_bol_grid,dset_id)
                        nu_c = c_double(dset_id)
                        redshift_c = c_double(redshift)
                        dtg_c = c_double(dtg)
                        input_c= np.power(10.,LF(L_bol_grid,parameters)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                        res = convolve_c(input_c,nu_c,redshift_c,dtg_c)
                        res = [i for i in res.contents]
                        PHI_model = np.array(res,dtype=np.float64)

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

def get_data_Xray(parameters,x,y,dataid=-4):
        alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z_TOT":np.array([]),"B":np.array([]),"ID":np.array([])}
        for key in dset_ids.keys():
                get_fit_data_Xray(alldata,parameters,zmins[key],zmaxs[key],key,dset_ids[key])

        bad = np.invert(np.isfinite(alldata["P_PRED"]))
        if (np.count_nonzero(bad) > 0): alldata["P_PRED"][bad] = -40.0

	select = alldata["ID"]==dataid
	logLbol = 0.0*alldata["L_OBS"][select]
	for i,Lband in enumerate(alldata["L_OBS"][select]):
		logLbol[i] = bolometric_correction_inverse(Lband,dataid)

	L1450 = bolometric_correction(logLbol,-5)

	M1450 =  -2.5*( L1450 + L_solar - np.log10(Fab*(con.c.value/1450e-10)))
	
	y = y[np.argsort(x)]
	x = np.sort(x)

	phi_fit_pts = np.interp(M1450, x, y)
        return M1450, phi_fit_pts+(alldata["P_OBS"][select]-alldata["P_PRED"][select]),alldata["D_OBS"][select]


#####
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

'''
L_1450 = bolometric_correction(L_bol_grid,-5)
nu_c = c_double(-5)
redshift_c = c_double(redshift)
dtg_c = c_double(dtg)
input_c= np.power(10.,LF(L_bol_grid,parameters_special)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
res = convolve_c(input_c,nu_c,redshift_c,dtg_c)
res = [i for i in res.contents]
PHI_1450 = np.array(res,dtype=np.float64)
x = -2.5*( L_1450 + L_solar - np.log10(Fab*(con.c.value/1450e-10)) )
y = np.log10(PHI_1450) - np.log10(2.5)
ax.plot(x,y,'-',lw=4,c='gold',label=r'$\rm Special$ $\rm fit$')
'''

#plot the predicted QLF in the UV
L_1450 = bolometric_correction(L_bol_grid,-5)
nu_c = c_double(-5)
redshift_c = c_double(redshift)
dtg_c = c_double(dtg)
input_c= np.power(10.,LF(L_bol_grid,parameters_global)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
res = convolve_c(input_c,nu_c,redshift_c,dtg_c)
res = [i for i in res.contents]
PHI_1450 = np.array(res,dtype=np.float64)  
x = -2.5*( L_1450 + L_solar - np.log10(Fab*(con.c.value/1450e-10)) )  # convert lum to mag
xm = x.copy()
ym = np.log10(PHI_1450) - np.log10(2.5) # convert to number density per mag
ax.plot(xm,ym,'-',lw=4,c='darkorchid',label=r'$\rm Global$ $\rm fit$')

y = return_kk18_lf_fitted(x,redshift) 
ax.plot(x,y,'--',dashes=(25,15),c='seagreen',label=r'$\rm Kulkarni+$ $\rm 2018$')

'''
x,y,dy,yfit=get_data()
#ax.errorbar(x,y,yerr=dy,capsize=6,linestyle='',lw=2,c='gray',mec='gray',marker='x', ms=10, capthick=2 ,label=r'$\rm UV$ ($\rm not$ $\rm rescaled$)')
ax.plot(x,y,linestyle='',lw=3,mew=2,c='gray',mec='gray',marker='x',ms=15,label=r'$\rm UV$ ($\rm original$)')
'''

x,y,dy,yfit=get_data_g15()
ax.plot(x,y,linestyle='',lw=3,mew=2,c='chocolate',mec='chocolate',marker='x',ms=15)

x,y,dy,yfit=get_data(rescale=True)
ax.errorbar(x,y,yerr=dy,capsize=6,linestyle='',lw=2,c='black',mec='black',alpha=0.8,marker='o', ms=10, capthick=2 ,label=r'$\rm UV$')

x,y,dy=get_data_Xray(parameters_global,xm,ym)
ax.errorbar(x,y,yerr=dy,capsize=9,linestyle='',lw=2,c='royalblue',mec='royalblue',marker='*', alpha=0.8, ms=15, capthick=3 ,label=r'$\rm Hard$ $\rm X-ray$')

if (redshift<6.5) and (redshift>=4):
	x = np.linspace(-14,-32,100) 
	y = return_giallongo15_lf_fitted(x,redshift)
	ax.plot(x,y,'--',dashes=(25,15),c='chocolate',label=r'$\rm Giallongo+$ $\rm 2015$')

prop = matplotlib.font_manager.FontProperties(size=25.0)
if redshift == 4.2:
	ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
#ax.set_xlabel(r'$\log{(L_{\rm B}/{\rm L}_{\odot})}$',fontsize=40,labelpad=2.5)
ax.set_xlabel(r'$M_{\rm 1450}$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\phi[{\rm mag}^{-1}{\rm cMpc}^{-3}])}$',fontsize=40,labelpad=5)
ax.text(0.88, 0.92, r'${\rm z\sim'+str(redshift)+'}$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=40)

ax.set_xlim(-17.0,-30.5)
ax.set_ylim(-12.0,-4.)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/tension_"+str(redshift)+".pdf",fmt='pdf')
#plt.show()

