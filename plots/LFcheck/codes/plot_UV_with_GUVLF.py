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

from new_load_giallongo15_lf_shape import *

obdataPath = '../GUVLF_data/'
def plt_obdata(fname,papername,color,label=True):
	obdata=np.genfromtxt(obdataPath+fname,names=True, comments='#')
	x_ob=obdata['m']
	phi=obdata['phi']
	id1=phi>0
	id2= (phi<=0)
	y_ob=0.0*x_ob
	y_ob[id1]=np.log10(phi[id1])
	y_ob[id2]=phi[id2]
	uperr=0.0*x_ob
	uperr[id1]=np.log10(phi[id1]+obdata['uperr'][id1])-y_ob[id1]
	uperr[id2]=obdata['uperr'][id2]
	low=phi-obdata['lowerr']
	low[low<=0]=1e-20
	lowerr=0.0*x_ob
	lowerr[id1]=-np.log10(low[id1])+y_ob[id1]
	lowerr[id2]=obdata['lowerr'][id2]
	#if label==True:
	#	ax.errorbar(x_ob,y_ob,yerr=(lowerr,uperr),c=color,lw=3,linestyle='',marker='o',markersize=9,capsize=4.5,label=papername)
	#else: 
	ax.errorbar(x_ob,y_ob,yerr=(lowerr,uperr),c=color,mec=color,linestyle='none',marker='o',lw=2,markersize=10,capsize=6,capthick=2)


redshift=float(sys.argv[1])
dtg = return_dtg(redshift)

parameters_init = np.array([0.41698725, 2.17443860, -4.82506430, 13.03575300, 0.63150872, -11.76356000, -14.24983300, -0.62298947, 1.45993930, -0.79280099])

fit_res=np.genfromtxt("../../../codes/lf_fit/output/fit_at_z_nofix.dat",names=True)
id=fit_res["z"]==redshift
parameters=np.array([ fit_res["gamma1"][id],fit_res["gamma2"][id],fit_res["phi_s"][id],fit_res["L_s"][id]])

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

def get_fit_data(alldata,parameters,zmin,zmax,dset_name,dset_id):
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

	'''
        if return_LF[dset_name]!=None:
		phi_fit_tmp = return_LF[dset_name](L_tmp, redshift)
		phi_fit_pts = np.interp(L_data ,L_tmp, phi_fit_tmp)
		PHI_data = PHI_data + (np.mean((phi_fit_pts))-np.mean((PHI_data)))
	'''
	if len(L_data)>0:
			alldata["L_OBS"]  = np.append(alldata["L_OBS"]  , L_data)
			alldata["P_OBS"]  = np.append(alldata["P_OBS"]  , PHI_data)
			alldata["D_OBS"]  = np.append(alldata["D_OBS"]  , DPHI_data)# + 0.01)
	print "NAME:",dset_name

def get_data(parameters=parameters_init):
        alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z_TOT":np.array([]),"B":np.array([]),"ID":np.array([])}
        for key in dset_ids.keys():
                get_fit_data(alldata,parameters,zmins[key],zmaxs[key],key,dset_ids[key])

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

L_1450 = bolometric_correction(L_bol_grid,-5)

nu_c = c_double(-5)
redshift_c = c_double(redshift)
dtg_c = c_double(dtg)
input_c= np.power(10.,LF(L_bol_grid,parameters)).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
res = convolve_c(input_c,nu_c,redshift_c,dtg_c)
res = [i for i in res.contents]
PHI_1450 = np.array(res,dtype=np.float64)
x = -2.5*( L_1450 + L_solar - np.log10(Fab*(con.c.value/1450e-10)) )
y = np.log10(PHI_1450) - np.log10(2.5)
ax.plot(x,y,'-',lw=4,c='crimson',label=r'$\rm Modelled$ $\rm quasar$ $\rm UVLF$')

data = np.genfromtxt("../GUVLF_data/schechter_fits.dat",names=True)
id_z = data['z']==redshift
def single_sch(M,phi_s,M_s,alpha):
        return np.log10( 0.4*np.log(10.)*phi_s*np.power(10.,-0.4*(M-M_s)*(alpha+1.))*np.exp(-np.power(10.,-0.4*(M-M_s))) )
x = np.linspace(-13,-25,100)
y = single_sch(x, 10**data['LogPhis'][id_z], data["Ms"][id_z], data["alpha"][id_z])
ax.plot(x,y,'-',lw=4,c='royalblue',label=r'$\rm Galaxy$ $\rm UVLF$ $\rm fit$')

x,y,dy,yfit=get_data()
ax.errorbar(x,y,yerr=dy,capsize=6,linestyle='',lw=2,c='gray',mec='gray',marker='o', ms=10, capthick=2 ,label=r'$\rm Quasar$ $\rm UVLF$ $\rm compilation$')
'''
if redshift==6:
	x = np.linspace(-14,-32,100) 
	y = return_giallongo15_lf_fitted(x,6.)
	ax.plot(x,y,'--',dashes=(25,15),c='seagreen',label=r'$\rm Giallongo+$ $\rm 2015$ $z=5-6.5$')
'''
if redshift==2:
	snapnum=33
	plt_obdata('obdata'+str(snapnum).zfill(3)+'_ala.dat',r'${\rm Alavi+}$ ${\rm 2014}$','darkorchid')
	plt_obdata('obdata'+str(snapnum).zfill(3)+'_par.dat',r'${\rm Parsa+}$ ${\rm 2016}$','darkorchid')
	plt_obdata('obdata'+str(snapnum).zfill(3)+'_oes.dat',r'${\rm Oesch+}$ ${\rm 2010}$','darkorchid')
	plt_obdata('obdata'+str(snapnum).zfill(3)+'_hat.dat',r'${\rm Hathi+}$ ${\rm 2010}$','darkorchid')
	plt_obdata('obdata'+str(snapnum).zfill(3)+'_met.dat',r'${\rm Mehta+}$ ${\rm 2017}$','darkorchid')
elif redshift==4:
	snapnum=21
	plt_obdata('obdata'+str(snapnum).zfill(3)+'.dat',r'${\rm Finkelstein+}$ ${\rm 2016}$'+'\n'+r'${\rm Compilation}$','darkorchid')
	plt_obdata('obdata'+str(snapnum).zfill(3)+'_par.dat','','darkorchid',label=False)
	plt_obdata('obdata'+str(snapnum).zfill(3)+'_ono.dat',r'${\rm Ono+}$ ${\rm 2018}$','darkorchid')
elif redshift==6:
	snapnum=13
	plt_obdata('obdata'+str(snapnum).zfill(3)+'.dat','','darkorchid',label=False)
	plt_obdata('obdata'+str(snapnum).zfill(3)+'_bou.dat',r'${\rm Bouwens+}$ ${\rm 2017}$','darkorchid')
	plt_obdata('obdata'+str(snapnum).zfill(3)+'_ate.dat',r'${\rm Atek+}$ ${\rm 2018}$','darkorchid')
	plt_obdata('obdata'+str(snapnum).zfill(3)+'_ono.dat','','darkorchid',label=False)

ax.errorbar([],[],yerr=([],[]),c='darkorchid',mec='darkorchid',linestyle='none',marker='o',lw=2,markersize=10,capsize=6,capthick=2,label=r'$\rm Galaxy$ $\rm UVLF$ $\rm compilation$')

prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
#ax.set_xlabel(r'$\log{(L_{\rm B}/{\rm L}_{\odot})}$',fontsize=40,labelpad=2.5)
ax.set_xlabel(r'$M_{\rm 1450}$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\phi[{\rm mag}^{-1}{\rm Mpc}^{-3}])}$',fontsize=40,labelpad=5)
ax.text(0.88, 0.92, r'${\rm z\sim'+str(redshift)+'}$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=40)

ax.set_xlim(-15.5,-30.5)
ax.set_ylim(-12.2,-0.3)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
#plt.savefig("../figs/compare_"+str(redshift)+".pdf",fmt='pdf')
plt.show()

