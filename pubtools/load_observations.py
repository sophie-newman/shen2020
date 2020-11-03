from data_copy import *
import numpy as np
from lf_shape_copy import *
import scipy.interpolate as inter
from utilities import *
from lf_fitter_data import *
from ctypes import *
import ctypes
import sys

# load observational data at a certain redshift and move them to the bolometric plane

redshift=float(sys.argv[1])
############################
def get_fit_data(alldata,dset_name,dset_id):
        # load the observational data of a specific dataset at the redshift
        # alldata: a container for observational data
        # dset_name: name of the dataset
        # dset_id: band id of the dataset, see utilities.py for definitions
        # the loaded observational data will be put into alldata

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
                phi_fit_tmp1 = return_LF[dset_name](L_tmp, redshift)
                phi_fit_pts1 = np.interp(L_data ,L_tmp, phi_fit_tmp1)
                PHI_data = PHI_data + (np.mean((phi_fit_pts1))-np.mean((PHI_data)))

        if (len(L_data) > 0):
                # load the best-fit bolometric QLF model and get the LF in the band
                L_model, PHI_model = return_qlf_in_band(redshift, dset_id, model='A', output_mag_for_UV=False)
                L_model = L_model - L_solar
                PHI_model = 10**PHI_model

                if dset_id==-5:
                        L_1450 = (-0.4*L_data) + np.log10(Fab*(con.c.value/1450e-10)) - L_solar
                        phi_i = np.interp(L_1450, L_model, np.log10(PHI_model))
                        phi_i = phi_i - np.log10(2.5)
                else:
                        phi_i = np.interp(L_data, L_model, np.log10(PHI_model))

                alldata["P_PRED"] = np.append(alldata["P_PRED"] , phi_i)  # number density predicted by the best-fit model
                alldata["L_OBS"]  = np.append(alldata["L_OBS"]  , L_data) # luminosity of the data point
                alldata["P_OBS"]  = np.append(alldata["P_OBS"]  , PHI_data) # number density from the data point
                alldata["D_OBS"]  = np.append(alldata["D_OBS"]  , DPHI_data + 0.01) # uncertainty in number density
                alldata["Z_TOT"]  = np.append(alldata["Z_TOT"]  , np.ones(len(L_data)) * redshift) # redshift
                alldata["ID"]     = np.append(alldata["ID"]     , np.ones(len(L_data)) * dset_id) # band id

                #determine the length of the luminosity bins
                #very rough estimate here, based on the distance to the next luminosity
                #only work for the bright end
                if dset_id==-5: d_Ldata = (L_data[1:] - L_data[:-1]) * 0.4
                else: d_Ldata = L_data[1:] - L_data[:-1]
                d_Ldata = np.append(np.array([0]), d_Ldata)/2.
                alldata["D_LOBS"] = np.append(alldata["D_LOBS"], d_Ldata)

def get_data(dataid):
        # main function to get the observational data in a certain band, and move them onto the bolometric plane
        alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z_TOT":np.array([]),"B":np.array([]),"ID":np.array([]),"D_LOBS":np.array([])}
        for key in dset_ids.keys():
                get_fit_data(alldata,key,dset_ids[key])

        bad = np.invert(np.isfinite(alldata["P_PRED"]))
        if (np.count_nonzero(bad) > 0): alldata["P_PRED"][bad] = -40.0

        if (dataid!=-5):
            select = alldata["ID"]==dataid
            logLbol = 0.0*alldata["L_OBS"][select]
            for i,Lband in enumerate(alldata["L_OBS"][select]):
                logLbol[i] = bolometric_correction_inverse(Lband,dataid)  # get the corresponding bolometric luminosities
        elif (dataid==-5):
            select = alldata["ID"]==dataid
            L_1450 = -0.4*alldata["L_OBS"][select] + np.log10(Fab*(con.c.value/1450e-10)) - L_solar
            logLbol = 0.0*L_1450
            for i,Lband in enumerate(L_1450):
                logLbol[i] = bolometric_correction_inverse(Lband,dataid)

        x,y = return_bolometric_qlf(redshift, model='A')
        x = x - L_solar
        phi_fit_pts = np.interp(logLbol, x, y)

        # the relative (relative to the best-fit model) postion on the bolometric plane is kept the same as the relative postion on the band LF plane
        return logLbol+L_solar, phi_fit_pts+(alldata["P_OBS"][select]-alldata["P_PRED"][select]),alldata["D_OBS"][select], alldata["D_LOBS"][select]


# below is a plot that shows the best-fit bolometric QLF along with the data points moved onto the bolometric plane
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

x,y = return_bolometric_qlf(redshift, model='A')
ax.plot(x,y,'--',dashes=(25,15),c='darkorchid',label=r'$\rm Global$ $\rm fit$ $\rm A$')

x,y,dy,dx=get_data(dataid=-5)
id1 = (x>48.5)
ax.errorbar(x,y,yerr=dy,fillstyle='none',linestyle='none',c='seagreen',mec='seagreen',marker='o',ms=10,capsize=6,capthick=2,lw=2,mew=2,label=r'$\rm UV$ $\rm 1450\AA$')
ax.errorbar(x[id1],y[id1],yerr=dy[id1],xerr=dx[id1],fillstyle='none',linestyle='none',c='seagreen',mec='seagreen',marker='o',ms=10,capsize=6,capthick=2,lw=2,mew=2)

x,y,dy,dx=get_data(dataid=-1)
id1 = (x>48.5)
ax.errorbar(x,y,yerr=dy,fillstyle='none',linestyle='none',c='pink',mec='pink',marker='o',ms=10,capsize=6,capthick=2,lw=2,mew=2,label=r'$\rm B$ $\rm Band$')
ax.errorbar(x[id1],y[id1],yerr=dy[id1],xerr=dx[id1],fillstyle='none',linestyle='none',c='pink',mec='pink',marker='o',ms=10,capsize=6,capthick=2,lw=2,mew=2)

x,y,dy,dx=get_data(dataid=-4)
id1 = (x>48.5) #& (dx>0.5)
ax.errorbar(x,y,yerr=dy,fillstyle='none',linestyle='none',c='royalblue',mec='royalblue',marker='o',ms=10,capsize=6,capthick=2,lw=2,mew=2,label=r'$\rm Hard$ $\rm X-ray$')
ax.errorbar(x[id1],y[id1],yerr=dy[id1],xerr=dx[id1],fillstyle='none',linestyle='none',c='royalblue',mec='royalblue',marker='o',ms=10,capsize=6,capthick=2,lw=2,mew=2)

x,y,dy,dx=get_data(dataid=-3)
id1 = (x>48.5) #& (dx>0.5)
ax.errorbar(x,y,yerr=dy,fillstyle='none',linestyle='none',c='gray',mec='gray',marker='o',ms=10,capsize=6,capthick=2,lw=2,mew=2,label=r'$\rm Soft$ $\rm X-ray$')
ax.errorbar(x[id1],y[id1],yerr=dy[id1],xerr=dx[id1],fillstyle='none',linestyle='none',c='gray',mec='gray',marker='o',ms=10,capsize=6,capthick=2,lw=2,mew=2)

x,y,dy,dx=get_data(dataid=-2)
id1 = (x>48.5)
ax.errorbar(x,y,yerr=dy,fillstyle='none',linestyle='none',c='olive',mec='olive',marker='o',ms=10,capsize=6,capthick=2,lw=2,mew=2,label=r'$\rm Mid$ $\rm IR$')
ax.errorbar(x[id1],y[id1],yerr=dy[id1],xerr=dx[id1],fillstyle='none',linestyle='none',c='olive',mec='olive',marker='o',ms=10,capsize=6,capthick=2,lw=2,mew=2)

prop = matplotlib.font_manager.FontProperties(size=22.0)
ax.legend(prop=prop, numpoints=1, borderaxespad=0.5, loc=3, ncol=1, frameon=False)
ax.set_xlabel(r'$\log{(L_{\rm bol}[{\rm erg}\,{\rm s}^{-1}])}$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\phi[{\rm dex}^{-1}{\rm cMpc}^{-3}])}$',fontsize=40,labelpad=5)

ax.text(0.88, 0.92, r'${\rm z\sim'+str(redshift)+'}$' ,horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=40)

ax.set_xlim(42.5,50.3)
ax.set_ylim(-11.2,-2.3)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
#plt.savefig("../figs/bol_"+str(redshift)+".pdf",fmt='pdf')
plt.show()
