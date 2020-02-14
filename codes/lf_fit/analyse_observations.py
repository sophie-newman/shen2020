from data import *
import numpy as np
from lf_shape import *
import scipy.interpolate as inter
from convolve import *

import matplotlib.pyplot as plt
from lf_fitter_data import *

def get_data(alldata,zmin,zmax,dset_name,dset_id):
        z_lis = 0.5*(zmin + zmax)

        alldata_tem={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([])}

        for iz in range(len(z_lis)):
                redshift = z_lis[iz]
                L_data, PHI_data, DPHI_data = load_LF_data[dset_name](redshift)

                if dset_id==-5:
                        L_1450 = bolometric_correction(L_bol_grid,dset_id) + L_solar
                        M_1450 = -2.5*( L_1450 - np.log10(Fab*con.c.value/1450e-10) )
                        L_tmp  = np.sort(M_1450)
                else: L_tmp=bolometric_correction(L_bol_grid,dset_id)

                if return_LF[dset_name]!=None:
                        #if dset_id != -4:
                        phi_fit_tmp = return_LF[dset_name](L_tmp, redshift)
                        phi_fit_pts = np.interp(L_data ,L_tmp, phi_fit_tmp)
                        PHI_data = PHI_data + (np.mean((phi_fit_pts))-np.mean((PHI_data)))

                        if dset_id==-5:
                                L_1450 = (-0.4*L_data) + np.log10(Fab*(con.c.value/1450e-10)) - L_solar

                        alldata["L_OBS"]  = np.append(alldata["L_OBS"]  , L_data)
                        alldata["P_OBS"]  = np.append(alldata["P_OBS"]  , PHI_data)
                        alldata["D_OBS"]  = np.append(alldata["D_OBS"]  , DPHI_data + 0.01)
                        alldata["Z"]  = np.append(alldata["Z"]  , np.ones(len(L_data)) * redshift)
                        alldata["ID"] = np.append(alldata["ID"] , np.ones(len(L_data)) * dset_id)

                        alldata_tem["L_OBS"]  = np.append(alldata_tem["L_OBS"]  , L_data)
                        alldata_tem["P_OBS"]  = np.append(alldata_tem["P_OBS"]  , PHI_data)
                        alldata_tem["D_OBS"]  = np.append(alldata_tem["D_OBS"]  , DPHI_data + 0.01)

		if (dset_id==-5) and (zmax[iz]>=4.):
			print "-------------------------------------------------------------------"
	        	print "NAME:",dset_name
			print "z   :",zmin[iz]," - ",zmax[iz]
			print "L   :",L_data
			print "Phi :",PHI_data
			print "DPhi:",DPHI_data
			print "-------------------------------------------------------------------"


alldata={"P_PRED":np.array([]),"L_OBS":np.array([]),"P_OBS":np.array([]),"D_OBS":np.array([]),"Z":np.array([]),"WEIGHT":np.array([]),"ID":np.array([])}
for key in dset_ids.keys():
	get_data(alldata,zmins[key],zmaxs[key],key,dset_ids[key])

