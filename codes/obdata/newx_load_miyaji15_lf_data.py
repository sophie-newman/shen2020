# Miyaji et al. 2015
# measured 2-10 keV
from data import *
import numpy as np

def load_miyaji15_lf_data(z): 
	if (z <= 0.015) or (z > 4.60): return False
	else:	
		filename = datapath+'miyaji2015_table4.dat'
		data = np.genfromtxt(filename,names=['z_min', 'z_max', 'z_mean', 
			'L_min', 'L_max', 'L_c', 'nqso', 'phi', 'phierr'])

		id = (data["z_min"]< z) & (data["z_max"]>=z)

		L_X  = (data['L_min'][id]+data['L_max'][id])/2.
		phi  = data['phi'][id]
		phierr  = data['phierr'][id]
			 
		PHI_X = np.log10(phi)
		DPHI_X = ((np.log10(phi+phierr)-np.log10(phi)) + (np.log10(phi)-np.log10(phi-phierr)))/2.
		return L_X, PHI_X, DPHI_X

def return_miyaji15_lf_fitted(Llist,z):
	p1_44,p2_44,p3 = 5.29, -0.35, -5.6
	zb1_44,zb2 = 1.1, 2.7
	alpha = 0.18
	beta1, beta2 = 1.2, 1.5

	gamma1 = 1.17
	gamma2 = 2.80
	lbreak = 44.04
	phi_s = 1.56e-6

	def ed(Lx,z):
		p1 = p1_44 + beta1*(Lx - 44.)
		p2 = p2_44 + beta2*(Lx - 44.)

		zb1_0 = zb1_44 / np.power( 10**44/10**44.5, alpha )
		if Lx <= 44.5: zb1 = zb1_0 * np.power( 10**Lx/10**44.5, alpha )
		else: zb1 = zb1_0

		if z < zb1:  ed = (1.+z)**p1
		elif z < zb2: ed = (1.+zb1)**p1 * np.power( (1+z)/(1+zb1) ,p2)
		else: ed = (1.+zb1)**p1 * np.power( (1+zb2)/(1+zb1) ,p2) * np.power( (1+z)/(1+zb2) ,p3)

		return ed

	x = np.power(10., Llist-lbreak)
	eds = np.zeros(len(Llist))
	for i in range(len(eds)):
		eds[i] = ed(Llist[i],z)
	return np.log10( phi_s/( np.power(x,gamma1) + np.power(x,gamma2))*eds )
