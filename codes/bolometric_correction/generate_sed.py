from data import *
import numpy as np 
import astropy.constants as con
from scipy.interpolate import interp1d

def returnXray_H07(sed2500):
	f2500 = np.power(10.,sed2500)/(con.c.value/2500./1e-10)

	alphaox = -0.107*np.log10(f2500)+1.739

	f2kev=10**(alphaox/0.384)*f2500
	data3=np.genfromtxt(datapath+"XRAY_SED.dat",names=["lognu","logall"])	
	freq3=10**data3["lognu"]
	sed3 = data["logall"]
	lamb3= con.c.value/freq3*1e10

	lamb2kev = con.c.value/ (2.*1000.*con.e.value/con.h.value) *1e10
	freq2kev = 2.*1000.*con.e.value/con.h.value

	scale = np.log10( f2kev*(2.*1000.*con.e.value/con.h.value) ) - interp1d(np.log10(freq3),sed3)(np.log10(freq2kev))
	sed3 = sed3 + scale
	return lamb3[lamb3<50], sed3[lamb3<50]

def returnIR_to_UV(sed2500):
	data1=np.genfromtxt(datapath+"K13_SED.dat",names=["lognu","logall"],)
	freq1=10**data1['lognu']
	lamb1=con.c.value/(10**data1['lognu'])*1e10
	sed1 =data1['logall']
	
	#IR construction
	IR_end=lamb1[0]
	data2=np.genfromtxt(datapath+"R06_SED.dat",names=["lognu","logall","sigall","blue"],)
	lamb2=con.c.value/(10**data2['lognu'])*1e10
	sed2 =data2['blue']
	
	lamb2_ir= lamb2[lamb2>IR_end]
	sed2_ir = sed2[lamb2>IR_end]
	
	scale = ( (sed1[1]-sed1[0])/(np.log10(lamb1[1])-np.log10(lamb1[0])) * (np.log10(lamb2_ir[-1])-np.log10(lamb1[0])) + sed1[0] ) - sed2_ir[-1]
	
	lamb1= np.append( lamb2_ir, lamb1)
	sed1 = np.append( sed2_ir + scale, sed1)
	
	#UV construction
	UV_end  =912.
	UVslope =1.70
	
	lamb1,sed1 = lamb1[lamb1>UV_end],sed1[lamb1>UV_end]
	
	lambUV= np.linspace(UV_end, 600.,100)
	sedUV = sed1[-1] + (UVslope-1) * (np.log10(lambUV) - np.log10(lamb1[-1]))
	
	lamb1= np.append( lamb1, lambUV)
	sed1 = np.append( sed1 , sedUV)
	freq1= con.c.value/(lamb1*1e-10)

	totscale = sed2500-interp1d(np.log10(freq1),sed1)(np.log10(con.c.value/2500e-10))

	return lamb1, sed1+totscale

def returnXray(sed2500):
	f2500 = np.power(10.,sed2500)/(con.c.value/2500./1e-10)

	alphaox = -0.107*np.log10(f2500)+1.739

	f2kev=10**(alphaox/0.384)*f2500
	data3=np.genfromtxt("xspec_lib/Xspec_1_9.dat",names=["E","f"],)	
	freq3=data3['E']*1000.*con.e.value/con.h.value
	lamb3=con.c.value/freq3 * 1e10
	sed3 =np.log10(data3['f']*data3['E'])

	lamb2kev = con.c.value/ (2.*1000.*con.e.value/con.h.value) *1e10
	freq2kev = 2.*1000.*con.e.value/con.h.value

	scale = np.log10( f2kev*(2.*1000.*con.e.value/con.h.value) ) - interp1d(np.log10(freq3),sed3)(np.log10(freq2kev))
	sed3 = sed3 + scale
	return lamb3[lamb3<50], sed3[lamb3<50]

def returnall(sed2500):
	lambNX, sedNX = returnIR_to_UV(sed2500)
	lambX, sedX = returnXray(sed2500)

	lamball = np.append(lambNX, lambX)
	sedall  = np.append(sedNX, sedX)
	freqall = con.c.value/(lamball*1e-10)	
	return lamball, sedall

lamb1, sed1 = returnall(45.55)
#import matplotlib.pyplot as plt 
#plt.loglog(lamb1,sed1)
#plt.show()
np.savetxt(datapath+"MySED.dat",np.c_[lamb1,sed1])
