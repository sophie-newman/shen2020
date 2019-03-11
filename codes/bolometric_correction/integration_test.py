from data import *
import numpy as np 
import astropy.constants as con
import sedpy
from scipy.integrate import romberg
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import time

def returnIR_to_UV_H07(sed2500):
	data1=np.genfromtxt(datapath+"OpticalSED.dat",names=["lognu","logall"],)
	lamb1=con.c.value/(10**data1['lognu'])*1e10
	sed1 =data1['logall']

	lamb1, sed1 = lamb1[(lamb1>1300) & (lamb1<10000)], sed1[(lamb1>1300) & (lamb1<10000)]
	overlay = interp1d(lamb1,sed1)

	################
	data2=np.genfromtxt(datapath+"R06_SED.dat",names=["lognu","logall","sigall","blue"],)
	lamb2=con.c.value/(10**data2['lognu'])*1e10
	sed2 =data2['blue']

	################
	UV_end  =1300.
	UVslope =1.76
	
	lamb2, sed2 = lamb2[lamb2>UV_end], sed2[lamb2>UV_end]
	
	UVlamb2= np.linspace(UV_end,500.,100)
	UVsed2 = sed2[-1] + (UVslope-1) * (np.log10(UVlamb2) - np.log10(lamb2[-1]))
	
	lamb2= np.append( lamb2, UVlamb2)
	sed2 = np.append( sed2 , UVsed2+ (sed2[-1]-UVsed2[0]))
	
	####
	#sed2[(lamb2>1300) & (lamb2<lamb1[0])] += overlay(lamb2[(lamb2>1300) & (lamb2<lamb1[0])]) - 1.

	id2500 = np.arange(0,len(lamb2),dtype=np.int32)[(lamb2-2500)<0][0]
	totscale = sed2500-sed2[id2500]

	return lamb2[lamb2>500], sed2[lamb2>500]+totscale


def returnIR_to_UV(sed2500):
	data1=np.genfromtxt(datapath+"K13_SED.dat",names=["lognu","logall"],)
	freq1=data1['lognu']
	lamb1=con.c.value/(10**data1['lognu'])*1e10
	sed1 =data1['logall']
	
	#IR construction
	IR_end=lamb1[0]
	data2=np.genfromtxt(datapath+"R06_SED.dat",names=["lognu","logall","sigall","blue"],)
	lamb2=con.c.value/(10**data2['lognu'])*1e10
	sed2 =data2['blue']
	
	lamb2_ir= lamb2[lamb2>IR_end]
	sed2_ir = sed2[lamb2>IR_end]
	
	scale = ((sed1[1]-sed1[0])/(lamb1[1]-lamb1[0]) * (lamb2_ir[-1]-lamb1[0]) + sed1[0])/sed2_ir[-1]
	
	lamb1= np.append( lamb2_ir, lamb1)
	sed1 = np.append( sed2_ir*scale, sed1)
	
	#UV construction
	UV_end  =912.
	UVslope =1.70
	
	lamb1,sed1 = lamb1[lamb1>UV_end],sed1[lamb1>UV_end]
	
	lambUV= np.linspace(UV_end,500.,10)
	sedUV = sed1[-1] + (UVslope-1) * (np.log10(lambUV) - np.log10(lamb1[-1]))
	
	lamb1= np.append( lamb1, lambUV)
	sed1 = np.append( sed1 , sedUV)
	freq1= np.log10(con.c.value/(lamb1*1e-10))

	totscale = sed2500-interp1d(lamb1,sed1)(2500.)

	return lamb1, sed1+totscale

def returnXray(sed2500):
	f2500 = np.power(10.,sed2500)/(con.c.value/2500./1e-10)
	#alphaox = -0.137*np.log10(f2500)+2.638
	alphaox = -0.107*np.log10(f2500)+1.739
	#alphaox = -0.154*np.log10(f2500)+3.176
	f2kev=10**(alphaox/0.384)*f2500
	data3=np.genfromtxt(datapath+"XRAY_SED.dat",names=["lognu","logall"],)
	lamb3=con.c.value/(10**data3['lognu'])*1e10
	
	lamb2kev = con.c.value/ (2.*1000.*con.e.value/con.h.value) *1e10
	#id2kev = np.arange(0,len(lamb3),dtype=np.int32)[(lamb3-lamb2kev)<0][0]

	scale = np.log10( f2kev*(2.*1000.*con.e.value/con.h.value) ) - interp1d(lamb3,data3['logall'])(lamb2kev)
	sed3 = data3['logall'] + scale
	return lamb3[lamb3<50], sed3[lamb3<50]

def returnall(sed2500):
	lambNX, sedNX = returnIR_to_UV(sed2500)
	lambX, sedX = returnXray(sed2500)

	lamball = np.append(lambNX, lambX)
	sedall  = np.append(sedNX, sedX)
	freqall = con.c.value/(lamball*1e-10)
	#plt.plot(lamball,sedall)
	#plt.xscale('log')
	#plt.show()
	#exit()
	
	funcint1 = interp1d(freqall,sedall)
	funcint2 = interp1d(np.log10(freqall),sedall)
	def func1(x):
		return np.power(10.,funcint1(x))/x
	def func2(logx):
		x=10**logx
		return np.power(10.,funcint1(x))
	def func3(logx):
		return np.power(10.,funcint2(logx))

	t0 = time.time()
	Lbol_trial=quad(func1, con.c.value/94.8e-6, 500.*1000.*con.e.value/con.h.value)[0]
	print time.time()-t0, np.log10(Lbol_trial)
	t0 = time.time()

	Lbol_1= romberg(func1, con.c.value/94.8e-6, 500.*1000.*con.e.value/con.h.value,
		divmax=50, tol=1e-2*Lbol_trial , rtol=1e-2)
	print time.time()-t0, np.log10(Lbol_1)
	t0 = time.time()
	Lbol_2= romberg(func2, np.log10(con.c.value/94.8e-6), np.log10(500.*1000.*con.e.value/con.h.value),
		divmax=50, tol=1e-2*Lbol_trial , rtol=1e-2)
	print time.time()-t0, np.log10(Lbol_2)
	t0 = time.time()
	Lbol_3= romberg(func3, np.log10(con.c.value/94.8e-6), np.log10(500.*1000.*con.e.value/con.h.value),
		divmax=50, tol=1e-2*Lbol_trial , rtol=1e-2)
	print time.time()-t0, np.log10(Lbol_3)
	t0 = time.time()

	id = (freqall>con.c.value/94.8e-6) & (freqall<500.*1000.*con.e.value/con.h.value)
	Lbol_4 = np.trapz(10**sedall[id]/freqall[id],freqall[id])
	print time.time()-t0, np.log10(Lbol_4)
	t0 = time.time()
	Lbol_5 = np.trapz(10**sedall[id], np.log10(freqall[id]))
	print time.time()-t0, np.log10(Lbol_5)
	t0 = time.time()

returnall(10)
