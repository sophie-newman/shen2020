from data import *
import numpy as np 
import astropy.constants as con
import sedpy
from scipy.integrate import romberg
from scipy.integrate import quad
from scipy.interpolate import interp1d

def returnIR_to_UV_H07(sed2500):
	data1=np.genfromtxt(datapath+"R06_SED.dat",names=["lognu","logall","sigall","blue"],)
	lamb1=con.c.value/(10**data1['lognu'])*1e10
	sed1 =data1['blue']

	
	UV_end  =1200.
	UVslope =1.76
	lamb1,sed1 = lamb1[lamb1>UV_end],sed1[lamb1>UV_end]
	
	lambUV= np.linspace(UV_end,500.,10)
	sedUV = sed1[-1] + (UVslope-1) * (np.log10(lambUV) - np.log10(lamb1[-1]))
	
	lamb1= np.append( lamb1, lambUV)
	sed1 = np.append( sed1 , sedUV)
	
	freq1= np.log10(con.c.value/(lamb1*1e-10))


	id2500 = np.arange(0,len(lamb1),dtype=np.int32)[(lamb1-2500)<0][0]
	totscale = sed2500-sed1[id2500]

	return lamb1[lamb1>500], freq1[lamb1>500], sed1[lamb1>500]+totscale


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

	id2500 = np.arange(0,len(lamb1),dtype=np.int32)[(lamb1-2500)<0][0]
	totscale = sed2500-sed1[id2500]

	return lamb1, freq1, sed1+totscale

def returnXray(sed2500):
	f2500 = np.power(10.,sed2500)/(con.c.value/2500./1e-10)
	#alphaox = -0.137*np.log10(f2500)+2.638
	alphaox = -0.107*np.log10(f2500)+1.739
	f2kev=10**(alphaox/0.384)*f2500
	data3=np.genfromtxt(datapath+"XRAY_SED.dat",names=["lognu","logall"],)
	lamb3=con.c.value/(10**data3['lognu'])*1e10
	#ratio= data3['logall'][np.abs(10**data3['lognu']*con.h.value/con.e.value/1000-2)<=0.05]-(-1.47406)
	#L_HX = np.log10( f2kev*(2.*1000.*con.e.value/con.h.value) / 10**ratio )
	#sed3 = data3['logall']-(-1.47406)+L_HX
	scale = np.log10( f2kev*(2.*1000.*con.e.value/con.h.value)) - data3['logall'][np.abs(10**data3['lognu']*con.h.value/con.e.value/1000-2)<=0.05]
	sed3 = data3['logall'] + scale
	return lamb3[lamb3<50], sed3[lamb3<50]

def returnall(sed2500):
	lambNX, freqNX, sedNX = returnIR_to_UV_H07(sed2500)
	lambX, sedX = returnXray(sed2500)

	lamball = np.append(lambNX, lambX)
	sedall  = np.append(sedNX, sedX)
	freqall = con.c.value/(lamball*1e-10)

	func = interp1d(freqall,(10**sedall)/freqall)
	LHX =quad(func, 2.*1000.*con.e.value/con.h.value,10.*1000.*con.e.value/con.h.value)[0]
	LSX =quad(func, 0.5*1000.*con.e.value/con.h.value,2.*1000.*con.e.value/con.h.value)[0]
	logLB  =sedall[(lamball-4450.)<0][0]
	logLIR =sedall[(lamball-150000.)<0][0]
	#Lbol=romberg(func, con.c.value/30e-6, 10.*1000.*con.e.value/con.h.value,divmax=20, tol=1e-2*10**(sed2500) , rtol=1e-2)
	Lbol=quad(func, con.c.value/100e-6, 500.*1000.*con.e.value/con.h.value)[0]
	return  np.log10(Lbol), np.log10(LHX), np.log10(LSX), logLB, logLIR

sed2500s = np.linspace(5,15,100)+L_solar
Lbols = 0*sed2500s
LHXs = 0*sed2500s
LSXs = 0*sed2500s
LBs  = 0*sed2500s
LIRs = 0*sed2500s
for i in range(len(sed2500s)):
	Lbols[i],LHXs[i],LSXs[i],LBs[i],LIRs[i] = returnall(sed2500s[i])

np.savetxt("bolcorr.dat", np.c_[Lbols,LHXs,LSXs,LBs,LIRs], header='Lbols,LHXs,LSXs,LBs,LIRs')




