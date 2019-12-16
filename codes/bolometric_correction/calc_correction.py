from data import *
import numpy as np 
import astropy.constants as con
import sedpy
from scipy.integrate import romberg
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
def returnIR_to_UV_H07(sed2500):
	''' # overlay optical SEDs with emission lines
	data1=np.genfromtxt(datapath+"OpticalSED.dat",names=["lognu","logall"],)
	lamb1=con.c.value/(10**data1['lognu'])*1e10
	sed1 =data1['logall']

	lamb1, sed1 = lamb1[(lamb1>1300) & (lamb1<10000)], sed1[(lamb1>1300) & (lamb1<10000)]
	overlay = interp1d(lamb1,sed1)
	'''
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
	freq2= con.c.value/(lamb2*1e-10)

	####
	#sed2[(lamb2>1300) & (lamb2<lamb1[0])] += overlay(lamb2[(lamb2>1300) & (lamb2<lamb1[0])]) - 1.

	totscale = sed2500-interp1d(np.log10(freq2),sed2)(np.log10(con.c.value/2500e-10))

	return lamb2[lamb2>500], sed2[lamb2>500]+totscale

def returnXray_H07(sed2500):
	f2500 = np.power(10.,sed2500)/(con.c.value/2500./1e-10)

	alphaox = -0.107*np.log10(f2500)+1.739

	f2kev=10**(alphaox/0.384)*f2500
	data3=np.genfromtxt(datapath+"XRAY_SED.dat",names=["lognu","logall"])	
	freq3=10**data3["lognu"]
	sed3 = data3["logall"]
	lamb3= con.c.value/freq3*1e10

	lamb2kev = con.c.value/ (2.*1000.*con.e.value/con.h.value) *1e10
	freq2kev = 2.*1000.*con.e.value/con.h.value

	scale = np.log10( f2kev*(2.*1000.*con.e.value/con.h.value) ) - interp1d(np.log10(freq3),sed3)(np.log10(freq2kev))
	sed3 = sed3 + scale
	return lamb3[lamb3<50], sed3[lamb3<50]

def returnIR_to_UV(sed2500):  #return the UV to IR SED
	data1=np.genfromtxt(datapath+"K13_SED.dat",names=["lognu","logall"],)
	freq1=10**data1['lognu']
	lamb1=con.c.value/(10**data1['lognu'])*1e10
	sed1 =data1['logall']
	
	#IR construction, using the R06 SED at the wavelengths beyond the coverage of K13
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
	sedUV = sed1[-1] + (UVslope-1) * (np.log10(lambUV) - np.log10(lamb1[-1])) #extrapolate to EUV with a power-law
	
	lamb1= np.append( lamb1, lambUV)
	sed1 = np.append( sed1 , sedUV)
	freq1= con.c.value/(lamb1*1e-10)

	totscale = sed2500-interp1d(np.log10(freq1),sed1)(np.log10(con.c.value/2500e-10))

	return lamb1, sed1+totscale

def returnXray(sed2500): #return the X-ray SED
	f2500 = np.power(10.,sed2500)/(con.c.value/2500./1e-10)

	#beta, C = 0.721, 4.531	
	#A, Cprime = 0.384*(1-beta), 0.384*C
	#alphaox = -A*np.log10(f2500)+Cprime
	alphaox = -0.107*np.log10(f2500)+1.739 #fiducial model

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
	'''
	lambNX, sedNX = returnIR_to_UV_H07(sed2500)
	lambX, sedX = returnXray(sed2500)

	lamball = np.append(lambNX, lambX)
	sedall  = np.append(sedNX, sedX)
	freqall = con.c.value/(lamball*1e-10)	
	#plt.loglog(lamball, sedall,'r-')
	'''
	lambNX, sedNX = returnIR_to_UV(sed2500)
	lambX, sedX = returnXray(sed2500)

	lamball = np.append(lambNX, lambX)
	sedall  = np.append(sedNX, sedX)
	freqall = con.c.value/(lamball*1e-10)	
	
	#plt.loglog(lamball, sedall,'b--')
	#plt.show()
	#exit()

	def integrate(sed, freq, fmin, fmax):
		idmin= np.arange(0,len(freq),dtype=np.int32)[freq>fmin][0]
		idmax= np.arange(0,len(freq),dtype=np.int32)[freq<fmax][-1]

		sed_at_fmin= sed[idmin-1] + (np.log10(fmin)-np.log10(freq[idmin-1])) * (sed[idmin]-sed[idmin-1])/(np.log10(freq[idmin])-np.log10(freq[idmin-1]))
		sed_at_fmax= sed[idmax]   + (np.log10(fmax)-np.log10(freq[idmax]))   * (sed[idmax+1]-sed[idmax])/(np.log10(freq[idmax+1])-np.log10(freq[idmax]))

		freq = np.append(np.append(fmin,freq[idmin:idmax+1]),fmax)
		sed  = np.append(np.append(sed_at_fmin,sed[idmin:idmax+1]),sed_at_fmax)

		return np.trapz(np.log(10)*10**sed, np.log10(freq))

	def tophat(sed, freq, fmin, fmax):
		x_target = np.linspace( np.log10(fmin), np.log10(fmax), 100)
		y_target = interp1d(np.log10(freq), sed)(x_target)
		return np.mean(y_target)

	LHX =integrate( sedall, freqall, 2.*1000.*con.e.value/con.h.value, 10.*1000.*con.e.value/con.h.value)

	LSX =integrate( sedall, freqall, 0.5*1000.*con.e.value/con.h.value,2.*1000.*con.e.value/con.h.value)
	
	logLB  = tophat( sedall, freqall, con.c.value/((4450+470)*1e-10), con.c.value/((4450-470)*1e-10))

	logL1450 = tophat( sedall, freqall, con.c.value/(1500e-10), con.c.value/(1400e-10))

	logLIR = tophat( sedall, freqall, con.c.value/((15+1.)*1e-6), con.c.value/((15-1.)*1e-6))

	Lbol= integrate( sedall, freqall, con.c.value/(30.*1e-6), 500*1000.*con.e.value/con.h.value)

	#check the wavelength range of bolometric luminosity
	Lbol2= integrate( sedall, freqall, con.c.value/(2.*1e-6), 500*1000.*con.e.value/con.h.value)
	print np.log10(Lbol) - np.log10(Lbol2)

	#check the correction factor in IR
	#logLNIR = tophat( sedall, freqall, con.c.value/((1.25+0.15)*1e-6), con.c.value/((1.25-0.15)*1e-6))
	#print logLIR - logLNIR
	#MB    = -2.5*( logLB - np.log10(Fab*con.c.value/4450e-10) ) 
	#M1450 = -2.5*( logL1450 - np.log10(Fab*con.c.value/1450e-10) )
	#print (M1450-MB)
	#exit()

	#check the correction facotr of Ultra hard xray band
	#LUHX = integrate( sedall, freqall, 10.*1000.*con.e.value/con.h.value, 40.*1000.*con.e.value/con.h.value)
	#print np.log10(LUHX) - np.log10(LHX)
	#exit()

	return  np.log10(Lbol), np.log10(LHX), np.log10(LSX), logLB, logL1450, logLIR

#sed2500s = np.linspace(5,15,200)+L_solar
sed2500s = np.linspace(10,14,10)+L_solar
Lbols = 0*sed2500s
LHXs  = 0*sed2500s
LSXs  = 0*sed2500s
LBs   = 0*sed2500s
L1450s= 0*sed2500s
LIRs  = 0*sed2500s
for i in range(len(sed2500s)):
	Lbols[i],LHXs[i],LSXs[i],LBs[i],L1450s[i],LIRs[i] = returnall(sed2500s[i])

exit()	
np.savetxt("bolcorr.dat", np.c_[Lbols,LHXs,LSXs,LBs,L1450s,LIRs], header='Lbols,LHXs,LSXs,LBs,L1450s,LIRs')
