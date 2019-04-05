from data import *
import numpy as np 
import astropy.constants as con
import sedpy
from scipy.integrate import romberg
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def returnIR_to_UV(sed2500, slope):
	optical_slope = slope
	lamb1= np.linspace(1e4,700.,100)
	sed1 = 1. + (optical_slope-1) * (np.log10(lamb1) - np.log10(lamb1[-1]))
	
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
	
	lambUV= np.linspace(UV_end,600.,100)
	sedUV = sed1[-1] + (UVslope-1) * (np.log10(lambUV) - np.log10(lamb1[-1]))
	
	lamb1= np.append( lamb1, lambUV)
	sed1 = np.append( sed1 , sedUV)
	freq1= con.c.value/(lamb1*1e-10)

	totscale = sed2500-interp1d(np.log10(freq1),sed1)(np.log10(con.c.value/2500e-10))

	return lamb1, sed1+totscale

def returnXray(sed2500, gamma, alphaox):
	f2500 = np.power(10.,sed2500)/(con.c.value/2500./1e-10)

	f2kev=10**(alphaox/0.384)*f2500
	#data3=np.genfromtxt("xspec_lib/Xspec_"+str(gamma).replace('.','_')+".dat",names=["lognu","logall"])
	#lamb3=con.c.value/(10**data3['lognu'])*1e10
	#freq3=10**data3['lognu']
	data3=np.genfromtxt("xspec_lib/Xspec_"+str(gamma).replace('.','_')+".dat",names=["E","f"])
	sed3 =np.log10(data3['E']*data3['f'])
	freq3 = data3['E']*1000.*con.e.value/con.h.value
	lamb3 = con.c.value/freq3*1e10

	lamb2kev = con.c.value/ (2.*1000.*con.e.value/con.h.value) *1e10
	freq2kev = 2.*1000.*con.e.value/con.h.value

	scale = np.log10( f2kev*(2.*1000.*con.e.value/con.h.value) ) - interp1d(np.log10(freq3),sed3)(np.log10(freq2kev))
	sed3 = sed3 + scale
	return lamb3[lamb3<50], sed3[lamb3<50]

def returnall(sed2500, slope, gamma, alphaox):
	lambNX, sedNX = returnIR_to_UV(sed2500,slope)
	lambX, sedX = returnXray(sed2500,gamma,alphaox)

	lamball = np.append(lambNX, lambX)
	sedall  = np.append(sedNX, sedX)
	freqall = con.c.value/(lamball*1e-10)	
	
	plt.loglog(lamball,sedall)
	plt.show()
	exit()

	def integrate(sed, freq, fmin, fmax):
		idmin= np.arange(0,len(freq),dtype=np.int32)[freq>fmin][0]
		idmax= np.arange(0,len(freq),dtype=np.int32)[freq<fmax][-1]

		sed_at_fmin= sed[idmin-1] + (fmin-freq[idmin-1]) * (sed[idmin]-sed[idmin-1])/(freq[idmin]-freq[idmin-1])
		sed_at_fmax= sed[idmax]   + (fmax-freq[idmax])   * (sed[idmax+1]-sed[idmax])/(freq[idmax+1]-freq[idmax])

		freq = np.append(np.append(fmin,freq[idmin:idmax+1]),fmax)
		sed  = np.append(np.append(sed_at_fmin,sed[idmin:idmax+1]),sed_at_fmax)

		return np.trapz(np.log(10)*10**sed, np.log10(freq))

	LHX =integrate( sedall, freqall, 2.*1000.*con.e.value/con.h.value, 10.*1000.*con.e.value/con.h.value)

	LSX =integrate( sedall, freqall, 0.5*1000.*con.e.value/con.h.value,2.*1000.*con.e.value/con.h.value)
	
	logLB  =np.mean(sedall[np.abs(lamball-4450.)<940])
	logLIR =np.mean(sedall[np.abs(lamball-150000.)<75000])

	Lbol= integrate( sedall, freqall, con.c.value/(30.*1e-6), 500*1000.*con.e.value/con.h.value)
	return  np.log10(Lbol), np.log10(LHX), np.log10(LSX), logLB, logLIR

def weighted_std(values, weights):
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)
    return np.sqrt(variance)

slope_uncertainty, slope_best = 0.125, 0.44
pindex_uncertainty, pindex_best = 0.3, 1.9
alphaox_uncertainty, alphaox_best = 1., 0.107

pindexlist = np.genfromtxt("./xspec_lib/pindex_list.dat")
pindex_weight = 1./np.sqrt(2.*np.pi*pindex_uncertainty)*np.exp(-((pindexlist-pindex_best)/pindex_uncertainty)**2)
slopelist = 0.44 + np.linspace(-0.25,0.25,1)
slope_weight = 1./np.sqrt(2.*np.pi*slope_uncertainty)*np.exp(-((slopelist-slope_best)/slope_uncertainty)**2)
alphaoxlist = 0.107 + np.linspace(-1.2,1.2,1)
alphaox_weight = 1./np.sqrt(2.*np.pi*alphaox_uncertainty)*np.exp(-((alphaoxlist-alphaox_best)/alphaox_uncertainty)**2)

sed2500s = np.linspace(5,15,3)+L_solar
sigmaHX = 0*sed2500s
sigmaSX = 0*sed2500s
sigmaB  = 0*sed2500s
sigmaIR = 0*sed2500s
for i in range(len(sed2500s)):
	HXcorrs =np.array([])
	SXcorrs =np.array([])
	Bcorrs  =np.array([])
	IRcorrs =np.array([])
	weights =np.array([])
	for j in range(len(slopelist)):
		for k in range(len(alphaoxlist)):
			for m in range(len(pindexlist)): 
				print i,j,k,m
				Lbol,LHX,LSX,LB,LIR = returnall(sed2500s[i],slopelist[j],pindexlist[m],alphaoxlist[k])
				print Lbol, LHX, Lbol-LHX
				HXcorrs=np.append(HXcorrs, Lbol-LHX)
				weights=np.append(weights, slope_weight[j]*alphaox_weight[k]*pindex_weight[m])
	print HXcorrs
	sigmaHX[i]=weighted_std(HXcorrs, weights)

print sigmaHX
#np.savetxt("sigmas.dat", np.c_[sigmaHX,sigmaSX,sigmaB,sigmaIR], header='HX,SX,B,IR')


