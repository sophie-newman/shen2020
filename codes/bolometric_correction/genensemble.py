from data import *
import numpy as np 
import astropy.constants as con
import sedpy
from scipy.integrate import romberg
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

from scipy.stats import binned_statistic

def returnIR_to_UV(sed2500, slope):
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

	sed1 = sed1[lamb1 > 1e4]
	lamb1= lamb1[lamb1 > 1e4]

	optical_slope = slope
	lamb_op= np.linspace(1e4, 700., 100)
	sed_op = sed1[-1] + (optical_slope-1) * (np.log10(lamb_op) - np.log10(lamb1[-1]))
	
	lamb1 = np.append(lamb1, lamb_op)
	sed1  = np.append(sed1,  sed_op)
	
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

def returnXray(sed2500, gamma, alphaox_disp):
	f2500 = np.power(10.,sed2500)/(con.c.value/2500./1e-10)

	alphaox = -0.107*np.log10(f2500) + 1.739 + alphaox_disp

	f2kev=10**(alphaox/0.384)*f2500
	
	if (gamma<photon_index_list.max()) and (gamma>=photon_index_list.min()):
		pup = photon_index_list[photon_index_list>gamma][0]
		plo = photon_index_list[photon_index_list<=gamma][-1]
	if gamma>=photon_index_list.max():
		pup = photon_index_list[-1]
		plo = photon_index_list[-2]
	if gamma<photon_index_list.min():
		pup = photon_index_list[1]
		plo = photon_index_list[0]

	dataup = np.genfromtxt("xspec_lib/Xspec_"+str(pup).replace('.','_')+".dat",names=["E","f"])
	datalo = np.genfromtxt("xspec_lib/Xspec_"+str(plo).replace('.','_')+".dat",names=["E","f"])
	data3 = {
		'E': datalo['E'] + (gamma - plo)*(dataup['E']-datalo['E'])/(pup-plo),
		'f': datalo['f']
	}

	sed3 =np.log10(data3['E']*data3['f'])
	freq3 = data3['E']*1000.*con.e.value/con.h.value
	lamb3 = con.c.value/freq3*1e10

	lamb2kev = con.c.value/ (2.*1000.*con.e.value/con.h.value) *1e10
	freq2kev = 2.*1000.*con.e.value/con.h.value

	scale = np.log10( f2kev*(2.*1000.*con.e.value/con.h.value) ) - interp1d(np.log10(freq3),sed3)(np.log10(freq2kev))
	sed3 = sed3 + scale
	return lamb3[lamb3<50], sed3[lamb3<50]

def returnall(sed2500, slope, gamma, alphaox_disp):
	lambNX, sedNX = returnIR_to_UV(sed2500,slope)
	lambX, sedX = returnXray(sed2500, gamma, alphaox_disp)

	lamball = np.append(lambNX, lambX)
	sedall  = np.append(sedNX, sedX)
	freqall = con.c.value/(lamball*1e-10)	

	def integrate(sed, freq, fmin, fmax):
		idmin= np.arange(0,len(freq),dtype=np.int32)[freq>fmin][0]
		idmax= np.arange(0,len(freq),dtype=np.int32)[freq<fmax][-1]

		sed_at_fmin= sed[idmin-1] + (np.log10(fmin)-np.log10(freq[idmin-1])) * (sed[idmin]-sed[idmin-1])/(np.log10(freq[idmin])-np.log10(freq[idmin-1]))
		sed_at_fmax= sed[idmax]   + (np.log10(fmax)-np.log10(freq[idmax]))   * (sed[idmax+1]-sed[idmax])/(np.log10(freq[idmax+1])-np.log10(freq[idmax]))

		freq = np.append(np.append(fmin,freq[idmin:idmax+1]),fmax)
		sed  = np.append(np.append(sed_at_fmin,sed[idmin:idmax+1]),sed_at_fmax)

		return np.trapz(np.log(10)*10**sed, np.log10(freq))

	LHX =integrate( sedall, freqall, 2.*1000.*con.e.value/con.h.value, 10.*1000.*con.e.value/con.h.value)

	LSX =integrate( sedall, freqall, 0.5*1000.*con.e.value/con.h.value,2.*1000.*con.e.value/con.h.value)
	
	logLB  =np.mean(sedall[np.abs(lamball-4450.)<940./2])

	logLIR =np.mean(sedall[np.abs(lamball-150000.)<15000.])

	Lbol= integrate( sedall, freqall, con.c.value/(30.*1e-6), 500*1000.*con.e.value/con.h.value)
	return  np.log10(Lbol), np.log10(LHX), np.log10(LSX), logLB, logLIR

photon_index_list = np.genfromtxt('./xspec_lib/pindex_list.dat')

sigma_sl, sl_best = 0.125, 0.44
sigma_pi, pi_best = 0.2, 1.9
sigma_ox, ox_best = 0.1, 0

np.random.seed(4367)
Nsamples = 100000

par1=np.random.normal(sl_best, sigma_sl, size=Nsamples)
par2=np.random.normal(pi_best, sigma_pi, size=Nsamples)
par3=np.random.normal(ox_best, sigma_ox, size=Nsamples)
par5=np.random.uniform(-5,15,size=Nsamples)+L_solar
#par5=np.random.normal(10, 5 ,size=Nsamples)+L_solar

Lbols   = np.zeros( Nsamples )
HXcorrs = np.zeros( Nsamples )
SXcorrs = np.zeros( Nsamples )
Bcorrs  = np.zeros( Nsamples )
IRcorrs = np.zeros( Nsamples )

for i in range(Nsamples):
	if i%100==0: print i,'/',Nsamples
	Lbol,LHX,LSX,LB,LIR = returnall(par5[i],par1[i],par2[i],par3[i])				
	Lbols[i]   = Lbol
	HXcorrs[i] = Lbol-LHX
	SXcorrs[i] = Lbol-LSX
	Bcorrs[i]  = Lbol-LB
	IRcorrs[i] = Lbol-LIR

np.savetxt('./corr_ensemble/ensemble.dat', np.c_[Lbols, HXcorrs, SXcorrs, Bcorrs, IRcorrs])
np.savetxt('./corr_ensemble/ensemble_pars.dat', np.c_[par5, par1, par2, par3])