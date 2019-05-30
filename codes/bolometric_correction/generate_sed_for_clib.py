from data import *
import numpy as np 
import astropy.constants as con
from scipy.interpolate import interp1d

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

	return freq1, sed1+totscale

def returnXray():
	data3=np.genfromtxt("xspec_lib/Xspec_1_9.dat",names=["E","f"],)	
	freq3=data3['E']*1000.*con.e.value/con.h.value
	lamb3=con.c.value/freq3 * 1e10
	sed3 =np.log10(data3['f']*data3['E'])

	lamb2kev = con.c.value/ (2.*1000.*con.e.value/con.h.value) *1e10
	freq2kev = 2.*1000.*con.e.value/con.h.value

	return freq3, sed3

freq_xray=np.array([ 16.00, 16.02, 16.04, 16.06, 16.08, 16.10, 16.12, 16.14, 16.16, 16.18, 16.20, 16.22, 16.24, 16.26, 16.28, 16.30, 16.32, 16.34, 16.36, 16.38, 16.40, 16.42, 16.44, 16.46, 16.48, 16.50,
16.52, 16.54, 16.56, 16.58, 16.60, 16.62, 16.64, 16.66, 16.68, 16.70, 16.72, 16.74, 16.76, 16.78, 16.80, 16.82, 16.84, 16.86, 16.88, 16.90, 16.92, 16.94, 16.96, 16.98, 17.00, 17.02,
17.04, 17.06, 17.08, 17.10, 17.12, 17.14, 17.16, 17.18, 17.20, 17.22, 17.24, 17.26, 17.28, 17.30, 17.32, 17.34, 17.36, 17.38, 17.40, 17.42, 17.44, 17.46, 17.48, 17.50, 17.52, 17.54,
17.56, 17.58, 17.60, 17.62, 17.64, 17.66, 17.68, 17.70, 17.72, 17.74, 17.76, 17.78, 17.80, 17.82, 17.84, 17.86, 17.88, 17.90, 17.92, 17.94, 17.96, 17.98, 18.00, 18.02, 18.04, 18.06,
18.08, 18.10, 18.12, 18.14, 18.16, 18.18, 18.20, 18.22, 18.24, 18.26, 18.28, 18.30, 18.32, 18.34, 18.36, 18.38, 18.40, 18.42, 18.44, 18.46, 18.48, 18.50, 18.52, 18.54, 18.56, 18.58,
18.60, 18.62, 18.64, 18.66, 18.68, 18.70, 18.72, 18.74, 18.76, 18.78, 18.80, 18.82, 18.84, 18.86, 18.88, 18.90, 18.92, 18.94, 18.96, 18.98, 19.00, 19.02, 19.04, 19.06, 19.08, 19.10,
19.12, 19.14, 19.16, 19.18, 19.20, 19.22, 19.24, 19.26, 19.28, 19.30, 19.32, 19.34, 19.36, 19.38, 19.40, 19.42, 19.44, 19.46, 19.48, 19.50, 19.52, 19.54, 19.56, 19.58, 19.60, 19.62,
19.64, 19.66, 19.68, 19.70, 19.72, 19.74, 19.76, 19.78, 19.80, 19.82, 19.84, 19.86, 19.88, 19.90, 19.92, 19.94, 19.96, 19.98, 20.00, 20.02, 20.04, 20.06, 20.08, 20.10, 20.12, 20.14,
20.16, 20.18, 20.20, 20.22, 20.24, 20.26, 20.28, 20.30, 20.32, 20.34, 20.36, 20.38, 20.40, 20.42, 20.44, 20.46, 20.48, 20.50, 20.52, 20.54, 20.56, 20.58, 20.60, 20.62, 20.64, 20.66,
20.68])
freq_temp, sed_temp = returnXray()
print np.log10(integrate( sed_temp, freq_temp, 2.*1000.*con.e.value/con.h.value, 10.*1000.*con.e.value/con.h.value))
print len(freq_xray)
print interp1d(np.log10(freq_temp), sed_temp)(freq_xray)

freq_bband=np.array([12.50, 12.52, 12.54, 12.56, 12.58, 12.60, 12.62, 12.64, 12.66, 12.68, 12.70, 12.72, 12.74, 12.76, 12.78, 12.80, 12.82, 12.84, 12.86, 12.88, 12.90, 12.92, 12.94, 12.96, 12.98, 13.00,
13.02, 13.04, 13.06, 13.08, 13.10, 13.12, 13.14, 13.16, 13.18, 13.20, 13.22, 13.24, 13.26, 13.28, 13.30, 13.32, 13.34, 13.36, 13.38, 13.40, 13.42, 13.44, 13.46, 13.48, 13.50, 13.52,
13.54, 13.56, 13.58, 13.60, 13.62, 13.64, 13.66, 13.68, 13.70, 13.72, 13.74, 13.76, 13.78, 13.80, 13.82, 13.84, 13.86, 13.88, 13.90, 13.92, 13.94, 13.96, 13.98, 14.00, 14.02, 14.04,
14.06, 14.08, 14.10, 14.12, 14.14, 14.16, 14.18, 14.20, 14.22, 14.24, 14.26, 14.28, 14.30, 14.32, 14.34, 14.36, 14.38, 14.40, 14.42, 14.44, 14.46, 14.48, 14.50, 14.52, 14.54, 14.56,
14.58, 14.60, 14.62, 14.64, 14.66, 14.68, 14.70, 14.72, 14.74, 14.76, 14.78, 14.80, 14.82, 14.84, 14.86, 14.88, 14.90, 14.92, 14.94, 14.96, 14.98, 15.00, 15.02, 15.04, 15.06, 15.08,
15.10, 15.12, 15.14, 15.16, 15.18, 15.20, 15.22, 15.24, 15.26, 15.28, 15.30, 15.32, 15.34, 15.36, 15.38, 15.40, 15.42, 15.44, 15.46, 15.48, 15.50, 15.52, 15.54, 15.56, 15.58, 15.60,
15.62, 15.64, 15.66, 15.68])
freq_temp, sed_temp = returnIR_to_UV(45.)
print tophat( sed_temp, freq_temp, con.c.value/((4450+470)*1e-10), con.c.value/((4450-470)*1e-10))
print len(freq_bband)
print interp1d(np.log10(freq_temp), sed_temp)(freq_bband)

