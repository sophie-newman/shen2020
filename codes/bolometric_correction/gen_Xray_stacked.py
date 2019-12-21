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

photon_index_list = np.genfromtxt('./xspec_lib/pindex_list.dat')

def returnXray_sub(f2kev, gamma):
        #interpolate between precalculated X-ray SEDs/ extrapolate if outside
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

	return freq3, 10**sed3
        #return lamb3[lamb3<50], sed3[lamb3<50]

def returnXray():
	np.random.seed(4367)
	Nsamples = 1000
	pindexes = np.random.normal(1.9 , 0.2, size=Nsamples) #photon index
	for i in range(Nsamples):
		freq, sed = returnXray_sub( 1e-3 ,  pindexes[i])
		if i==0:
			sed_stacked = sed.copy()
		else:
			sed_stacked = sed_stacked + sed

	return freq, np.log10(sed_stacked) - 15.

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
