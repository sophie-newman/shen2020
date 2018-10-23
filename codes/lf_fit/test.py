from data import *
from ctypes import *
import numpy as np
import ctypes

#load the shared object file
c_extenstion = CDLL(homepath+'codes/c_lib/convolve.so')

convolve_c = c_extenstion.convolve
convolve_c.restype = ctypes.POINTER(ctypes.c_double * 100)

redshift=c_double(2.)
nu=c_double(-1.)
input_c= L_bol_grid.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

res = convolve_c(input_c,redshift,nu)
res = [i for i in res.contents]
res = np.array(res,dtype=np.float64)
print res
