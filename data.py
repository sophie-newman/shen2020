import numpy as np
import os, errno
import __main__
import sys

homepath="/Users/xuejianshen/Desktop/QuasarLF/git/"
sys.path.append(homepath+"codes/bolometric_correction/")
sys.path.append(homepath+"codes/convolution/")
sys.path.append(homepath+"codes/lf_fit/")
sys.path.append(homepath+"codes/possessing_obdata/")
sys.path.append(homepath+"codes/")

#COSMOLOGY
hubble      = 0.7
Omega0      = 0.3
OmegaBaryon = 0.05
OmegaLambda = 0.7

#KEYS
KEYS={   
	"LF_model": "Fiducial",
	"Extinction": "Fiducial",
	"Fit_method": "Chisq",
	"FIT_KEY": 0
}
#Z_NORM_SET
#FIT_LOG_Z

#LUMINOSITY GRID
L_bol_grid = np.linspace(8.,18.,11)+0.008935
d_log_l_bol = L_bol_grid[1]-L_bol_grid[0]

#functions
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

