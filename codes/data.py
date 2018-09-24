import numpy as np
import os, errno
import __main__
import sys

homepath=""
sys.path.append()
sys.path.append()
sys.path.append()
sys.path.append()

#COSMOLOGY
hubble      = 0.7
Omega0      = 0.3
OmegaBaryon = 0.05
OmegaLambda = 0.7

#KEYS
KEYS={   
	"LF_model": "Fiducial",
	"Extinction": "Fiducial",
	"Fit_method": "Chisq"
}
Z_NORM_SET
FIT_LOG_Z

#LUMINOSITY GRID
L_bol_grid = np.linspace(8.,18.,101)+0.008935
d_log_l_bol = L_bol[1]-L_bol[0]

#
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

