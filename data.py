import numpy as np
import os, errno
import __main__
import sys
from astropy.cosmology import FlatLambdaCDM

homepath="/Users/xuejianshen/Desktop/QuasarLF/git/"
sys.path.append(homepath+"codes/bolometric_correction/")
sys.path.append(homepath+"codes/convolution/")
sys.path.append(homepath+"codes/lf_fit/")
sys.path.append(homepath+"codes/possessing_obdata/")
sys.path.append(homepath+"codes/")

datapath=homepath+"data/"

#COSMOLOGY
hubble      = 0.7
Omega0      = 0.3
OmegaBaryon = 0.04
OmegaLambda = 0.7
cosmo        = FlatLambdaCDM(H0=hubble*100, Om0=Omega0)

#CONSTANTS
M_sun_Bband_AB = 4.77  #4.72

#KEYS
KEYS={   
	"LF_model": "Fiducial",
	"Extinction_model": "Fiducial",
	"Fit_method": "Chisq",
	"Fit_key": 0
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

# function to correct QSO luminosities from the old-school 
#    Omega_M = 1, q=0.5 cosmology to the modern 0.7/0.3 cosmology
#
cosmo_old        = FlatLambdaCDM(H0=0.5*100, Om0=1)
def lum_correct_cosmology(redshift):
    dlum_old = cosmo_old.luminosity_distance(redshift).value
    dlum_new = cosmo.luminosity_distance(redshift).value
    distance_corr_factor = (dlum_new/dlum_old)**2
    return distance_corr_factor

# function to correct QSO number densities from the old-school 
#    Omega_M = 1, q=0.5 cosmology to the modern 0.7/0.3 cosmology
#
def phi_correct_cosmology(redshift):
    vol_old = cosmo_old.comoving_volume(redshift).value
    vol_new = cosmo.comoving_volume(redshift).value
    volume_corr_factor = vol_old/vol_new
    return volume_corr_factor
