import numpy as np
import os, errno
import __main__
import sys
from astropy.cosmology import FlatLambdaCDM
import astropy.constants as con

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
L_solar=np.log10(3.9e33)
#L_solar=np.log10(con.L_sun.value*1e7) #log10 of solar luminosity in erg/s

#KEYS
KEYS={   
	"LF_model": "Fiducial",
	"Extinction_model": "Fiducial",
	"Fit_method": "Chisq",
	"FIT_KEY": 0
}
FIT_LOG_Z = 1
Z_NORM_SET = 2.0

#LUMINOSITY GRID
L_bol_grid = np.linspace(8.,18.,101)[:-1]+0.008935
d_log_l_bol = L_bol_grid[1]-L_bol_grid[0]

#functions
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

# returns the distance modulus for the cosmology adopted in HRH06 :: 
#   Omega_Matter = 0.3, Omega_Lambda = 0.7, h = 0.7
#   -- this function can easily be changed to return the appopriate distance 
#      modulus for any cosmology, but remember that the *observed* magnitudes 
#      are the unchanged quantities (i.e. changing this, one needs to change 
#      the luminosities and densities of the QLF accordingly)
#   -- for convenience, use a fitted function to the distance modulus for 
#      this cosmology, which is good to better than a factor of 0.0001 
#      at z<=10 and 0.001 at 10 < z < 30 (i.e. no significant source of error here)
def distance_modulus(z):
    P0,P1,P2,P3,P4,P5,P6,P7,P8=-36.840151,-13.167313,0.39803874,-0.33694804,2.5504298,0.56714463,-0.50247992,-0.044725351,-0.80476843
    x = np.log10(1.+z)
    d = P0 + P1*np.power(x,P2) + P3*np.power(x,P4) + P5*np.power(x,P6) + P7*np.power(x,P8)
    return d

# function to correct QSO luminosities from the old-school 
#    Omega_M = 1, q=0.5 cosmology to the modern 0.7/0.3 cosmology
# note that the hubble constants change has already been corrected in observation data loading part
cosmo_old        = FlatLambdaCDM(H0=hubble*100, Om0=1)
def lum_correct_cosmology(redshift):
    dlum_old = cosmo_old.luminosity_distance(redshift).value
    dlum_new = cosmo.luminosity_distance(redshift).value
    distance_corr_factor = (dlum_new/dlum_old)**2
    return distance_corr_factor

# function to correct QSO number densities from the old-school 
#    Omega_M = 1, q=0.5 cosmology to the modern 0.7/0.3 cosmology
#
def phi_correct_cosmology(redshift):
    vol_old = cosmo_old.differential_comoving_volume(redshift).value
    vol_new = cosmo.differential_comoving_volume(redshift).value
    volume_corr_factor = vol_old/vol_new
    return volume_corr_factor
