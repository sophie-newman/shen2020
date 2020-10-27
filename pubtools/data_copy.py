import numpy as np
import os, errno
import __main__
import sys
from astropy.cosmology import FlatLambdaCDM
import astropy.constants as con

# path to the observational data, required only if you want to load observational data
homepath="/home/xuejian/works/quasarLF/git/"
sys.path.append(homepath+"codes/obdata/")

#COSMOLOGY
hubble      = 0.7
Omega0      = 0.3
OmegaBaryon = 0.04
OmegaLambda = 0.7
cosmo        = FlatLambdaCDM(H0=hubble*100, Om0=Omega0)

#CONSTANTS
M_sun_Bband_AB = 4.77  #4.72
L_solar = np.log10(3.9e33)  #4e33
Fab = 3631.*1e-23*4*np.pi*(10.*con.pc.value*100)**2 #erg/s/Hz
#L_solar=np.log10(con.L_sun.value*1e7) #log10 of solar luminosity in erg/s

#KEYS (redundant)
KEYS={
	"LF_model": "Fiducial",  #LDDE  PLE   Modified_Schechter
	"Extinction_model": "Fiducial",  #Fixed   H07
	"SED_model": "Fiducial",
	"Fit_method": "Chisq", #Modified_Chisq
	"FIT_KEY": 0
}

#LUMINOSITY GRID
L_bol_grid = np.linspace(8.,18.,101)[:-1]#+0.008935
N_bol_grid = len(L_bol_grid)
d_log_l_bol = L_bol_grid[1]-L_bol_grid[0]

L_bol_grid_for_CXB = np.linspace(6.,18.,121)[:-1]
N_bol_grid_for_CXB = len(L_bol_grid_for_CXB)

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

def lum_correct_cosmo_flexible(redshift,h_old,Om0_old):
    cosmo_to_correct = FlatLambdaCDM(H0=h_old*100, Om0=Om0_old)
    dlum_old = cosmo_to_correct.luminosity_distance(redshift).value
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

def phi_correct_cosmo_flexible(redshift,h_old,Om0_old):
    cosmo_to_correct = FlatLambdaCDM(H0=h_old*100, Om0=Om0_old)
    vol_old = cosmo_to_correct.differential_comoving_volume(redshift).value
    vol_new = cosmo.differential_comoving_volume(redshift).value
    volume_corr_factor = vol_old/vol_new
    return volume_corr_factor

# return the dust-to-gas ratio inferred by the MZR in Ma2016
# the dtg here is only a normalized value that should be feeded in the convolution C code
def return_dtg(z):
	dtg_mw=0.78
	def MZR_Ma2016(z):
		logM = 11
		logZ = 0.35*(logM-10) + 0.93*np.exp(-0.43*z)  - 1.05
		return 10**logZ

	return MZR_Ma2016(z)/MZR_Ma2016(0) * dtg_mw

def return_miyaji15_lf_fitted(Llist,z):
        Llist = Llist + L_solar
        p1_44,p2_44,p3 = 5.29, -0.35, -5.6
        zb1_44,zb2 = 1.1, 2.7
        alpha = 0.18
        beta1, beta2 = 1.2, 1.5

        gamma1 = 1.17
        gamma2 = 2.80
        lbreak = 44.04
        phi_s = 1.56e-6

        def ed(Lx,z):
                p1 = p1_44 + beta1*(Lx - 44.)
                p2 = p2_44 + beta2*(Lx - 44.)

                zb1_0 = zb1_44 / np.power( 10**44/10**44.5, alpha )
                if Lx <= 44.5: zb1 = zb1_0 * np.power( 10**Lx/10**44.5, alpha )
                else: zb1 = zb1_0

                if z < zb1:  ed = (1.+z)**p1
                elif z < zb2: ed = (1.+zb1)**p1 * np.power( (1+z)/(1+zb1) ,p2)
                else: ed = (1.+zb1)**p1 * np.power( (1+zb2)/(1+zb1) ,p2) * np.power( (1+z)/(1+zb2) ,p3)

                return ed

        x = np.power(10., Llist-lbreak)
        eds = np.zeros(len(Llist))
        for i in range(len(eds)):
                eds[i] = ed(Llist[i],z)
        return np.log10( phi_s/( np.power(x,gamma1) + np.power(x,gamma2))*eds )
