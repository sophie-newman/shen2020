from data import *
from load_2slaq_lf_data import *
from load_barger_lf_data import *
from load_beckmann_lf_data import *#
from load_brown_lf_data import *
from load_combo17_lf_data import *
from load_cristiani_lf_data import *
from load_hao_lf_data import *
from load_hasinger_lf_data import *#
from load_huchraburg_lf_data import *#
from load_hunt_lf_data import *
from load_jiang_lf_data import * #
from load_kdc_lf_data import *
from load_lafranca_lf_data import *
from load_matute_lf_data import *
from load_miyaji_lf_data import *
from load_nandra_lf_data import *
from load_sazrev_lf_data import *
from load_sdss_dr3_lf_data import *
from load_sdss_fan_lf_data import *
from load_shinozaki_lf_data import *#
from load_siana_lf_data import *
from load_silverman_sx_lf_data import *
from load_silverman_hx_lf_data import *
from load_ssg_lf_data import *
from load_ueda_lf_data import *

from new_load_bongiorno07_lf_data import *
from new_load_giallongo15_lf_data import *
from new_load_gilkman11_lf_data import *
from new_load_jiang09_lf_data import *
from new_load_kashikawa15_lf_data import *
from new_load_master12_lf_data import *
from new_load_mcgreer13_lf_data import *
from new_load_mcgreer18_lf_data import *
from new_load_palanque13_lf_data import *
from new_load_ross13_lf_data import *
from new_load_siana08_lf_data import *
from new_load_wang18_lf_data import *
from new_load_willott10_lf_data import *
from new_load_yang18_lf_data import *

dset_ids={
	#Bband
	"DR3": -1,
	"2SLAQ": -1,
	"HUNT": -1,
	"SIANA": -1,
	"CRISTIANI": -1,
	"KDC": -1,
	"FAN": -1,
	"SSG": -1,
	"COMBO17": -1,
	#soft xray
	"HASINGER": -3,
	"MIYAJI": -3,
	"SILVERMAN_SX": -3,
	#hard xray
	"UEDA": -4,
	"LAFRANCA": -4,
	"SILVERMAN_HX": -4,
	"BARGER": -4,
	"NANDRA": -4,
	"SAZREV": -4,
	#INFRARED
	"BROWN": -2,
	"MATUTE": -2,
	#narrow lines
	"HAO": -4
}

zmins={
	"DR3": np.array([0.40, 0.68, 1.06, 1.44, 1.82, 2.21, 2.60, 3.03, 3.50, 4.00, 4.50]),
	"2SLAQ":np.array([0.40, 0.68, 0.97, 1.25, 1.53, 1.81]),
	"HUNT":np.array([2.0]),
	"SIANA":np.array([2.9]),
	"CRISTIANI":np.array([4.0]),
	"KDC":np.array([4.0]),
	"FAN":np.array([3.6, 3.9, 4.4, 5.8]),
	"SSG":np.array([2.75]),
	"COMBO17":np.array([1.2, 1.8, 2.4, 3.0, 3.6, 4.2]),
	#soft xray
	"HASINGER":np.array([0.0, 0.2, 0.4, 0.8, 1.6, 3.2]),
	"MIYAJI":np.array([0.0, 0.2, 0.4, 0.8, 1.6, 2.3]),
	#hard xray
	"UEDA":np.array([0.0, 0.2, 0.4, 0.8, 1.6]),
	"LAFRANCA":np.array([0.0, 0.5, 1.0, 1.5, 2.5]),
	"SILVERMAN_SX":np.array([0.5, 1.0, 2.0, 4.0]),
	"SILVERMAN_HX":np.array([0.2, 1.5, 3.0]),
	"BARGER":np.array([0.1, 0.4, 0.8, 1.5, 3.0]),
	"NANDRA":np.array([2.75]),
	"SAZREV":np.array([0.0]),
	#INFRARED
	"BROWN":np.array([1.5]),
	"MATUTE":np.array([0.0, 0.5]),
	#narrow lines
	"HAO":np.array([0.00])
}
zmaxs={
	"DR3": np.array([0.68, 1.06, 1.44, 1.82, 2.21, 2.60, 3.03, 3.50, 4.00, 4.50, 5.00]),
	"2SLAQ":np.array([0.68, 0.97, 1.25, 1.53, 1.81, 2.10]),
	"HUNT":np.array([4.0]),
	"SIANA":np.array([3.4]),
	"CRISTIANI":np.array([5.2]),
	"KDC":np.array([4.5]),
	"FAN":np.array([3.9, 4.4, 5.0, 6.2]),
	"SSG":np.array([4.75]),
	"COMBO17":np.array([1.8, 2.4, 3.0, 3.6, 4.2, 4.8]),
	"HASINGER":np.array([0.2, 0.4, 0.8, 1.6, 3.2, 4.8]),
	"MIYAJI":np.array([0.2, 0.4, 0.8, 1.6, 2.3, 4.6]),
	"UEDA":np.array([0.2, 0.4, 0.8, 1.6, 3.0]),
	"LAFRANCA":np.array([0.5, 1.0, 1.5, 2.5, 3.5]),
	"SILVERMAN_SX":np.array([1.0, 1.5, 3.0, 5.5]),
	"SILVERMAN_HX":np.array([0.5, 2.0, 4.0]),
	"BARGER":np.array([0.4, 0.8, 1.2, 3.0, 5.0]),
	"NANDRA":np.array([3.25]),
	"SAZREV":np.array([0.1]),
	"BROWN":np.array([2.5]),
	"MATUTE":np.array([0.2, 2.0]),
	"HAO":np.array([0.10])
}
load_LF_data={
	"DR3":load_sdss_dr3_lf_data,
	"2SLAQ":load_2slaq_lf_data,
	"HUNT":load_hunt_lf_data,
	"SIANA":load_siana_lf_data,
	"CRISTIANI":load_cristiani_lf_data,
	"KDC":load_kdc_lf_data,
	"FAN":load_sdss_fan_lf_data,
	"SSG":load_ssg_lf_data,
	"COMBO17":load_combo17_lf_data,
	"HASINGER":load_hasinger_lf_data,
	"MIYAJI":load_miyaji_lf_data,
	"UEDA":load_ueda_lf_data,
	"LAFRANCA":load_lafranca_lf_data,
	"SILVERMAN_SX":load_silverman_sx_lf_data,
	"SILVERMAN_HX":load_silverman_hx_lf_data,
	"BARGER":load_barger_lf_data,
	"NANDRA":load_nandra_lf_data,
	"SAZREV":load_sazrev_lf_data,
	"BROWN":load_brown_lf_data,
	"MATUTE":load_matute_lf_data,
	"HAO":load_hao_lf_data
}
return_LF={
	"DR3":return_sdss_dr3_lf_fitted,
	"2SLAQ":return_2slaq_lf_fitted,
	"HUNT":return_combo17_lf_fitted,
	"SIANA":return_combo17_lf_fitted,
	"CRISTIANI":return_combo17_lf_fitted,
	"KDC":return_sdss_dr3_lf_fitted,
	"FAN":return_sdss_fan_lf_fitted,
	"SSG":return_sdss_dr3_lf_fitted,
	"COMBO17":return_combo17_lf_fitted,
	"HASINGER":return_hasinger_lf_fitted,
	"MIYAJI":return_miyaji_lf_fitted,
	"UEDA":return_ueda_lf_fitted,
	"LAFRANCA":return_lafranca_lf_fitted,
	"SILVERMAN_SX":return_hasinger_lf_fitted,
	"SILVERMAN_HX":return_ueda_lf_fitted,
	"BARGER":return_ueda_lf_fitted,
	"NANDRA":None,
	"SAZREV":return_ueda_lf_fitted,
	"BROWN":return_brown_lf_fitted,
	"MATUTE":return_matute_lf_fitted,
	"HAO":return_ueda_lf_fitted
}
