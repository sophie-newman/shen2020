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

from new_load_akiyama18_lf_data import *
from new_load_bongiorno07_lf_data import *
from new_load_croom09_lf_data import *
from new_load_fontanot07_lf_data import *
from new_load_giallongo15_lf_data import *
from new_load_gilkman11_lf_data import *
from new_load_ikeda12_lf_data import *
from new_load_jiang09_lf_data import *
from new_load_jiang16_lf_data import *
from new_load_kashikawa15_lf_data import *
from new_load_master12_lf_data import *
from new_load_matsuoka18_lf_data import *
from new_load_mcgreer13_lf_data import *
from new_load_mcgreer18_lf_data import *
from new_load_palanque13_lf_data import *
from new_load_palanque16_lf_data import *
from new_load_ross13_lf_data import *
from new_load_shen12_lf_data import *
from new_load_siana08_lf_data import *
from new_load_wang18_lf_data import *
from new_load_willott10_lf_data import *
from new_load_yang16_lf_data import *
from new_load_yang18_lf_data import *

from new_load_kk18_lf_shape import *

from newx_load_aird08_lf_data import *
from newx_load_aird10_lf_data import *
from newx_load_aird15_lf_data import *
from newx_load_aird15_b_lf_data import *
from newx_load_embero09_lf_data import *
from newx_load_fiore12_lf_data import *
from newx_load_khorunzhev18_lf_data import *
from newx_load_miyaji15_lf_data import *

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

new_dset_ids={
	#1450
	"AKIYAMA18":  -1,
	"BONGIORNO07":-1,
	"CROOM09":    -1,
	"FONTANOT07": -1,
	"GIALLONGO15":-1,
	"GILKMAN11":  -1,
	"IKEDA12":    -1,
	"JIANG09":    -1,
	"JIANG16":    -1,
	"KASHIKAWA15":-1,
	"MASTER12":   -1,
	"MATSUOKA18": -1,
	"MCGREER13":  -1,
	"MCGREER18":  -1,
	"PALANQUE13": -1,
	"PALANQUE16": -1,
	"ROSS13":     -1,
	"ROSS13_S82": -1,
	"SHEN12":     -1,
	"SIANA08":    -1,
	"WANG18":     -1,
	"WILLOTT10":  -1,
	"YANG16":     -1,
	"YANG18":     -1,
	#Xray H
	"AIRD08":      -4,
	"AIRD10":      -4,
	"AIRD15_SX":   -4,
	"AIRD15_HX":   -4,
	"AIRD15_b":    -4,
	"EBRERO09_SX": -4,
	"EBRERO09_HX": -4,
	"FIORE12":     -4,
	"KHORUNZHEV18":-4,
	"MIYAJI15":    -4
}

kk_dset_ids={
	#1450
	"BONGIORNO07":-1,
	"GIALLONGO15":-1,
	"GILKMAN11":  -1,
	"JIANG09":    -1,
	"KASHIKAWA15":-1,
	"MASTER12":   -1,
	"MCGREER13":  -1,
	"PALANQUE13": -1,
	"ROSS13":     -1,
	"ROSS13_S82": -1,
	"SIANA08":    -1,
	"WILLOTT10":  -1
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

new_zmins={
	#1450
	"AKIYAMA18":  np.array([3.6]),
	"BONGIORNO07":np.array([2.0,2.5,3.0]),
	"CROOM09":    np.array([0.4,0.68,1.06,1.44,1.82,2.20]),
	"FONTANOT07": np.array([3.5,4.0]),
	"GIALLONGO15":np.array([4.0,4.5,5.0]),
	"GILKMAN11":  np.array([3.8]),
	"IKEDA12":    np.array([4.57]),
	"JIANG09":    np.array([5.7]),
	"JIANG16":    np.array([5.7]),
	"KASHIKAWA15":np.array([5.85]),
	"MASTER12":   np.array([3.1,3.5]),
	"MATSUOKA18": np.array([5.7]),
	"MCGREER13":  np.array([4.7]),
	"MCGREER18":  np.array([4.7]),
	"PALANQUE13": np.array([2.2]),
	"PALANQUE16": np.array([0.68,1.06,1.44,1.82,2.20,2.60,3.00,3.50]),
	"ROSS13":     np.array([2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 3.0, 3.25]),
	"ROSS13_S82": np.array([2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 3.0, 3.25]),
	"SHEN12":     np.array([0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.4,2.9,3.5,4.0,4.5]),
	"SIANA08":    np.array([2.8]),
	"WANG18":     np.array([6.45]),
	"WILLOTT10":  np.array([5.75]),
	"YANG16":     np.array([4.7]),
	"YANG18":     np.array([0.5,1.0,1.5,2.0,2.5,3.0]),#,3.5,4.0]),
	#HX
	"AIRD08":      np.array([2.5]),
	"AIRD10":      np.array([0.0,0.2,0.5,0.8,1.0,1.2,1.5,2.0,2.5]),
	"AIRD15_SX":   np.array([0.01,0.20,0.40,0.60,0.80,1.00,1.20,1.50,2.00,2.50,3.50,5.00]),
	"AIRD15_HX":   np.array([0.01,0.20,0.40,0.60,0.80,1.00,1.20,1.50,2.00,2.50,3.50,5.00]),
	"AIRD15_b":    np.array([0.1,0.5,1.0]),
	"EBRERO09_SX": np.array([0.01,0.5,1.0,2.0]),
	"EBRERO09_HX": np.array([0.01,0.5,1.0,2.0]),
	"FIORE12":     np.array([3.,4.,5.8]),
	"KHORUNZHEV18":np.array([3.,3.19,3.47,3.90,4.30]),
	"MIYAJI15":    np.array([0.015,0.20,0.40,0.80,1.60,2.30])
}
new_zmaxs={
	#1450
	"AKIYAMA18":  np.array([4.3]),
	"BONGIORNO07":np.array([2.5,3.0,4.0]),
	"CROOM09":    np.array([0.68,1.06,1.44,1.82,2.20,2.60]),
	"FONTANOT07": np.array([4.0,5.2]),
	"GIALLONGO15":np.array([4.5,5.0,6.5]),
	"GILKMAN11":  np.array([5.2]),
	"IKEDA12":    np.array([5.57]),
	"JIANG09":    np.array([6.6]),
	"JIANG16":    np.array([6.4]),
	"KASHIKAWA15":np.array([6.45]),
	"MASTER12":   np.array([3.5,5.0]),
	"MATSUOKA18": np.array([6.5]),
	"MCGREER13":  np.array([5.1]),
	"MCGREER18":  np.array([5.4]),
	"PALANQUE13": np.array([2.6]),
	"PALANQUE16": np.array([1.06,1.44,1.82,2.20,2.60,3.00,3.50,4.00]),
	"ROSS13":     np.array([2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 3.0, 3.25, 3.5]),
	"ROSS13_S82": np.array([2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 3.0, 3.25, 3.5]),
	"SHEN12":     np.array([0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.4,2.9,3.5,4.0,4.5,5.0]),
	"SIANA08":    np.array([3.5]),
	"WANG18":     np.array([7.05]),
	"WILLOTT10":  np.array([6.45]),
	"YANG16":     np.array([5.4]),
	"YANG18":     np.array([1.0,1.5,2.0,2.5,3.0,3.5]),#,4.0,4.5])
	#HX
	"AIRD08":      np.array([3.5]),
	"AIRD10":      np.array([0.2,0.5,0.8,1.0,1.2,1.5,2.0,2.5,3.0]),
	"AIRD15_SX":   np.array([0.20,0.40,0.60,0.80,1.00,1.20,1.50,2.00,2.50,3.50,5.00,7.00]),
	"AIRD15_HX":   np.array([0.20,0.40,0.60,0.80,1.00,1.20,1.50,2.00,2.50,3.50,5.00,7.00]),
	"AIRD15_b":    np.array([0.5,1.0,3.0]),
	"EBRERO09_SX": np.array([0.5,1.0,2.0,3.0]),
	"EBRERO09_HX": np.array([0.5,1.0,2.0,3.0]),
	"FIORE12":     np.array([4.,5.,7.5]),
	"KHORUNZHEV18":np.array([3.19,3.47,3.90,4.30,5.10]),
	"MIYAJI15":    np.array([0.20,0.40,0.80,1.60,2.30,4.60])
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

new_load_LF_data={
	#1450
	"AKIYAMA18":  load_akiyama18_lf_data,
	"BONGIORNO07":load_bongiorno07_lf_data,
	"CROOM09":    load_croom09_lf_data,
	"FONTANOT07": load_fontanot07_lf_data,
	"GIALLONGO15":load_giallongo15_lf_data,
	"GILKMAN11":  load_gilkman11_lf_data,
	"IKEDA12":    load_ikeda12_lf_data,
	"JIANG09":    load_jiang09_lf_data,
	"JIANG16":    load_jiang16_lf_data,
	"KASHIKAWA15":load_kashikawa15_lf_data,
	"MASTER12":   load_master12_lf_data,
	"MATSUOKA18": load_matsuoka18_lf_data,
	"MCGREER13":  load_mcgreer13_lf_data,
	"MCGREER18":  load_mcgreer18_lf_data,
	"PALANQUE13": load_palanque13_lf_data,
	"PALANQUE16": load_palanque16_lf_data,
	"ROSS13":     load_ross13_lf_data,
	"ROSS13_S82": load_ross13_s82_lf_data,
	"SHEN12":     load_shen12_lf_data,
	"SIANA08":    load_siana08_lf_data,
	"WANG18":     load_wang18_lf_data,
	"WILLOTT10":  load_willott10_lf_data,
	"YANG16":     load_yang16_lf_data,
	"YANG18":     load_yang18_lf_data,
	#Xray H
	"AIRD08":      load_aird08_lf_data,
	"AIRD10":      load_aird10_lf_data,
	"AIRD15_SX":   load_aird15_lf_data_softsel,
	"AIRD15_HX":   load_aird15_lf_data_hardsel,
	"AIRD15_b":    load_aird15_b_lf_data,
	"EBRERO09_SX": load_ebrero09_SX_lf_data,
	"EBRERO09_HX": load_ebrero09_HX_lf_data,
	"FIORE12":     load_fiore12_lf_data,
	"KHORUNZHEV18":load_khorunzhev18_lf_data_vitobins,
	"MIYAJI15":    load_miyaji15_lf_data
}
