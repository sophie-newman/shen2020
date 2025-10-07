import sys 

######### CONFIG #######
with_clib = True      # needed if you want to convert bolometric luminosities to band luminosities
with_obs_data = True  # needed if you want to load observational data
if with_obs_data: with_clib = True

# home directory path,  PLEASE modify it to the path on your machine
homepath="/mnt/home/snewman/quasarlf/pubtools"
datapath=homepath+"/data/"

# path to the observational data, required only if you want to load observational data
if with_obs_data:
	sys.path.append(homepath+"/obdata_copy/")
