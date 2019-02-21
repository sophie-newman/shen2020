#!/bin/sh
#SBATCH -p productionQ
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=6
#SBATCH --exclusive
#SBATCH --mail-user=xshen@caltech.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mem-per-cpu=8000
#SBATCH -t 1-00:00           # Runtime in D-HH:MM

mpirun -n 48 python lf_fitter.py > main_log.txt
