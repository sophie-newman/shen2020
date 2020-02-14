#!/bin/sh
#SBATCH -p productionQ
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=12
#SBATCH --exclusive
#SBATCH --mail-user=xshen@caltech.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mem-per-cpu=4000
#SBATCH -t 1-00:00           # Runtime in D-HH:MM

cd ..
mpirun -n 96 python lf_fitter_shallowfaint.py &> gfit_shallowfaint.log
