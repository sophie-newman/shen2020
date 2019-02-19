#!/bin/sh
#SBATCH -p productionQ
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --exclusive
#SBATCH --mail-user=xshen@caltech.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mem-per-cpu=4000
#SBATCH -t 1-00:00           # Runtime in D-HH:MM

mpirun -n 16 python test.py 3 > output/test_log.txt

