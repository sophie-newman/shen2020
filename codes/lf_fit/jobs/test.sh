#!/bin/sh
#SBATCH -p productionQ
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --mem-per-cpu=4000
#SBATCH -t 1-00:00           # Runtime in D-HH:MM

module purge

module load python/anaconda2-4.1.1
source activate customenv
module load gcc/7.3.0
module load openmpi/2.0.1
module load hdf5/1.8.17

module list
which srun
which mpiexec
which mpirun
