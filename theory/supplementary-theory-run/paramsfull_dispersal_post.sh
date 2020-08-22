#!/bin/bash -l

#SBATCH --account=coexistence
#SBATCH --time=04:00:00
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lash1937@gmail.com
#SBATCH --job-name=rich_synch_all

# Change to the relevant working directory
cd /project/coexistence/lshoema1/richness_synchrony/

# Load R and MPI
module load gcc/7.3.0 r/3.5.3 openmpi/3.1.0 r-rmpi/0.6-9-r353-py27

mpirun -np 1 R CMD BATCH --no-restore --no-save --quiet paramsfull_dispersal_post.R paramsfull_dispersal_post.log 
