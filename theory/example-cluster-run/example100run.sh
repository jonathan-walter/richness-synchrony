#!/bin/bash -l
#PBS -N GeoSynch100
#PBS -l walltime=02:00:00,nodes=1:ppn=24,mem=250gb
#PBS -m abe
#PBS -M lshoema1@uwyo.edu
#PBS -j oe

# Change to the relevant working directory
cd ~/GeoSynchSimMesabi

# Load R and MPI
module load R/3.4.4
module load ompi/3.0.0/gnu-7.2.0-centos7

export RMPI_TYPE=OPENMPI
export OMPI_MCA_mpi_warn_on_fork=0

# Run the R script itself, saving the output to
#	a log file
mpirun -np 1 R CMD BATCH --no-restore --no-save --quiet example100run.R example100run.log 
