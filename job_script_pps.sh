#!/bin/bash

#SBATCH --time=10:00:00   # walltime
#SBATCH --ntasks=64   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8192M   # memory per CPU core
#SBATCH -J "convergence_pps"   # job name
#SBATCH --mail-user=timothy@daymachine.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_ON_NODE


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load julia
module load python

julia -t 64 --project="." ./examples/convergence_study_pps.jl
