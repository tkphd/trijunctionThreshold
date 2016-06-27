#!/bin/bash
# SLURM batch script for trijunction drag grain growth simulations:
# generate initial condition for ideal and drag comparison

# Make sure your binary has been compiled before launching!
# Usage: sbatch AMOS_ideal_medium_----.sh

# Declare common SLURM settings

# Working directory and output file (datafiles  and stdout redirect)
#SBATCH -D /gpfs/u/scratch/GGST/GGSTkllt/trijunctionThreshold/ideal/medium
#SBATCH -o /gpfs/u/barn/GGST/GGSTkllt/trijunctionThreshold/dat/ideal/AMOS_medium_init.log

# Cluster partition and job size
#SBATCH --partition small
#SBATCH --overcommit
#SBATCH --time 1440
#SBATCH --nodes 64
#SBATCH --ntasks-per-node 32
#SBATCH -J medium_4096

# When and to whom notifications should be sent
#SBATCH --mail-type=END
#SBATCH --mail-user=kellet@rpi.edu

if [[ ! -d $SLURM_SUBMIT_DIR ]]
then
	mkdir -p $SLURM_SUBMIT_DIR
fi

SRCDIR=/gpfs/u/barn/GGST/GGSTkllt/trijunctionThreshold/src/ideal/medium

if [[ ! -f $SRCDIR/q_GG.out ]]
then
	echo "Error: ${SRCDIR}/q_GG.out not found: cd ${SRCDIR} && make bgq"
	exit
fi

srun --runjob-opts="--mapping TEDCBA" $SRCDIR/./q_GG.out --example 2 qmedium.dat
