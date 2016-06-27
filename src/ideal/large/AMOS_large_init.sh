#!/bin/bash
# SLURM batch script for trijunction drag grain growth simulations:
# generate initial condition for ideal and drag comparison

# Make sure your binary has been compiled before launching!
# Usage: sbatch AMOS_ideal_large_----.sh

# Declare common SLURM settings

# Working directory and output file (datafiles  and stdout redirect)
#SBATCH -D /gpfs/u/scratch/GGST/GGSTkllt/trijunctionThreshold/ideal/large
#SBATCH -o /gpfs/u/barn/GGST/GGSTkllt/trijunctionThreshold/dat/ideal/AMOS_large_init.log

# Cluster partition and job size
#SBATCH --partition medium
#SBATCH --overcommit
#SBATCH --time 720
#SBATCH --nodes 256
#SBATCH --ntasks-per-node 32
#SBATCH -J large_8192

# When and to whom notifications should be sent
#SBATCH --mail-type=END
#SBATCH --mail-user=kellet@rpi.edu

if [[ ! -d $SLURM_SUBMIT_DIR ]]
then
	mkdir -p $SLURM_SUBMIT_DIR
fi

SRCDIR=/gpfs/u/barn/GGST/GGSTkllt/trijunctionThreshold/src/ideal/large

if [[ ! -f $SRCDIR/q_GG.out ]]
then
	echo "Error: ${SRCDIR}/q_GG.out not found: cd ${SRCDIR} && make bgq"
	exit
fi

srun --runjob-opts="--mapping TEDCBA" $SRCDIR/./q_GG.out --example 2 qlarge.dat
