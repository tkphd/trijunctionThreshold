#!/bin/bash
# SLURM batch script for trijunction drag grain growth simulations

# Make sure your binary has been compiled before launching!
# Usage: sbatch -d singleton AMOS_ideal_small_----.sh

# Declare common SLURM settings

# Working directory and output file (datafiles  and stdout redirect)
#SBATCH -D /gpfs/u/scratch/GGST/GGSTkllt/trijunctionThreshold/ideal/small/run10
#SBATCH -o /gpfs/u/barn/GGST/GGSTkllt/trijunctionThreshold/dat/ideal/AMOS_small_run10.log

# Cluster partition and job size
#SBATCH --partition small
#SBATCH --overcommit
#SBATCH --time 1440
#SBATCH --nodes 16
#SBATCH --ntasks-per-node 32
#SBATCH -J small_2048

# When and to whom notifications should be sent
#SBATCH --mail-type=END
#SBATCH --mail-user=kellet@rpi.edu

if [[ ! -d $SLURM_SUBMIT_DIR ]]
then
	mkdir -p $SLURM_SUBMIT_DIR
fi

SRCDIR=/gpfs/u/barn/GGST/GGSTkllt/trijunctionThreshold/src/ideal/small

if [[ ! -f $SRCDIR/q_GG.out ]]
then
	echo "Error: ${SRCDIR}/q_GG.out not found: cd ${SRCDIR} && make bgq"
	exit
fi

cp ../qsmall.dat ./
srun --runjob-opts="--mapping TEDCBA" $SRCDIR/./q_GG.out qsmall.dat 100000 5000
