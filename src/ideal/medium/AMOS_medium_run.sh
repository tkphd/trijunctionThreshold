#!/bin/bash
# SLURM batch script for trijunction drag grain growth simulations

# Make sure your binary has been compiled before launching!
# Usage: sbatch -d singleton AMOS_ideal_medium_----.sh

# Declare common SLURM settings

# Working directory and output file (datafiles  and stdout redirect)
#SBATCH -o /gpfs/u/barn/GGST/GGSTlwsd/trijunctionThreshold/ideal/AMOS_medium_run10.log

# Cluster partition and job size
#SBATCH --partition small
#SBATCH --overcommit
#SBATCH --time 1440
#SBATCH --nodes 64
#SBATCH --ntasks-per-node 32
#SBATCH -J medium_4096

# When and to whom notifications should be sent
#SBATCH --mail-type=END
#SBATCH --mail-user=lewisd2@rpi.edu

SRCDIR=/gpfs/u/barn/GGST/GGSTlwsd/trijunctionThreshold/src/ideal/medium

if [[ ! -f $SRCD$INIDIR/q_GG.out ]]
then
	echo "Error: ${SRCDI$INIDIR/q_GG.out not found: cd ${SRCDIR} && make bgq"
	exit
fi

DATDIR=/gpfs/u/scratch/GGST/GGSTlwsd/trijunctionThreshold/ideal/medium/run10
if [[ ! -d $DATDIR ]]
then
	mkdir -p $DATDIR
fi

INIDIR=/gpfs/u/scratch/GGST/GGSTlwsd/trijunctionThreshold/ideal/medium

cp $INIDIR/qmedium.dat $DATDIR/
srun -D $DATDIR --runjob-opts="--mapping TEDCBA" $SRCDIR/./q_GG.out qmedium.dat 100000 5000
