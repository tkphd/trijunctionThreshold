#!/bin/bash
# SLURM batch script for trijunction drag grain growth simulations

# Make sure your binary has been compiled before launching!
# Usage: sbatch -d singleton AMOS_ideal_extra_----.sh

# Declare common SLURM settings

# Working directory and output file (datafiles  and stdout redirect)
#SBATCH -o /gpfs/u/barn/GGST/GGSTlwsd/trijunctionThreshold/dat/ideal/AMOS_extra_run10.log

# Cluster partition and job size
#SBATCH --partition medium
#SBATCH --overcommit
#SBATCH --time 720
#SBATCH --nodes 512
#SBATCH --ntasks-per-node 32
#SBATCH -J extra_16384

# When and to whom notifications should be sent
#SBATCH --mail-type=END
#SBATCH --mail-user=lewisd2@rpi.edu

SRCDIR=/gpfs/u/barn/GGST/GGSTlwsd/trijunctionThreshold/src/ideal/extra
if [[ ! -f $SRCDIR/q_GG.out ]]
then
	echo "Error: $SRCDIR/q_GG.out not found: cd ${SRCDIR} && make bgq"
	exit
fi

DATDIR=/gpfs/u/scratch/GGST/GGSTlwsd/trijunctionThreshold/ideal/extra/run10
if [[ ! -d $DATDIR ]]
then
	mkdir -p $DATDIR
fi

INIDIR=/gpfs/u/scratch/GGST/GGSTlwsd/trijunctionThreshold/ideal/extra

cp $INIDIR/qextra.dat $DATDIR/
srun --runjob-opts="--mapping TEDCBA" $SRCDIR/./q_GG.out $DATDIR/qextra.dat 100000 5000
