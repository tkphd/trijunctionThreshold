#!/bin/bash
# SLURM batch script for trijunction drag grain growth simulations

# Make sure your binary has been compiled before launching!
# Usage: sbatch -d singleton AMOS_drag_small_----.sh

# Declare common SLURM settings

# Working directory and output file (datafiles  and stdout redirect)
#SBATCH -o /gpfs/u/barn/GGST/GGSTlwsd/trijunctionThreshold/dat/drag/AMOS_small_run10.log

# Cluster partition and job size
#SBATCH --partition small
#SBATCH --overcommit
#SBATCH --time 1440
#SBATCH --nodes 16
#SBATCH --ntasks-per-node 32
#SBATCH -J small_2048

# When and to whom notifications should be sent
#SBATCH --mail-type=END
#SBATCH --mail-user=lewisd2@rpi.edu

SRCDIR=/gpfs/u/barn/GGST/GGSTlwsd/trijunctionThreshold/src/drag/small
if [[ ! -f $SRCDIR/q_GG.out ]]
then
	echo "Error: ${SRCDIR}/q_GG.out not found: cd ${SRCDIR} && make bgq"
	exit
fi

DATDIR=/gpfs/u/scratch/GGST/GGSTlwsd/trijunctionThreshold/drag/small/run10
if [[ ! -d $DATDIR ]]
then
	mkdir -p $DATDIR
fi

INIDIR=/gpfs/u/scratch/GGST/GGSTlwsd/trijunctionThreshold/ideal/small

cp $INIDIR/qsmall.dat $DATDIR/
srun --runjob-opts="--mapping TEDCBA" $SRCDIR/./q_GG.out $DATDIR/qsmall.dat 100000 5000
