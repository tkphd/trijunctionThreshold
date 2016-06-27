#!/bin/bash
# Usage: qsub -cwd -M tnk10 -pe nodal 64 ideal_medium_----.sh

cd /data/tnk10/trijunctionThreshold/ideal/medium/run1
SRCDIR=/users/tnk10/research/projects/trijunctionThreshold/src/ideal/medium/

MPIRUN=/usr/bin/mpirun.openmpi
SCRIPT=parallel_GG.out

rsync -a ../medium.dat ./
$MPIRUN -np $NSLOTS $SRCDIR/./$SCRIPT medium.dat 100000 1000
