#!/bin/bash
# Usage: qsub -cwd -M tnk10 -pe nodal 64 ideal_large_----.sh

cd /data/tnk10/trijunctionThreshold/ideal/large/run1
SRCDIR=/users/tnk10/research/projects/trijunctionThreshold/src/ideal/large/

MPIRUN=/usr/bin/mpirun.openmpi
SCRIPT=parallel_GG.out

rsync -a ../large.dat ./
$MPIRUN -np $NSLOTS $SRCDIR/./$SCRIPT large.dat 100000 1000
