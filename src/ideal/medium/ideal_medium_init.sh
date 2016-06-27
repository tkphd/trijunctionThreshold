#!/bin/bash
# Usage: qsub -cwd -M tnk10 -pe nodal 64 ideal_medium_----.sh

cd /data/tnk10/trijunctionThreshold/ideal/medium
SRCDIR=/users/tnk10/research/projects/trijunctionThreshold/src/ideal/medium/

MPIRUN=/usr/bin/mpirun.openmpi
SCRIPT=parallel_GG.out

$MPIRUN -np $NSLOTS $SRCDIR/./$SCRIPT --example 2 medium.dat
