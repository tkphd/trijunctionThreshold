#!/bin/bash
# Usage: qsub -cwd -M tnk10 -pe nodal 64 ideal_large_----.sh

cd /data/tnk10/trijunctionThreshold/ideal/large
SRCDIR=/users/tnk10/research/projects/trijunctionThreshold/src/ideal/large/

MPIRUN=/usr/bin/mpirun.openmpi
SCRIPT=parallel_GG.out

$MPIRUN -np $NSLOTS $SRCDIR/./$SCRIPT --example 2 large.dat
