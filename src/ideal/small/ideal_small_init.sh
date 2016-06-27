#!/bin/bash
# Usage: qsub -cwd -M tnk10 -pe nodal 64 ideal_small_----.sh

cd /data/tnk10/trijunctionThreshold/ideal/small
SRCDIR=/users/tnk10/research/projects/trijunctionThreshold/src/ideal/small/

MPIRUN=/usr/bin/mpirun.openmpi
SCRIPT=parallel_GG.out

$MPIRUN -np $NSLOTS $SRCDIR/./$SCRIPT --example 2 small.dat
