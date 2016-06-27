#!/bin/bash
# Usage: qsub -cwd -M tnk10 -pe nodal 64 drag_small_----.sh

cd /data/tnk10/trijunctionThreshold/drag/small/run1
SRCDIR=/users/tnk10/research/projects/trijunctionThreshold/src/drag/small/

MPIRUN=/usr/bin/mpirun.openmpi
SCRIPT=parallel_GG.out

rsync -a ../../../ideal/small/small.dat ./
$MPIRUN -np $NSLOTS $SRCDIR/./$SCRIPT small.dat 100000 1000
