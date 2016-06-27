#!/bin/bash
# Usage: qsub -cwd -M tnk10 -pe nodal 64 drag_large_----.sh

cd /data/tnk10/trijunctionThreshold/drag/large/run1
SRCDIR=/users/tnk10/research/projects/trijunctionThreshold/src/drag/large/

MPIRUN=/usr/bin/mpirun.openmpi
SCRIPT=parallel_GG.out

rsync -a ../../../ideal/large/large.dat ./
$MPIRUN -np $NSLOTS $SRCDIR/./$SCRIPT large.dat 100000 1000
