#!/bin/bash
# Usage: qsub -cwd -M tnk10 -pe nodal 64 drag_medium_----.sh

cd /data/tnk10/trijunctionThreshold/drag/medium/run1
SRCDIR=/users/tnk10/research/projects/trijunctionThreshold/src/drag/medium/

MPIRUN=/usr/bin/mpirun.openmpi
SCRIPT=parallel_GG.out

rsync -a ../../../ideal/medium/medium.dat ./
$MPIRUN -np $NSLOTS $SRCDIR/./$SCRIPT medium.dat 100000 1000
