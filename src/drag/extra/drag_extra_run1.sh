#!/bin/bash
# Usage: qsub -cwd -M tnk10 -pe nodal 64 drag_extra_----.sh

cd /data/tnk10/trijunctionThreshold/drag/extra/run1
SRCDIR=/users/tnk10/research/projects/trijunctionThreshold/src/drag/extra/

MPIRUN=/usr/bin/mpirun.openmpi
SCRIPT=parallel_GG.out

rsync -a ../../../ideal/extra/extra.dat ./
$MPIRUN -np $NSLOTS $SRCDIR/./$SCRIPT extra.dat 100000 1000
