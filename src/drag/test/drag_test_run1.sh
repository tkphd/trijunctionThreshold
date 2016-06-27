#!/bin/bash
# Usage: qsub -cwd -M tnk10 -pe nodal 64 drag_test_----.sh

cd /data/tnk10/trijunctionThreshold/drag/test/run1
SRCDIR=/users/tnk10/research/projects/trijunctionThreshold/src/drag/test/

MPIRUN=/usr/bin/mpirun.openmpi
SCRIPT=parallel_GG.out

rsync -a ../../../ideal/test/test.dat ./
$MPIRUN -np $NSLOTS $SRCDIR/./$SCRIPT test.dat 100000 1000
