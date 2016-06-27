#!/bin/bash
# Usage: qsub -cwd -M tnk10 -pe nodal 64 ideal_test_----.sh

cd /data/tnk10/trijunctionThreshold/ideal/test/run2
SRCDIR=/users/tnk10/research/projects/trijunctionThreshold/src/ideal/test/

MPIRUN=/usr/bin/mpirun.openmpi
SCRIPT=parallel_GG.out

rsync -a ../test.dat ./
$MPIRUN -np $NSLOTS $SRCDIR/./$SCRIPT test.dat 100000 1000
