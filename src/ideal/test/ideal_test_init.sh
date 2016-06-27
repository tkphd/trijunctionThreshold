#!/bin/bash
# Usage: qsub -cwd -M tnk10 -pe nodal 64 ideal_test_----.sh

cd /data/tnk10/trijunctionThreshold/ideal/test
SRCDIR=/users/tnk10/research/projects/trijunctionThreshold/src/ideal/test/

MPIRUN=/usr/bin/mpirun.openmpi
SCRIPT=parallel_GG.out

$MPIRUN -np $NSLOTS $SRCDIR/./$SCRIPT --example 2 test.dat
