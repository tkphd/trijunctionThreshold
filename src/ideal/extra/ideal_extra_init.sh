#!/bin/bash
# Usage: qsub -cwd -M tnk10 -pe nodal 64 ideal_extra_----.sh

cd /data/tnk10/trijunctionThreshold/ideal/extra
SRCDIR=/users/tnk10/research/projects/trijunctionThreshold/src/ideal/extra/

MPIRUN=/usr/bin/mpirun.openmpi
SCRIPT=parallel_GG.out

$MPIRUN -np $NSLOTS $SRCDIR/./$SCRIPT --example 2 extra.dat
