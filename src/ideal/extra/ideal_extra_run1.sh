#!/bin/bash
# Usage: qsub -cwd -M tnk10 -pe nodal 64 ideal_extra_----.sh

cd /data/tnk10/trijunctionThreshold/ideal/extra/run1
SRCDIR=/users/tnk10/research/projects/trijunctionThreshold/src/ideal/extra/

MPIRUN=/usr/bin/mpirun.openmpi
SCRIPT=parallel_GG.out

rsync -a ../extra.dat ./
$MPIRUN -np $NSLOTS $SRCDIR/./$SCRIPT extra.dat 100000 1000
