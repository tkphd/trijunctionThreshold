#!/bin/bash
# Usage: ./qsub.sh /full/path/to/script_to_launch.sh

SCRIPT=$1
CORES=64

cd /users/tnk10/research/projects/trijunctionThreshold/dat

qsub -cwd -m be -M tnk10 -pe nodal $CORES -q rack4,rack5 $SCRIPT
