#!/bin/bash
# Usage: ./qsub.sh script_to_launch.sh

SCRIPT=$1
CORES=64

qsub -cwd -m be -M tnk10 -pe nodal $CORES -q rack4,rack5 $SCRIPT
