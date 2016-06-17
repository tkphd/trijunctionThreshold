#!/bin/bash
# Usage: ./qsub.sh PREREQ_ID script_to_launch.sh

PRIOR=$1
SCRIPT=$2
CORES=64

qsub -cwd -m be -M tnk10 -pe nodal $CORES -q rack4,rack5 -hold_jid $PRIOR $SCRIPT
