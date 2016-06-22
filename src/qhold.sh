#!/bin/bash
# Usage: ./qsub.sh PREREQ_JOB_ID /full/path/to/script_to_launch.sh

PRIOR=$1
SCRIPT=$2
CORES=64

findIdeal=$(echo $SCRIPT | grep -q ideal)
findDrag=$(echo $SCRIPT | grep -q drag)

cd ../dat

if [[ $findIdeal ]]
then
	cd ideal
elif [[ $findDrag ]]
then
	cd drag
fi

qsub -cwd -m be -M tnk10 -pe nodal $CORES -q rack4,rack5 -hold_jid $PRIOR $SCRIPT
