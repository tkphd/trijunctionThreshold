#!/bin/bash

# Run make in each directory.

# Usage:
# Serial:       ./build_all.sh
# Parallel:     ./build_all.sh parallel
# Blue Gene/Q : ./build_all.sh bgq
# Cleanup:      ./build_all.sh clean

for d in {ideal,drag}
do
	for e in {test,small,medium,large,extra}
	do
		cd $d/$e
		make $1
		cd ../..
	done
done
