#!/bin/bash
#
# Usage: ./AMOSsub.sh /full/path/to/ideal_init.sh && ./AMOSsub.sh /full/path/to/ideal_run.sh && ./AMOSsub.sh /full/path/to/drag_run.sh
#
# -d singleton tells SLURM to hold the job until all previous jobs with the same name have finished.
# This allows us to queue up both the init and run scripts without worrying about the starting a job
# before the initial Voronoi tessellation is in place.

sbatch -d singleton $1
