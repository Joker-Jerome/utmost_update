#!/bin/bash
#SBATCH --output dsq-adjustment_task-%A_%1a-%N.out
#SBATCH --array 0-3
#SBATCH --job-name adjustment
#SBATCH --mem-per-cpu=8g -t 1-00:00:00 -p bigmem,pi_zhao,general

# DO NOT EDIT LINE BELOW
/ysm-gpfs/apps/software/dSQ/0.96/dSQBatch.py /gpfs/ysm/project/zhao/zy92/utmost_update/code/submit/adjustment_task.txt /gpfs/ysm/project/zhao/zy92/utmost_update/code/submit

