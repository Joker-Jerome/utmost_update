#!/bin/bash
#SBATCH --output dsq-adjustment_task_12_22-%A_%2a-%N.out
#SBATCH --array 0-10
#SBATCH --job-name adjustment_farnam
#SBATCH --mem-per-cpu=32G -t 10:00:00 -p bigmem,general,pi_zhao

# DO NOT EDIT LINE BELOW
/ysm-gpfs/apps/software/dSQ/0.96/dSQBatch.py /gpfs/loomis/pi/zhao2/zy92/projects/utmost_update/utmost_update/code/submit/adjustment_task_12_22.txt /gpfs/loomis/pi/zhao2/zy92/projects/utmost_update/utmost_update/code/submit

