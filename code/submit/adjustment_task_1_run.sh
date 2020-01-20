#!/bin/bash
#SBATCH --output dsq-adjustment_task_1-%A_%1a-%N.out
#SBATCH --array 0
#SBATCH --job-name adjustment_1
#SBATCH --mem-per-cpu=32G -t 10:00:00 -p bigmem,week,day

# DO NOT EDIT LINE BELOW
/gpfs/loomis/apps/avx/software/dSQ/0.96/dSQBatch.py /gpfs/loomis/pi/zhao2/zy92/projects/utmost_update/utmost_update/code/submit/adjustment_task_1.txt /gpfs/loomis/pi/zhao2/zy92/projects/utmost_update/utmost_update/code/submit

