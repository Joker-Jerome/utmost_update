#!/bin/bash
#SBATCH --output dsq-adjustment_task_1_22-%A_%2a-%N.out
#SBATCH --array 0-21
#SBATCH --job-name adjustment_task
#SBATCH --mem-per-cpu=16G -t 10:00:00 -p bigmem,day

# DO NOT EDIT LINE BELOW
/ysm-gpfs/apps/software/dSQ/1.0/dSQBatch.py --job-file /gpfs/loomis/pi/zhao2/zy92/projects/utmost_update/utmost_update/code/submit/adjustment_task_1_22.txt --status-dir /gpfs/loomis/pi/zhao2/zy92/projects/utmost_update/utmost_update/code/submit

