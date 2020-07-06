#!/bin/bash
#SBATCH --output dsq-training_task_test-%A_%1a-%N.out
#SBATCH --array 0-4
#SBATCH --job-name ctimp
#SBATCH --mem-per-cpu=64G -t 1:00:00 -p bigmem,pi_zhao

# DO NOT EDIT LINE BELOW
/ysm-gpfs/apps/software/dSQ/1.04/dSQBatch.py --job-file /gpfs/loomis/project/zhao/zy92/utmost_update/utmost_update/code/submit/training_task_test.txt --status-dir /gpfs/loomis/project/zhao/zy92/utmost_update/utmost_update/code/submit

