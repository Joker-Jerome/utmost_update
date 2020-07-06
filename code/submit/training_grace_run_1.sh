#!/bin/bash
#SBATCH --output dsq-training_task_test_1-%A_%4a-%N.out
#SBATCH --array 0-9807
#SBATCH --job-name utmost_training
#SBATCH --mem-per-cpu=64G -t 2:00:00 -p bigmem,day

# DO NOT EDIT LINE BELOW
/gpfs/loomis/apps/avx/software/dSQ/1.05/dSQBatch.py --job-file /gpfs/loomis/project/zhao/zy92/utmost_update/utmost_update/code/submit/training_task_test_1.txt --status-dir /gpfs/loomis/project/zhao/zy92/utmost_update/utmost_update/code/submit

