#!/bin/bash
#SBATCH --output dsq-extract_normalized_expression_task-%A_%2a-%N.out
#SBATCH --array 0-48
#SBATCH --job-name extract_exp
#SBATCH --mem-per-cpu=64G -t 1:00:00 -p bigmem,pi_zhao,general,scavenge

# DO NOT EDIT LINE BELOW
/ysm-gpfs/apps/software/dSQ/1.05/dSQBatch.py --job-file /gpfs/loomis/project/zhao/zy92/utmost_update/utmost_update/code/submit/extract_normalized_expression_task.txt --status-dir /gpfs/loomis/project/zhao/zy92/utmost_update/utmost_update/code/submit

