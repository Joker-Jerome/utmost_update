#!/bin/bash
#SBATCH --output dsq-task_gene_info-%A_%2a-%N.out
#SBATCH --array 0-21
#SBATCH --job-name gene_info
#SBATCH --mem-per-cpu=8G -t 1:00:00 -p bigmem,day

# DO NOT EDIT LINE BELOW
/ysm-gpfs/apps/software/dSQ/1.03/dSQBatch.py --job-file /gpfs/loomis/project/zhao/zy92/utmost_update/utmost_update/code/submit/task_gene_info.txt --status-dir /gpfs/loomis/project/zhao/zy92/utmost_update/utmost_update/code/submit

