#!/bin/bash
#SBATCH --output dsq-summary_task-%A_%1a-%N.out
#SBATCH --array 0-8
#SBATCH --job-name snp_summary
#SBATCH --mem-per-cpu=16G -t 10:00:00 -p general,pi_zhao,bigmem

# DO NOT EDIT LINE BELOW
/ysm-gpfs/apps/software/dSQ/0.96/dSQBatch.py /gpfs/loomis/pi/zhao2/zy92/projects/utmost_update/utmost_update/code/submit/summary_task.txt /gpfs/loomis/pi/zhao2/zy92/projects/utmost_update/utmost_update/code/submit

