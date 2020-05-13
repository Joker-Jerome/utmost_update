#!/bin/bash
#SBATCH --output dsq-extraction_task-%A_%2a-%N.out
#SBATCH --array 0-21
#SBATCH --job-name extraction_exp
#SBATCH --mem-per-cpu=32G -t 20:00:00 -p bigmem,day

# DO NOT EDIT LINE BELOW
/ysm-gpfs/apps/software/dSQ/0.96/dSQBatch.py /gpfs/loomis/pi/zhao2/zy92/projects/utmost_update/utmost_update/code/submit/extraction_task.txt /gpfs/loomis/pi/zhao2/zy92/projects/utmost_update/utmost_update/code/submit

