#!/bin/bash -l
###############################################################################
# Scheduling PARAMETERS

#SBATCH --job-name=scggm
#SBATCH --mem-per-cpu=1120
#SBATCH --array=51-100
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH -o /home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/scggm/01-cv-xy-yy-scggm/logs/%j.out
#SBATCH -e /home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/scggm/01-cv-xy-yy-scggm/logs/%j.err
#SBATCH -D /home/cconley/repos/sim-spacemap/hub11-20/tune/hub11-20-250/scggm/
# End the job at the first sign of an error
set -e
#Print each command to stdout before executing it
set -v
#tell me the host of my job
hostname
module add matlab
/usr/bin/R --no-save --no-restore --no-site-file --no-init-file  --args ${SLURM_ARRAY_TASK_ID} 31 < 01-cv-xy-yy-scggm.R 
