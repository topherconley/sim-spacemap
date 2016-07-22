#!/bin/bash -l

# Scheduling PARAMETERS
#SBATCH --job-name=hsp01yy
#SBATCH --mem-per-cpu=4120
#SBATCH --array=1-100
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -o /home/cconley/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/spacemap/01-cv-yy-space/logs/d%j.out
#SBATCH -e /home/cconley/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/spacemap/01-cv-yy-space/logs/d%j.err
#SBATCH -D /home/cconley/repos/sim-spacemap/powlaw/tune/pl-mod-01-250/spacemap 
set -e
set -v
/usr/bin/R --no-save --no-restore --no-site-file --no-init-file  --args ${SLURM_ARRAY_TASK_ID} 8 < 01-cv-yy-space.R 
