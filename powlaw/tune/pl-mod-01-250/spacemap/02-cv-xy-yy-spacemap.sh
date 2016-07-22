#!/bin/bash -l

#SBATCH --job-name=02map
#SBATCH --mem-per-cpu=1120

# Array job specifications:
#SBATCH --array=1-100
#SBATCH --nodes=1
#SBATCH --cpus-per-task=31
#SBATCH -o /home/cconley/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/spacemap/02-cv-xy-yy-spacemap/logs/d%j.out
#SBATCH -e /home/cconley/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/spacemap/02-cv-xy-yy-spacemap/logs/d%j.err
#SBATCH -D /home/cconley/repos/sim-spacemap/powlaw/tune/pl-mod-01-250/spacemap 
set -e
set -v
/usr/bin/R --no-save --no-restore --no-site-file --no-init-file  --args ${SLURM_ARRAY_TASK_ID} 31 < 02-cv-xy-yy-spacemap.R 
