#!/bin/bash -l

#SBATCH --job-name=03map
#SBATCH --mem-per-cpu=2120

# Array job specifications:
#SBATCH --array=1-100
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -o /home/cconley/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/spacemap/03-cv-xy-yy-spacemap/logs/d%j.out
#SBATCH -e /home/cconley/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/spacemap/03-cv-xy-yy-spacemap/logs/d%j.err
#SBATCH -D /home/cconley/repos/sim-spacemap/powlaw/tune/pl-mod-01-250/spacemap 
set -e
set -v
/usr/bin/R --no-save --no-restore --no-site-file --no-init-file  --args ${SLURM_ARRAY_TASK_ID} 4 < 03-cv-xy-yy-spacemap.R 
