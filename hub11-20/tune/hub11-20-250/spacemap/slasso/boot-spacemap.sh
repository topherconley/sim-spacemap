#!/bin/bash -l
#SBATCH --job-name=boot
#SBATCH --mem-per-cpu=2120
#SBATCH --array=1-100
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH -o  /home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/boot-spacemap/logs/d%j.out
#SBATCH -e  /home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/boot-spacemap/logs/d%j.err
#SBATCH -D  /home/cconley/repos/sim-spacemap/hub11-20/tune/hub11-20-250/spacemap/slasso
set -e
set -v

/usr/bin/R --no-save --no-restore --no-site-file --no-init-file  --args ${SLURM_ARRAY_TASK_ID} 4 < boot-spacemap.R 
