#!/bin/bash -l

#SBATCH --job-name=-all-tune
#SBATCH --mem-per-cpu=1000
#SBATCH --nodes=1
#SBATCH --cpus-per-task=31
#SBATCH -o /home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/space/report-all-tunings.out
#SBATCH -e /home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/space/report-all-tunings.err
#SBATCH -D  /home/cconley/repos/sim-spacemap/hub11-20/tune/hub11-20-250/space
set -e
set -v
/usr/bin/R --no-save --no-restore --no-site-file --no-init-file  --args 31  < report-all-voted-space-perf.R
