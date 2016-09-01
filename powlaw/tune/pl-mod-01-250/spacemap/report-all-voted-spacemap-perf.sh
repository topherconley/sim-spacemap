#!/bin/bash -l

#SBATCH --job-name=map-all
#SBATCH --mem-per-cpu=1000
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -o /home/cconley/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/spacemap/report-all-tunings.out
#SBATCH -e /home/cconley/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/spacemap/report-all-tunings.err
#SBATCH -D /home/cconley/repos/sim-spacemap/powlaw/tune/pl-mod-01-250/spacemap
set -e
set -v
/usr/bin/R --no-save --no-restore --no-site-file --no-init-file  --args 31  < report-all-voted-spacemap-perf.R
