#!/bin/bash -l

###############################################################################
##
## NOTES:
##
## Submit as:
## 
##    sbatch ./<script_name>.sh
## 
## (1) When specifying --array as a range it must start from a positive
##     integer e.g.,
##       sbatch --array=0-9 
##     is not allowed.
##
## (2) Negative numbers are not allowed in --range
##     e.g.,
##      sbatch --array=-5,-4,-3,-2,-1,0,1,2,3,4,5
##     is NOT allowed.
##
## (3) Zero can be included if specified separately.
##    e.g., 
##       sbatch --array=0,1-9
##     is allowed.
##
## (4) Ranges can be combined with specified job numbers.
##    e.g., 
##       sbatch --array=0,1-4,6-10,50-100
##     is allowed.
##
###############################################################################
# Scheduling PARAMETERS

# Load R module:
## module load R 

# Name of the job - you'll probably want to customize this.
#SBATCH --job-name=03map
# Tell Gauss how much memory per CPU your job will use:
#SBATCH --mem-per-cpu=1120

# Array job specifications:
#SBATCH --array=1-100
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4

# Email notifications (optional), type=BEGIN, END, FAIL, ALL
# Uncomment, if desired:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cjconley@ucdavis.edu

# Standard out and Standard Error output files with the job number in the name.

#SBATCH -o /home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/03-cv-xy-yy-spacemap/logs/d%j.out
#SBATCH -e /home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/03-cv-xy-yy-spacemap/logs/d%j.err

#Project Directory
#SBATCH -D  /home/cconley/repos/sim-spacemap/hub11-20/tune/hub11-20-250/spacemap/slasso

###############################################################################
# Bash options
# End the job at the first sign of an error
set -e
#Print each command to stdout before executing it
set -v
#tell me the host of my job
hostname

#Load matlab into the name space.
#module load matlab/7.14

# Execute each of the jobs with a different index (the R script will then process
# this to do something different for each index):
/usr/bin/R --no-save --no-restore --no-site-file --no-init-file  --args ${SLURM_ARRAY_TASK_ID} 4 < 03-cv-xy-yy-spacemap.R 
