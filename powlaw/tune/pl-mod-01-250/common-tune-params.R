
#Common settings to tuning

#number of samples
N <- 250
#number of folds
kcv <- 10
#tolerance 
tol <- 1e-6
#number of iterations for estimating sigs
iter <- 3
## max coordinate descent iterations 
cd_iter <- 1e7
## ridge penalty for elastic net of space
sridge <- 0 

#partition
foldfile <- "/home/cconley/scratch-data/sim-spacemap/powlaw/pl-mod-01/datasets/n250/pl-mod-01-n250-train-test-partition.rds"
partition <- readRDS(file = foldfile)
testSetIds <- partition$testSetIds 
trainSetIds <- partition$trainSetIds
