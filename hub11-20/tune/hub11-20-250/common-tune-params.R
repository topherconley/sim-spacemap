
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

library(caret)
set.seed(95482)
testSetIds <- caret::createFolds(1:N, k = kcv)
trainSetIds <- lapply(testSetIds, function(ts) setdiff(1:N, ts))
#assure that the test set and train sets are non-overlapping
stopifnot(!any(sapply(1:10, function(i) length(intersect(trainSetIds[[i]], testSetIds[[i]]))) > 0))

#foldfile <- "/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/train-test-partition.rds"
#saveRDS(list(test = testSetIds, train = trainSetIds), file = foldfile)
