
#Common settings to tuning

#number of samples
N <- 250
#number of folds
kcv <- 10
#tolerance 
tol <- 1e-3
#number of iterations for estimating sigs
iter <- 3
## max coordinate descent iterations 
cd_iter <- 1e7
## ridge penalty for elastic net of space
sridge <- 0 

library(caret)
set.seed(95482)
testSetIds <- caret::createFolds(1:N, k = kcv)
trainSetids <- lapply(testSetIds, function(ts) setdiff(1:N, ts))
#assure that the test set and train sets are non-overlapping
stopifnot(!any(sapply(1:10, function(i) length(intersect(trainSetids[[i]], testSetIds[[i]]))) > 0))
