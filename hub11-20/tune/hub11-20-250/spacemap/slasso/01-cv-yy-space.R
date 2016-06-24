#Obtain data set id
args <- commandArgs(TRUE)
dataid <- as.integer(args[1])
ncores <- as.integer(args[2])

#load common settings
source("~/repos/sim-spacemap/hub11-20/tune/hub11-20-250/common-tune-params.R")
#load input data
dat <- readRDS(paste0("~/scratch-data/sim-spacemap/hub11-20/2016/datasets/n250/hub11-20-250-dataid_", sprintf("%03d", dataid), ".rds"))


#setup parallelism
suppressPackageStartupMessages(library(doParallel))
cl <- makeCluster(ncores)
registerDoParallel(cl)


#result path 
res_path <- paste0("/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/01-cv-yy-space/d", sprintf("%03d", dataid))
cmd <- paste("mkdir -p", res_path)
system(cmd)

################operate only on the Y's
#make sure the rss is on the Y--Y matrix
residConditional <- FALSE
dat$XY <- dat$Y
dat$Xindex <- dat$Yindex

#################################

#01 tuning space

lam1start <- function(n, q, alpha) { 
  sqrt(n) * qnorm(1 - (alpha/ (2*q^2)))
}
Q <- ncol(dat$Y)
lam0 <- lam1start(n = floor(N - N/kcv), q = Q, alpha = 0.005)
lam0

#neighborhood parameter about lam0. 
eps.50 <- 0.5
eps.10 <- 0.1
#initial grid size. 
ngrid <- 50
#grid
tsp <- expand.grid(slasso = sqrt(seq((lam0*(1 - eps.10))^2, (lam0*(1 + eps.50))^2, length = ngrid)))
#summary(tsp)

message("Running Space Y--Y cross-validation")
tictoc <- system.time({cvSpaceYY <- spacemap::crossValidation(data = dat, fold = kcv, 
                                                              trainIds = trainSetIds, testIds = testSetIds, 
                                                              method = "space", tuneGrid = tsp,
                                                              sridge = sridge, tol = tol, iter = iter, cd_iter = cd_iter, 
                                                              Xindex = Xindex, Yindex = Yindex, 
                                                              residConditional = residConditional, 
                                                              res_path = res_path)})

saveRDS(cvSpaceYY, file = "/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/01-cv-yy-space/01-cv_space_100_datasets.rds")

stopCluster(cl)
# 
# message(tictoc[3])
# Perf$space <- spacemap::cvPerf(cvOut = cvSpace, trueParCor = trueParCor, method = "space", 
#                                Xindex = Xindex, Yindex = Yindex, conditional = FALSE)