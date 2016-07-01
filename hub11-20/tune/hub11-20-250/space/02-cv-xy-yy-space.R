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
res_path <- paste0("/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/space/02-cv-xy-yy-space/d", sprintf("%03d", dataid))
cmd <- paste("mkdir -p", res_path)
system(cmd)

################operate only on the Y's
#make sure the rss is on the Y--Y matrix
residConditional <- TRUE

#################################

#02 tuning space
ngrid <- 500

cv01dir <- "~/scratch-data/sim-spacemap/hub11-20/2016/results/n250/space/01-cv-xy-yy-space"
cv01file <- file.path(cv01dir, paste0("d", sprintf("%03d", dataid), "/cv_object_dataid_", dataid,".rds"))
cv01 <- readRDS(cv01file)

slasso <- cv01$tuneGrid
#grid
tsp <- expand.grid(slasso = seq(slasso - 5, slasso + 5, length = ngrid))
#
message("Running Space (YY,XY) cross-validation")
library(spacemap)
tictoc <- system.time({cvSpace <- spacemap::crossValidation(data = dat, fold = kcv, 
                                                              trainIds = trainSetIds, testIds = testSetIds, 
                                                              method = "space", tuneGrid = tsp,
                                                              sridge = sridge, tol = tol, iter = iter, cd_iter = cd_iter, 
                                                              Xindex = dat$Xindex, Yindex = dat$Yindex, 
                                                              residConditional = residConditional, 
                                                              res_path = res_path, refitRidge = 0.0)})

cvSpace$tictoc <- tictoc[3]
cvSpace$dataid <- dataid
saveRDS(cvSpace, file = file.path(res_path, paste0("cv_object_dataid_", dataid, ".rds")))
stopCluster(cl)

cvperf <- spacemap::cvPerf(cvOut = cvSpace, trueParCor = dat$trueParCor, method = "space", Xindex = dat$Xindex, Yindex = dat$Yindex)
saveRDS(cvperf, file = file.path(res_path, paste0("cv_perf_dataid_", dataid, ".rds")))
