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
res_path <- paste0("/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/02-cv-xy-yy-spacemap/d", sprintf("%03d", dataid))
cmd <- paste("mkdir -p", res_path)
system(cmd)

#################################
#02 tuning spacemap
ngrid <- 10
#the median best slasso across the 100 data sets 
med.slasso <- 101.4032
tmap <- expand.grid(slasso = seq(med.slasso - 20, med.slasso + 20, length = ngrid),
                    rlasso = seq(0, 50, length = ngrid),
                    rgroup = seq(0, 50, length = ngrid))
tmap <- tmap[!(tmap$rlasso < 10 & tmap$rgroup < 10),]

message("Running Spacemap cross-validation")
library(spacemap)
tictoc <- system.time({cvSpacemap <- spacemap::crossValidation(data = dat, fold = kcv, 
                                                              trainIds = trainSetIds, testIds = testSetIds, 
                                                              method = "spacemap", tuneGrid = tmap,
                                                              tol = tol, iter = iter, cd_iter = cd_iter, sridge = sridge, 
                                                              res_path = res_path, refitRidge = 0.0)})
cvSpacemap$tictoc <- tictoc[3]
cvSpacemap$dataid <- dataid
saveRDS(cvSpacemap, file = file.path(res_path, paste0("cv_object_dataid_", dataid, ".rds")))
stopCluster(cl)
#performance 
cvperf <- spacemap::cvPerf(cvOut = cvSpacemap, trueParCor = dat$trueParCor, method = "spacemap", tol = tol)
saveRDS(cvperf, file = file.path(res_path, paste0("cv_perf_dataid_", dataid, ".rds")))
