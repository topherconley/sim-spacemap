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
res_path <- paste0("/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/03-cv-xy-yy-spacemap/d", sprintf("%03d", dataid))
cmd <- paste("mkdir -p", res_path)
system(cmd)

#################################
#03 tuning spacemap
ngrid <- 10

cv02dir <- "~/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/02-cv-xy-yy-spacemap"
cv02file <- file.path(cv02dir, paste0("d", sprintf("%03d", dataid), "/cv_object_dataid_", dataid,".rds"))
cv02 <- readRDS(cv02file)

slasso <- cv02$tuneGrid$slasso
rlasso <- cv02$tuneGrid$rlasso
rgroup <- cv02$tuneGrid$rgroup

tmap <- expand.grid(slasso = seq(max(slasso - 5,0), slasso + 5, length = ngrid),
                    rlasso = seq(max(rlasso - 5,0), rlasso + 5, length = ngrid),
                    rgroup = seq(max(rgroup - 10,0), rgroup + 10, length = ngrid))

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
