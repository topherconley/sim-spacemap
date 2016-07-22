#Obtain data set id
args <- commandArgs(TRUE)
dataid <- as.integer(args[1])
ncores <- as.integer(args[2])

#load common settings
source("~/repos/sim-spacemap/powlaw/tune/pl-mod-01-250/common-tune-params.R")
#load input data
dat <- readRDS(paste0("~/scratch-data/sim-spacemap/powlaw/pl-mod-01/datasets/n250/pl-mod-01-dataid_", sprintf("%03d", dataid), ".rds"))

#setup parallelism
suppressPackageStartupMessages(library(doParallel))
cl <- makeCluster(ncores)
registerDoParallel(cl)

#result path 
res_path <- paste0("~/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/spacemap/02-cv-xy-yy-spacemap/d", sprintf("%03d", dataid))
cmd <- paste("mkdir -p", res_path)
system(cmd)

#################################
#02 tuning spacemap
ngrid <- 10

cv01dir <- "~/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/spacemap/01-cv-yy-space"
cv01file <- file.path(cv01dir, paste0("d", sprintf("%03d", dataid), "/cv_object_dataid_", dataid,".rds"))
cv01 <- readRDS(cv01file)
slasso <- cv01$tuneGrid

tmap <- expand.grid(slasso = seq(slasso - 20, slasso + 20, length = ngrid),
                    rlasso = seq(0, 50, length = ngrid),
                    rgroup = seq(0, 50, length = ngrid))
tmap <- tmap[!(tmap$rlasso < 10 & tmap$rgroup < 10),]

saveRDS(tmap, file = file.path(res_path, paste0("tune_grid_", dataid, ".rds")))

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
