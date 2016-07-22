
#Obtain data set id
args <- commandArgs(TRUE)
ncores <- as.integer(args[1])

library(foreach)
library(doParallel)
cl <- makeCluster(ncores)
registerDoParallel(cl)

#get the true model selection
dat <- readRDS(paste0("~/scratch-data/sim-spacemap/hub11-20/2016/datasets/n250/hub11-20-250-dataid_", "001", ".rds"))
truth <- dat$trueParCor

#report all the tuning performance for each data set given a method
report_all_perf <- function(path, method = c("spacemap", "space", "scggm"), truth, fpat = "tuneid_") { 
  files <- list.files(path = path, 
                      pattern = fpat, full.names = TRUE)
  tuneid <- as.integer(regmatches(files, regexpr(pattern = "(?<=tuneid_)(.*?)(?=_fold)", 
                                                 text = files, perl = TRUE)))
  
  #assures that the tuning paramters were grouped appropriately.
  lambda_files <- split(files, tuneid)
  stopifnot(all(sapply(lambda_files, 
                       function(f) length(unique(
                         substring(
                           regmatches(f, regexpr(pattern = "(?<=tuneid_)(.*?)(?=_fold)",
                                                 text = f, perl = TRUE)), first = 3, last = 10000))) == 1)))
  
  if (tuneid[1] == 0) { 
    lambda_files <- lambda_files[-1]
  }
  
  #calculate cv-vote
  library(foreach)
  tol <- 1e-6
  library(spacemap)
  perf <- foreach(lf = lambda_files, .packages = "spacemap") %do% {
    fit <- cvVoteRDS(mod_lambda = lf, tol = 1e-6, major_thresh = 0.5, method = method)
    formatPerf(fit = fit, truth = truth, tol = tol)
  }
  library(data.table)
  voted <- rbindlist(perf)
  setkey(voted, comparison)
  voted
}

basedir01 <- "/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/scggm/01-cv-xy-yy-scggm"
ddirs01 <- list.dirs(path = basedir01, recursive = F, full.names = T)
basedir02 <- "/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/scggm/02-cv-xy-yy-scggm"
ddirs02 <- list.dirs(path = basedir02, recursive = F, full.names = T)
basedir03 <- "/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/scggm/03-cv-xy-yy-scggm"
ddirs03 <- list.dirs(path = basedir03, recursive = F, full.names = T)
ddirs <- c(ddirs01, ddirs02,ddirs03)
ddirs <- ddirs[!grepl(pattern = "logs", x = ddirs)]
mapv <- foreach(path = ddirs) %dopar% { 
  report_all_perf(path = path, method = "scggm", truth = truth, fpat = "scggm_cv_tuneid_")
}
library(data.table)
map_voted <- rbindlist(mapv)
setkey(map_voted, comparison)
saveRDS(map_voted, file = "/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/scggm/scggm_all_tuning_cv_vote_results.rds")
