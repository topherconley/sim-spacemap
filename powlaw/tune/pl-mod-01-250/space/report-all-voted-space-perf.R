
#Obtain data set id
args <- commandArgs(TRUE)

library(spacemap)
#each dataset has the same underlying truth
mod <- readRDS("/home/cconley/scratch-data/sim-spacemap/powlaw/pl-mod-01/pl-mod-01-model.rds")
truth <- list(xy = mod$model_sig_xy, yy = mod$model_sig_yy) 

#report all the tuning performance for each data set given a method
report_all_perf <- function(path, method = c("spacemap", "space", "scggm"), truth) { 
  files <- list.files(path = path, 
                      pattern = "tuneid_", full.names = TRUE)
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
  tol <- 0.0 
  library(spacemap)
  perf <- foreach(lf = lambda_files, .packages = "spacemap") %do% {
    fit <- cvVoteRDS(mod_lambda = lf, tol = tol, major_thresh = 0.5, method = method)
    formatPerf(fit = fit, truth = truth, tol = tol)
  }
  library(data.table)
  voted <- rbindlist(perf)
  setkey(voted, comparison)
  voted
}

basedir02 <- "~/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/space/01-cv-xy-yy-space/"
ddirs02 <- list.dirs(path = basedir02, recursive = F, full.names = T)
basedir03 <- "~/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/space/02-cv-xy-yy-space/"
ddirs03 <- list.dirs(path = basedir03, recursive = F, full.names = T)
ddirs <- c(ddirs02,ddirs03)
ddirs <- ddirs[!grepl(pattern = "logs", x = ddirs)]
library(foreach)
tol <- 0.0 
mapv <- foreach(path = ddirs) %do% { 
  report_all_perf(path = path, method = "space", truth = truth)
}
library(data.table)
map_voted <- rbindlist(mapv)
setkey(map_voted, comparison)
saveRDS(map_voted, file = "~/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/space/space_all_tuning_cv_vote_results.rds")
save.image(file = "~/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/space/space_all_tuning_cv_vote_results.rda")
