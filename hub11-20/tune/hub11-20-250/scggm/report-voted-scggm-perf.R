
basedir <- "/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/scggm"

#find best models
source("id-best-iteration.R")
#
source("cv_vote.R")
source("scggm_perf.R")

#find best hold outs per data id
get_best_holdout_files <- function(i, cv_selected) { 
  ddir <- file.path(basedir,  
                    paste0(sprintf("%02d", cv_selected$top_iter[i]), "-cv-xy-yy-scggm"), 
                    paste0("d", sprintf("%03d", cv_selected$dataid[i])))
  bestho <- list.files(path = ddir, pattern = paste0("scggm_cv_tuneid_", 
                                                     sprintf("%03d",cv_selected$top_tuneid[i])), 
                       full.names = TRUE)
}
#indexed by dataid
best_holdout_files <- lapply(cv_selected$dataid, get_best_holdout_files, cv_selected = cv_selected)


#each dataset has the same underlying truth
dataid <- cv_selected$dataid[1]
dat <- readRDS(paste0("~/scratch-data/sim-spacemap/hub11-20/2016/datasets/n250/hub11-20-250-dataid_", sprintf("%03d", dataid), ".rds"))
truth <- dat$trueParCor

#calculate cv-vote
library(foreach)
tol <- 1e-6
voted_scggm_perf <- foreach(bhf = best_holdout_files) %do% { 
  fit <- cvVote(mod_lambda = bhf, tol = 1e-6, major_thresh = 0.5, method = "scggm")
  scggm_perf(fit = fit, truth = truth, tol = tol)
} 
