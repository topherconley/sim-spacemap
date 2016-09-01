
basedir <- "~/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/scggm/"

#find best models : outputs the <cv_selected> object seen below. 
source("id-best-iteration.R")
#indexed by dataid
best_holdout_files <- lapply(seq_along(did), get_best_holdout_files, cv_selected = cv_selected)
top_cv_files  <- lapply(seq_along(did), set_top_cv_vote_path, cv_selected = cv_selected)

library(spacemap)
#each dataset has the same underlying truth
dataid <- cv_selected$dataid[1]
mod <- readRDS("/home/cconley/scratch-data/sim-spacemap/powlaw/pl-mod-01/pl-mod-01-model.rds")
truth <- list(xy = mod$model_sig_xy, yy = mod$model_sig_yy) 

#calculate cv-vote
library(foreach)
library(spacemap)
tol <- 1e-6
voted_scggm_perf <- foreach(bhf = best_holdout_files, tcvp = top_cv_files, .packages = "spacemap") %do% { 
  fit <- cvVoteRDS(mod_lambda = bhf, tol = tol, major_thresh = 0.5, method = "scggm")
  saveRDS(fit, file = tcvp)
  formatPerf(fit = fit, truth = truth, tol = tol)
} 
library(data.table)
voted <- rbindlist(voted_scggm_perf)
setkey(voted, comparison)
voted["(yy,xy)",]
voted[,list(mean_mcc = mean(mcc), mean_power = mean(power), mean_fdr = mean(fdr)),by = comparison]
voted[,list(sd_mcc = sd(mcc), sd_power = sd(power), sd_fdr = sd(fdr)),by = comparison]
saveRDS(voted, file = file.path(basedir, "scggm_top_cv_voting_performance.rds"))
save.image(file = file.path(basedir, "scggm_top_cv_voting.rda"))

