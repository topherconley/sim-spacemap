#find best models : outputs the <cv_selected> object seen below. 
source("id-best-iteration.R")
best_holdout_files <- lapply(cv_selected$dataid, get_best_holdout_files, cv_selected = cv_selected)

library(spacemap)
#each dataset has the same underlying truth
dataid <- cv_selected$dataid[1]
mod <- readRDS("/home/cconley/scratch-data/sim-spacemap/powlaw/pl-mod-01/pl-mod-01-model.rds")
truth <- list(xy = mod$model_sig_xy, yy = mod$model_sig_yy) 

#calculate cv-vote
library(foreach)
tol <- 0.0
voted_spacemap_perf <- foreach(bhf = best_holdout_files, .packages = "spacemap") %do% { 
  fit <- cvVoteRDS(mod_lambda = bhf, tol = tol, major_thresh = 0.5, method = "spacemap")
  formatPerf(fit = fit, truth = truth, tol = tol)
} 
library(data.table)
voted <- rbindlist(voted_spacemap_perf)
setkey(voted, comparison)
voted["(yy,xy)",]
voted[,list(mean_mcc = mean(mcc), mean_power = mean(power), mean_fdr = mean(fdr)),by = comparison]
voted[,list(sd_mcc = sd(mcc), sd_power = sd(power), sd_fdr = sd(fdr)),by = comparison]
saveRDS(voted, file = file.path(basedir, "spacemap_top_cv_voting_performance.rds"))
save.image(file = file.path(basedir, "spacemap_top_cv_voting.rda"))

