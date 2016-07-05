#find best models : outputs the <cv_selected> object seen below. 
source("id-best-iteration.R")
best_holdout_files <- lapply(cv_selected$dataid, get_best_holdout_files, cv_selected = cv_selected)

library(spacemap)
#each dataset has the same underlying truth
dataid <- cv_selected$dataid[1]
dat <- readRDS(paste0("~/scratch-data/sim-spacemap/hub11-20/2016/datasets/n250/hub11-20-250-dataid_", sprintf("%03d", dataid), ".rds"))
truth <- dat$trueParCor

#calculate cv-vote
library(foreach)
tol <- 1e-6
voted_spacemap_perf <- foreach(bhf = best_holdout_files, .packages = "spacemap") %do% { 
  fit <- cvVoteRDS(mod_lambda = bhf, tol = 1e-6, major_thresh = 0.5, method = "spacemap")
  formatPerf(fit = fit, truth = truth, tol = tol)
} 
library(data.table)
voted <- rbindlist(voted_spacemap_perf)
setkey(voted, comparison)
voted["(yy,xy)",]
voted[,list(median_mcc = median(mcc), median_power = median(power), median_fdr = median(fdr)),by = comparison]
voted[,list(sd_mcc = sd(mcc), sd_power = sd(power), sd_fdr = sd(fdr)),by = comparison]
saveRDS(voted, file = file.path(basedir, "spacemap_top_cv_voting_performance.rds"))
save.image(file = file.path(basedir, "spacemap_top_cv_voting.rda"))

