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

#  X Hub  Power and FDR calculation 
#true positive  = X  with at least one true edge and allows for some false positives
#false positive = X with all false positive edges
#fale negative = X with no predicted edges but indeed is a X hub

simple_xhub_spacemap_perf <- foreach(bhf = best_holdout_files, .packages = "spacemap") %do% { 
  bhf <- sub(pattern = "/home/cconley/", replacement = "~/", x = bhf)
  fit <- cvVoteRDS(mod_lambda = bhf, tol = 1e-6, major_thresh = 0.5, method = "spacemap")
  xh  <- which(rowSums(abs(truth$xy) > tol) > 0)
  nxh <- which(rowSums(abs(truth$xy) > tol) == 0)
  pxh <- which(rowSums(abs(fit$xy) > tol) > 0)
  pnxh <- which(rowSums(abs(fit$xy) > tol) == 0)
  tp <- length(intersect(pxh, xh))
  fp <- length(setdiff(pxh, xh))
  fn <- length(setdiff(xh, pxh))
  tn <- length(intersect(pnxh, nxh))
  stopifnot(sum(c(tp,fp,fn,tn)) == nrow(fit$xy))
  data.frame(power = algoPower(tp, fn), fdr = algoFDR(tp, fp), mcc = algoMCC(tp, fn, fp, tn))
} 

library(data.table)
simple_xhub_voted <- rbindlist(simple_xhub_spacemap_perf)
simple_xhub_voted[,list(mean_mcc = mean(mcc), mean_power = mean(power), mean_fdr = mean(fdr))]
basedir <- sub("/home/cconley/", "~/", basedir)
saveRDS(simple_xhub_voted, file = file.path(basedir, "simple_xhub_spacemap_top_cv_voting_performance.rds"))
save.image(file = file.path(basedir, "simple_xhub_spacemap_top_cv_voting.rda"))

