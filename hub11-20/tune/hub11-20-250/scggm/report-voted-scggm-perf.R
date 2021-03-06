
basedir <- "~/scratch-data/sim-spacemap/hub11-20/2016/results/n250/scggm"

#find best models : outputs the <cv_selected> object seen below. 
source("id-best-iteration.R")
#indexed by dataid
best_holdout_files <- lapply(cv_selected$dataid, get_best_holdout_files, cv_selected = cv_selected)
top_cv_files  <- lapply(cv_selected$dataid, set_top_cv_vote_path, cv_selected = cv_selected)

#each dataset has the same underlying truth
dataid <- cv_selected$dataid[1]
dat <- readRDS(paste0("~/scratch-data/sim-spacemap/hub11-20/2016/datasets/n250/hub11-20-250-dataid_", sprintf("%03d", dataid), ".rds"))
truth <- dat$trueParCor

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
voted[,list(median_mcc = median(mcc), median_power = median(power), median_fdr = median(fdr)),by = comparison]
voted[,list(sd_mcc = sd(mcc), sd_power = sd(power), sd_fdr = sd(fdr)),by = comparison]
saveRDS(voted, file = file.path(basedir, "scggm_top_cv_voting_performance.rds"))
save.image(file = file.path(basedir, "scggm_top_cv_voting.rda"))

#  XY  Power and FDR calculation with just X hubs (ignoring non-x hubs)

xhub_voted_scggm_perf <- foreach(bhf = best_holdout_files, .packages = "spacemap") %do% { 
  bhf <- sub(pattern = "/home/cconley/", replacement = "~/", x = bhf)
  fit <- cvVoteRDS(mod_lambda = bhf, tol = 1e-6, major_thresh = 0.5, method = "scggm")
  xhub_index <- which(rowSums(abs(truth$xy) > tol) > 0)
  xhub_fit <- fit
  xhub_fit$xy <- fit$xy[xhub_index,]
  xhub_truth <- truth
  xhub_truth$xy <- xhub_truth$xy[xhub_index,]
  formatPerf(fit = xhub_fit, truth = xhub_truth, tol = tol)[3,]
} 

xhub_voted <- rbindlist(xhub_voted_scggm_perf)
xhub_voted[,list(mean_mcc = mean(mcc), mean_power = mean(power), mean_fdr = mean(fdr))]
xhub_voted[,list(sd_mcc = sd(mcc), sd_power = sd(power), sd_fdr = sd(fdr)),by = comparison]
basedir <- sub("/home/cconley/", "~/", basedir)
saveRDS(xhub_voted, file = file.path(basedir, "xhub_scggm_top_cv_voting_performance.rds"))
save.image(file = file.path(basedir, "xhub_scggm_top_cv_voting.rda"))


#  X Hub  Power and FDR calculation 
#true positive  = X  with at least one true edge and allows for some false positives
#false positive = X with all false positive edges
#fale negative = X with no predicted edges but indeed is a X hub

simple_xhub_scggm_perf <- foreach(bhf = best_holdout_files, .packages = "spacemap") %do% { 
  bhf <- sub(pattern = "/home/cconley/", replacement = "~/", x = bhf)
  fit <- cvVoteRDS(mod_lambda = bhf, tol = 1e-6, major_thresh = 0.5, method = "scggm")
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
simple_xhub_voted <- rbindlist(simple_xhub_scggm_perf)
simple_xhub_voted[,list(mean_mcc = mean(mcc), mean_power = mean(power), mean_fdr = mean(fdr))]
basedir <- sub("/home/cconley/", "~/", basedir)
saveRDS(simple_xhub_voted, file = file.path(basedir, "simple_xhub_scggm_top_cv_voting_performance.rds"))
save.image(file = file.path(basedir, "simple_xhub_scggm_top_cv_voting.rda"))

