
basedir <- "~/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/"

source("~/repos/sim-spacemap/hub11-20/tune/hub11-20-250/common-tune-params.R")
dat <- readRDS(paste0("~/scratch-data/sim-spacemap/hub11-20/2016/datasets/n250/hub11-20-250-dataid_", sprintf("%03d", 1), ".rds"))
truth <- dat$trueParCor
library(foreach)
library(spacemap)
dataids <- seq_len(6)
lbp <- foreach(dataid = dataids) %do% { 
  res_path <- paste0("~/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/boot-spacemap/d", sprintf("%03d", dataid))
  bvSpmap <- readRDS(file = file.path(res_path, paste0("ninety_perc_boot_vote_dataid_", sprintf("%03d", dataid), ".rds")))
  fit <- list(yy = bvSpmap$bv$ParCor, xy = bvSpmap$bv$Gamma)
  formatPerf(fit = fit, truth = truth, tol = tol)
}
library(data.table)
perf_boot_vote <- rbindlist(lbp)
saveRDS(perf_boot_vote, 
        file = file.path(basedir, "spacemap_ninety_perc_boot_vote_performance.rds"))

xbarv <- perf_boot_vote[,list(mean_mcc = mean(mcc), mean_power = mean(power), mean_fdr = mean(fdr)),by = comparison]
sdv <- perf_boot_vote[,list(sd_mcc = sd(mcc), sd_power = sd(power), sd_fdr = sd(fdr)),by = comparison]
power <- c(xbarv$mean_mcc[1], xbarv$mean_power[2:3], xbarv$mean_fdr[2:3])
sd <- c(sdv$sd_mcc[1], sdv$sd_power[2:3], sdv$sd_fdr[2:3])
fmt <- paste0(signif(power, 4), " (", signif(sd, 2), ")")

#  XY  Power and FDR calculation with just X hubs (ignoring non-x hubs)

xhub_voted_spacemap_perf <- foreach(dataid = dataids) %do% { 
  res_path <- paste0("~/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/boot-spacemap/d", sprintf("%03d", dataid))
  bvSpmap <- readRDS(file = file.path(res_path, paste0("ninety_perc_boot_vote_dataid_", sprintf("%03d", dataid), ".rds")))
  fit <- list(yy = bvSpmap$bv$ParCor, xy = bvSpmap$bv$Gamma)
  xhub_index <- which(rowSums(abs(truth$xy) > tol) > 0)
  xhub_fit <- fit
  xhub_fit$xy <- fit$xy[xhub_index,]
  xhub_truth <- truth
  xhub_truth$xy <- xhub_truth$xy[xhub_index,]
  formatPerf(fit = xhub_fit, truth = xhub_truth, tol = tol)[3,]
} 

xhub_voted <- rbindlist(xhub_voted_spacemap_perf)
xhub_voted[,list(mean_mcc = mean(mcc), mean_power = mean(power), mean_fdr = mean(fdr))]
xhub_voted[,list(sd_mcc = sd(mcc), sd_power = sd(power), sd_fdr = sd(fdr)),by = comparison]
saveRDS(xhub_voted, file = file.path(basedir, "xhub_ninety_perc_spacemap_top_cv_voting_performance.rds"))
save.image(file = file.path(basedir, "xhub_ninety_perc_spacemap_top_cv_voting.rda"))


#  X Hub  Power and FDR calculation 
#true positive  = X  with at least one true edge and allows for some false positives
#false positive = X with all false positive edges
#fale negative = X with no predicted edges but indeed is a X hub

simple_xhub_voted_spacemap_perf <- foreach(dataid = dataids) %do% { 
  res_path <- paste0("~/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/boot-spacemap/d", sprintf("%03d", dataid))
  bvSpmap <- readRDS(file = file.path(res_path, paste0("ninety_perc_boot_vote_dataid_", sprintf("%03d", dataid), ".rds")))
  fit <- list(yy = bvSpmap$bv$ParCor, xy = bvSpmap$bv$Gamma)
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
simple_xhub_voted <- rbindlist(simple_xhub_voted_spacemap_perf)
simple_xhub_voted[,list(mean_mcc = mean(mcc), mean_power = mean(power), mean_fdr = mean(fdr))]
saveRDS(simple_xhub_voted, file = file.path(basedir, "simple_xhub_ninety_perc_spacemap_top_cv_voting_performance.rds"))
save.image(file = file.path(basedir, "simple_xhub_ninety_perc_spacemap_top_cv_voting.rda"))

