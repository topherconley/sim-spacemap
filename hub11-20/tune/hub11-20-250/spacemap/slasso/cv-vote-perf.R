
basedir <- "/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/03-cv-xy-yy-spacemap/"
datdir <- "~/scratch-data/sim-spacemap/hub11-20/2016/datasets/n250/"
library(spacemap)
spacemap_perf <- foreach(dataid = 1:100) %do% { 
  ddir <- file.path(basedir, paste0("d", sprintf("%03d", dataid)))
  datfile <- file.path(datdir, paste0("hub11-20-250-dataid_", sprintf("%03d", dataid), ".rds"))
  dat <- readRDS(datfile)
  cvint <- readRDS(file = paste0(ddir, "/cv_object_dataid_", dataid, ".rds"))
  tid <- cvint$minIndices["rss"]                                 
  tf <- list.files(path = ddir, pattern = paste0("tuneid_", sprintf("%03d", tid)))
  cvrds <- cvVoteRDS(mod_lambda = tf, tol = 1e-6, major_thresh= 0.5, method = "spacemap")
  saveRDS(cvrds, file = paste0(ddir, "/min_rss_cv_vote_dataid_", dataid, ".rds"))
  formatPerf(fit = cvrds, truth = dat$trueParCor, method = "Spacemap", tol = 1e-06)
}

library(data.table)
spacemap_perf_dt <- rbindlist(spacemap_perf)
spacemap_perf_dt[,list(median_mcc = median(mcc))]
spacemap_perf_dt[,list(median_power = median(mcc), median_fdr = median(fdr)),by = comparison]

