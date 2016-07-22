#####################
dataids <- 1:100
report_cv_score <- function(basedir, testSetLen, dataids = 1:100, tol = 1e-6, vote = TRUE) {
  library(foreach)
  lres <- foreach(dataid = dataids, .packages = "spacemap") %do% { 
    ddir <- file.path(basedir, paste0("d", sprintf("%03d", dataid)))
    cvOut <- readRDS(file = paste0(ddir, "/cv_object_dataid_", dataid, ".rds"))
    tid <- cvOut$minIndices["rss"]                                 
    am <- avgMetrics(cvOut, testSetLen)
    log_rss <- log(am[,"rss"])
    stopifnot(which.min(am[,"rss"]) == tid)
    if (vote)   { 
      tf <- list.files(path = ddir, pattern = paste0("tuneid_", sprintf("%03d", tid)), full.names = TRUE)
      cvrds <- cvVoteRDS(mod_lambda = tf, tol = tol, major_thresh= 0.5, method = "space")
      saveRDS(cvrds, file = paste0(ddir, "/min_rss_cv_vote_dataid_", dataid, ".rds"))
    }
    data.frame(dataid = dataid, log_rss = min(log_rss), tuneid = tid)
  }
  library(data.table)
  res <- rbindlist(lres)
}

foldfile <- "/home/cconley/scratch-data/sim-spacemap/powlaw/pl-mod-01/datasets/n250/pl-mod-01-n250-train-test-partition.rds"
partition <- readRDS(file = foldfile)
testSetLen <- sapply(partition$testSetIds, length)
library(spacemap)

basedir <- "~/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/space/"
dir01 <- file.path(basedir, "01-cv-xy-yy-space/")
dir02 <- file.path(basedir, "02-cv-xy-yy-space/")

space_scores01 = report_cv_score(basedir = dir01, testSetLen = testSetLen)
space_scores02 = report_cv_score(basedir = dir02, testSetLen = testSetLen)

iter_tuneids <- cbind(i01 = space_scores01$tuneid, i02 = space_scores02$tuneid)
top_iter <- ifelse(space_scores02$log_rss > space_scores01$log_rss, "i01", "i02")
top_tuneid <- sapply(dataids, function(i) iter_tuneids[i,top_iter[i]])
#integer version for file construction
itop_iter <- ifelse(space_scores02$log_rss > space_scores01$log_rss, 1, 2)
cv_selected <- data.frame(dataid = dataids, top_iter = itop_iter, top_tuneid = top_tuneid)
saveRDS(cv_selected, file = file.path(basedir, "top_cv_selected.rds"))

#find best hold outs per data id
get_best_holdout_files <- function(i, cv_selected) {
  ddir <- file.path(basedir,
                    paste0(sprintf("%02d", cv_selected$top_iter[i]), "-cv-xy-yy-space"),
                    paste0("d", sprintf("%03d", cv_selected$dataid[i])))
  bestho <- list.files(path = ddir, pattern = paste0("tuneid_",
                                                     sprintf("%03d",cv_selected$top_tuneid[i])),
                       full.names = TRUE)
}


