

#identify which iteration has the min cv score
library(data.table)
min_score <- function(dt) dt[,list(score = min(cv_scggm)), by = dataid]
library(foreach)

#
did <- c(1:24, 26:29, 31:100)
iters <- 1:2
iter_res <- foreach(i = iters) %do% { 
  dir <- paste0("~/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/scggm/", sprintf("%02d", i), "-cv-xy-yy-scggm")
  file <- "cv_score_scggm_summary.rds"
  readRDS(file.path(dir, file))
}
names(iter_res) <- paste0("iter_", iters)

mcv <- foreach(res = iter_res, .combine = 'cbind') %do% { 
  res$min_cv_scores <- min_score(res$cv_scores)
  #make sure index matches up with ordered data id
  stopifnot(identical(res$min_cv_scores$dataid, did))
  res$min_cv_scores$score 
}
if(!is.matrix(mcv)) { 
 mcv <- matrix(mcv, ncol = 1)
}
top_iter <- apply(mcv, 1, which.min)

wmcv <- foreach(res = iter_res, .combine = 'cbind') %do% { 
  #make sure index matches up with ordered data id
  stopifnot(identical(res$top$dataid, did))
  res$top$best_id
}
if(!is.matrix(wmcv)) { 
 wmcv <- matrix(wmcv, ncol = 1)
}

top_tuneid <- sapply(seq_along(did), function(i) wmcv[i,top_iter[i]])

cv_selected <- data.frame(dataid = did, top_iter, top_tuneid)

#find best hold outs per data id
get_best_holdout_files <- function(i, cv_selected) { 
  ddir <- file.path(basedir,  
                    paste0(sprintf("%02d", cv_selected$top_iter[i]), "-cv-xy-yy-scggm"), 
                    paste0("d", sprintf("%03d", cv_selected$dataid[i])))
  bestho <- list.files(path = ddir, pattern = paste0("scggm_cv_tuneid_", 
                                                     sprintf("%03d",cv_selected$top_tuneid[i])), 
                       full.names = TRUE)
}

set_top_cv_vote_path <- function(i, cv_selected) { 
  ddir <- file.path(basedir,  
                    paste0(sprintf("%02d", cv_selected$top_iter[i]), "-cv-xy-yy-scggm"), 
                    paste0("d", sprintf("%03d", cv_selected$dataid[i])))
  file.path(ddir, "scggm_top_cv_vote_object.rds")
}

