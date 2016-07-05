

#identify which iteration has the min cv score
library(data.table)
min_score <- function(dt) dt[,list(score = min(cv_scggm)), by = dataid]
library(foreach)

iter_res <- foreach(i = 1:3) %do% { 
  dir <- paste0("/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/scggm/", sprintf("%02d", i), "-cv-xy-yy-scggm")
  file <- "cv_score_scggm_summary.rds"
  readRDS(file.path(dir, file))
}
names(iter_res) <- paste0("iter_", 1:3)

mcv <- foreach(res = iter_res, .combine = 'cbind') %do% { 
  res$min_cv_scores <- min_score(res$cv_scores)
  #make sure index matches up with ordered data id
  stopifnot(identical(res$min_cv_scores$dataid, 1:100))
  res$min_cv_scores$score 
}
top_iter <- apply(mcv, 1, which.min)

wmcv <- foreach(res = iter_res, .combine = 'cbind') %do% { 
  #make sure index matches up with ordered data id
  stopifnot(identical(res$top$dataid, 1:100))
  res$top$best_id
}
top_tuneid <- sapply(1:100, function(i) wmcv[i,top_iter[i]])

cv_selected <- data.frame(dataid = 1:100, top_iter, top_tuneid)

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

