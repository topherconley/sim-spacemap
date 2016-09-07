
setwd("~/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap")

top <- new.env()
load(envir = top, file = "spacemap_top_cv_voting.rda")
best_path <- sapply(top$best_holdout_files, function(x) dirname(x[1]))
best_path <- gsub(pattern = "/home/cconley", "~", best_path)
best_cv_object_file <- file.path(best_path, paste0("cv_object_dataid_", 1:100, ".rds"))

#get data id
get_tune <- function(f) { 
  tmp <- readRDS(f) 
  res <- tmp$tune
  ddid <- regmatches(x = f, m = regexpr(pattern = "d[0-9]{3}", text = f))
  res$dataid <- as.integer(sub("d", "", ddid))
  res
}

best_tune <- do.call(rbind,lapply(best_cv_object_file, get_tune))
saveRDS(best_tune, file = "cv_selected_tune_params_per_dataid.rds")

