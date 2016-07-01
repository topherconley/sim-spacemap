#Description: Given a list of data directories, parse out the CV score predicted by scggm. 
#report the errors with the data IDs, the tune IDs, and report the top performing tune IDs.
#returns a data
report_cv_score <- function(ddirs, topn = 5) { 
  library(foreach)
  cve <- foreach(i = seq_along(ddirs)) %do% { 
    dir <- ddirs[i]
    dataid <- as.integer(substring(dir, first = nchar(dir) - 2, last = nchar(dir)))
    cve_meta_file <- list.files(path = dir, pattern = "cv_error_files_scggm_dataid_", full.names = T)
    cv_error_files <- readRDS(cve_meta_file)
    tuneid <- sapply(cv_error_files, function(cefs) {
      ti <- as.integer(regmatches(cefs, regexpr(pattern = "(?<=tuneid_)(.*?)(?=_fold)", 
                                                text = cefs, perl = TRUE)))
      uti <- unique(ti)
      stopifnot(length(uti) == 1)
      uti
    })
    cv_error_groups <- lapply(cv_error_files, function(cefs) sapply(cefs, readLines))
    
    mean_cv_error <- sapply(cv_error_groups, function(x) mean(as.numeric(x)))
    data.frame(cv_scggm = mean_cv_error, tuneid = tuneid, dataid = dataid)
  }
  library(data.table)
  cve.dt  <- data.table(rbindlist(cve))
  best_tuneid_by_data <- cve.dt[,list(best_id = which.min(cv_scggm)), by = dataid]
  best_tuneids_by_data <- cve.dt[,list(best_id = sort(cv_scggm, index.return = TRUE)$ix[1:topn]), by = dataid]
  list(top = best_tuneid_by_data, top_few = best_tuneids_by_data, cv_scores = cve.dt)
}

