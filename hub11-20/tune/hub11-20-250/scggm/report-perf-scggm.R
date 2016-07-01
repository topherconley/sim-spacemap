scggm_cv_perf <- function() { 
  comparison <- factor(c("(yy,xy)", "yy", "xy"), levels = c( "(yy,xy)", "yy", "xy"))
  #performance names to index
  pn <- list(power = c("power", "powerYY", "powerXY"),
             fdr = c("fdr", "fdrYY", "fdrXY"),
             mcc = c("mcc", "mccYY", "mccXY"),
             tp = c("tp", "tpYY", "tpXY"),
             fn = c("fn", "fnYY", "fnXY"),
             fp = c("fp", "fpYY", "fpXY"))
  #for each dataset
  perf <- foreach(i = seq_along(ddirs)) %do% {
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
  }:q
  :
  
  
}

scggm_perf <- function(sc, learn) {
  comparison <- factor(c("(yy,xy)", "yy", "xy"), levels = c( "(yy,xy)", "yy", "xy"))
  #performance names to index
  pn <- list(power = c("power", "powerYY", "powerXY"),
             fdr = c("fdr", "fdrYY", "fdrXY"),
             mcc = c("mcc", "mccYY", "mccXY"),
             tp = c("tp", "tpYY", "tpXY"),
             fn = c("fn", "fnYY", "fnXY"),
             fp = c("fp", "fpYY", "fpXY"))
  library(spacemap)
  library(foreach)
  ipf <- c("reportJointPerf")
  scp <- foreach(fit = sc$scggm, .export = ipf) %do% {
    fitParCor <- list(xy = as.matrix(fit$ThetaXY), yy = as.matrix(fit$ThetaYY))
    perf <- reportJointPerf(fitParCor, sc$trueParCor, tol = sc$tol, verbose = FALSE)
    data.frame(algorithm = "sCGGM",
               power = perf[pn$power], fdr =  perf[pn$fdr], mcc = perf[pn$mcc],
               tp = perf[pn$tp], fn = perf[pn$fn], fp = perf[pn$fp],
               comparison = comparison)
  }
  
  
  library(data.table)
  scperf <- as.data.frame(data.table::rbindlist(scp))
  scperf$lam1 <- rep(sc$tuneScggm$lam1, times = 3)
  scperf$lam2 <- rep(sc$tuneScggm$lam2, times = 3)
  scperf$learn <- learn
  scperf
}

