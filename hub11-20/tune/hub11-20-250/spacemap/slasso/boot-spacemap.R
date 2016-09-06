#Obtain data set id
args <- commandArgs(TRUE)
dataid <- as.integer(args[1])
ncores <- as.integer(args[2])

#load common settings
source("~/repos/sim-spacemap/hub11-20/tune/hub11-20-250/common-tune-params.R")
#load input data
dat <- readRDS(paste0("~/scratch-data/sim-spacemap/hub11-20/2016/datasets/n250/hub11-20-250-dataid_", sprintf("%03d", dataid), ".rds"))

#setup parallelism
suppressPackageStartupMessages(library(doParallel))
cl <- makeCluster(ncores)
registerDoParallel(cl)

#result path 
res_path <- paste0("/home/cconley/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/boot-spacemap/d", sprintf("%03d", dataid))
cmd <- paste("mkdir -p", res_path)
system(cmd)

#set up tuning grid 
best_tune_path <- "~/scratch-data/sim-spacemap/hub11-20/2016/results/n250/spacemap/cv_selected_tune_params_per_dataid.rds"
best_tune <- readRDS(best_tune_path)
tune_param_cols <- which(!(names(best_tune) %in% "dataid"))
tune_per_dataid <- best_tune[match(dataid, best_tune$dataid),tune_param_cols]

#Number of bootstraps 
B <- 1000L
boot <- TRUE
seed <-  73311L
refitRidge <- 0.0
perc_boot <- 0.9

library(spacemap)
library(foreach)
library(doRNG)

#fit to bootstrap replicates
tictoc <- system.time({boot_cv_fits <- spacemap::reprodEdge(Y.m = dat$Y, X.m = dat$X, tune = tune_per_dataid, 
                              boot = boot, B = B, seed = seed, 
                              p0 = perc_boot, iter = iter, tol = tol, cd_iter = cd_iter, 
                              refitRidge = refitRidge)})
saveRDS(boot_cv_fits, file = file.path(res_path, paste0("ninety_perc_", "boot_replicates_dataid_", sprintf("%03d", dataid), ".rds")))
saveRDS(tictoc[3], file = file.path(res_path, paste0("ninety_perc_", "boot_time_dataid_", sprintf("%03d", dataid), ".rds")))
save.image(file = file.path(res_path, paste0("ninety_perc_", "boot_dataid_", sprintf("%03d", dataid), ".rda")))

#Bootstrap vote
boot_vote <- function(boot_cv_fits, vote_thresh = 0.5) {
  P <- nrow(boot_cv_fits[[1]][[1]]$Gamma)
  Q <- ncol(boot_cv_fits[[1]][[1]]$Gamma)
  library(spacemap)
  library(foreach)
  boot_spmap <- foreach(bspmap = boot_cv_fits[[1]], .combine = spacemap::addmods) %do% { 
    if (bspmap$convergence) { 
      bxy <- abs(as.matrix(bspmap$Gamma)) > tol
      byy <- abs(as.matrix(bspmap$ParCor)) > tol
      list(Gamma = bxy, ParCor = byy, 
           dfGamma = nonZeroWhole(bxy, tol), 
           dfParCor <- nonZeroUpper(byy, tol), nuniq = bspmap$convergence)
    } else { 
      list("Gamma" = matrix(0.0,P,Q), "ParCor" = matrix(0.0,Q,Q),  
           "dfGamma" = "No Convergence", "dfParCor" = "No Convergence", "nuniq" = bspmap$convergence)
    }
  }
  diag(boot_spmap$ParCor) <- 0
  
  #adjust for the number that actually converged
  effB <- sum(boot_spmap$nuniq)
  
  #vote based on effective bootstraps replicates. 
  bv <- list(ParCor = boot_spmap$ParCor > vote_thresh*effB,
             Gamma = boot_spmap$Gamma > vote_thresh*effB)
  
  list(bv = bv,
       bc = boot_spmap)
}

bvSpmap <- boot_vote(boot_cv_fits)
saveRDS(bvSpmap, file = file.path(res_path, paste0("ninety_perc_", "boot_vote_dataid_", sprintf("%03d", dataid), ".rds")))

stopCluster(cl)
