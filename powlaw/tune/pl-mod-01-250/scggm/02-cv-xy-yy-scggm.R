#Obtain data set id
args <- commandArgs(TRUE)
dataid <- as.integer(args[1])
ncores <- as.integer(args[2])

#load common settings
source("~/repos/sim-spacemap/powlaw/tune/pl-mod-01-250/common-tune-params.R")
#load input data
dat <- readRDS(paste0("~/scratch-data/sim-spacemap/powlaw/pl-mod-01/datasets/n250/pl-mod-01-dataid_", sprintf("%03d", dataid), ".rds"))

#output base directory
data_part_dir <- "~/scratch-data/sim-spacemap/powlaw/pl-mod-01/datasets/n250/partition/"
xtrainfile <- vector(mode = "character", length = length(testSetIds))
xtestfile <- vector(mode = "character", length = length(testSetIds))
ytrainfile <- vector(mode = "character", length = length(testSetIds))
ytestfile <- vector(mode = "character", length = length(testSetIds))

#maintain double precision (1e-16)
dp <- function(x) t(apply(x, 1, function(x) sprintf("%.16f", x)))
for(i in seq_along(testSetIds)) { 
   xtrainfile[i] <- paste0(data_part_dir, "xtrain_dataid_", sprintf("%03d", dataid), "_fold_", sprintf("%02d", i), ".txt")
   xtestfile[i] <- paste0(data_part_dir, "xtest_dataid_", sprintf("%03d", dataid), "_fold_", sprintf("%02d", i), ".txt")
   ytrainfile[i] <- paste0(data_part_dir, "ytrain_dataid_", sprintf("%03d", dataid), "_fold_", sprintf("%02d", i), ".txt")
   ytestfile[i] <- paste0(data_part_dir, "ytest_dataid_", sprintf("%03d", dataid), "_fold_", sprintf("%02d", i), ".txt")
   #Already written to file from step 01
   #write.table(dp(dat$X[trainSetIds[[i]],]), file=xtrainfile[i], sep=",", row.names=FALSE, col.names=FALSE, quote = F)
   #write.table(dp(dat$X[testSetIds[[i]],]), file=xtestfile[i], sep=",", row.names=FALSE, col.names=FALSE, quote = F)
   #write.table(dp(dat$Y[trainSetIds[[i]],]), file=ytrainfile[i], sep=",", row.names=FALSE, col.names=FALSE, quote = F)
   #write.table(dp(dat$Y[testSetIds[[i]],]), file=ytestfile[i], sep=",", row.names=FALSE, col.names=FALSE, quote = F)
}


#setup parallelism
suppressPackageStartupMessages(library(doParallel))
cl <- makeCluster(ncores)
registerDoParallel(cl)

#result path 
res_path <- paste0("~/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/scggm/02-cv-xy-yy-scggm/d", sprintf("%03d", dataid))
cmd <- paste("mkdir -p", res_path)
system(cmd)

##########TUNING#######################

#prev tuning scggm
default_grid <- c(0.32, 0.16, 0.08, 0.04, 0.02, 0.01)
tsc_prev <- expand.grid(lam1 = default_grid,
                     lam2 = default_grid)
dir_prev <- "~/scratch-data/sim-spacemap/powlaw/pl-mod-01/results/n250/scggm/01-cv-xy-yy-scggm/"
res_prev <- readRDS(file  = file.path(dir_prev, "cv_score_scggm_summary.rds"))
#find the top tuning parameter for this specific dataset. 
top_tsc_prev_id <- res_prev$top$best_id[res_prev$top$dataid == dataid]
top_lam1 <- tsc_prev[top_tsc_prev_id,"lam1"]
top_lam2 <- tsc_prev[top_tsc_prev_id,"lam2"]
#create new grid
tsc <- expand.grid(lam1 = seq(top_lam1*.60, top_lam1*1.90, length  = 6),
                    lam2 = seq(top_lam2*.30, top_lam2*1.90, length  = 6))
saveRDS(tsc, file = file.path(res_path, paste0("scggm_tune_grid_dataid_", sprintf("%03d", dataid), ".rds"))) 
#################################
message("Running sCGGM (YY,XY) cross-validation")
#code for scggm is hosted here.
scggm_dir <- "~/repos/scggm/"
#tmp dir to write matlab to file
library(spacemap)
library(foreach)
library(Matrix)
pkgs <- "Matrix"
cve_scggm <- foreach(l = seq_len(nrow(tsc)), .packages = pkgs) %:% foreach(ff = 1:kcv, .combine = 'c') %dopar% { 
   xy_theta_file <- file.path(res_path, paste0("xy_theta", "_tuneid_", sprintf("%03d", l), "_foldid_", sprintf("%02d", ff), ".csv"))
   yy_theta_file <- file.path(res_path, paste0("yy_theta", "_tuneid_", sprintf("%03d", l), "_foldid_", sprintf("%02d", ff), ".csv"))
   error_file <- file.path(res_path, paste0("cv_error", "_tuneid_", sprintf("%03d", l), "_foldid_", sprintf("%02d", ff), ".csv"))
   matlabcmd <- paste0("matlab -nojvm -nodisplay -nosplash -singleCompThread -r \"cd '", scggm_dir, "'; ", "call_scggm_cv(",
                       "'", xtrainfile[ff], "'",",", "'", xtestfile[ff],"'", ",", 
                       "'", ytrainfile[ff], "'",",", "'", ytestfile[ff],"'", ",", 
                       tsc$lam1[l],",",tsc$lam2[l], ",", 
                       "'", xy_theta_file, "'",",", "'", yy_theta_file,"'", ",", 
                       "'", error_file, "'","); exit;\"")
   system(matlabcmd)
   ThetaXY = as.matrix(read.csv(xy_theta_file, header = FALSE))
   ThetaYY = as.matrix(read.csv(yy_theta_file, header = FALSE))
   system(paste("rm", xy_theta_file, yy_theta_file))
   out <- list(ThetaXY = Matrix(ThetaXY), ThetaYY = Matrix(ThetaYY), 
               tune_id = l, fold_id = ff)
   outfile <- file.path(res_path, paste0("scggm_cv_tuneid_", 
                        sprintf("%03d", l), "_foldid_", sprintf("%02d", ff), ".rds"))
   saveRDS(out, file = outfile)
   error_file
}
stopCluster(cl)

#######################################
saveRDS(cve_scggm, file = file.path(res_path, paste0("cv_error_files_scggm_dataid_", sprintf("%03d", dataid), ".rds")))
